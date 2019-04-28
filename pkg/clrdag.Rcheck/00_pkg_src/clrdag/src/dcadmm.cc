/*
ref: Li, Shen, Pan "Likelihood inference for a large directed acyclic graph"
author: Chunlin Li (li000007@umn.edu)
version 0.99.99: work in parallel, update auto initialization function
*/

#include <omp.h>
#include <Rmath.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

inline double SIGN(double x)
{
    return (x < 0 ? (-1) : 1);
} 

inline double MAX(double x, double y)
{
    return (x >= y ? x : y);
}

/*
obj: calculate the objective 
X: data matrix
A: adjacency matrix
W: weight matrix
D: index of hypothesized edges, no penalty is imposed
mu: sparsity par
*/
//inline double objective(const mat &X, const mat &A, const umat &W, const umat &D, double mu)
inline double objective(const mat &X, const mat &A)
{
    mat Y = X - X * A.t();
    return  0.5*accu(Y % Y) / X.n_rows;
} 


/*
DFS_cycle: depth-first search for cycle detection
           return 1 if the edge (start -> end) added to A forms a cycle 
           return 0 otherwise
A: adjacency matrix
start, end: indices of starting nodes
*/
int DFS_cycle(const mat &A, int start, int end)
{
    if(start == end)
    {
        return 1;
    }

    int p = (int) A.n_cols;
    int bottom = 0;
    int top = 0;
    int S_len = 1;
    int cycle = 0;

    uvec AN(p);
    AN.fill(0);
    uvec S(p);
    S.fill(0);

    AN(start) = 1;
    S(0) = start;

    while(S_len > 0)
    {
        int i = S(bottom);
        S_len--;
        bottom++;

        for(int j = 0; j < p; ++j)
        {
            if(A(i,j) != 0)
            {
                if(j == end)
                {
                    cycle = 1;
                    break;
                }
                else
                {
                    if(AN(j) == 0)
                    {
                        top++;
                        S(top) = j;
                        S_len++;
                        AN(j) = 1;
                    }
                }
            }
        }

        if(cycle == 1)
        {
            break;
        }
    }

    return cycle;
}


/*
force_DAG: truncate the minimum strength edges to satisfy the DAG constraint
XTX: p by p matrix, t(X)*X
A: adjacency matrix
W: weight matrix
D: p by p integer, index of hypothesized edges, no penalty is imposed
tau: real, threshold par in TLP > 0
*/
void force_DAG(const mat &XTX, mat &A, umat &W, umat &D, double tau)
{
    int p = (int) A.n_cols;
    mat AABS = abs(A);

    #pragma omp parallel for
    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(AABS(i,j) < tau && D(i,j) == 0)
            {
                A(i,j) = 0.0;
                AABS(i,j) = INFINITY;
                W(i,j) = 1;
            }
            else
            {
                W(i,j) = 0;
            }
        }
    }

    uvec idx = sort_index(AABS);
    for(int i = 0; i < p*p; ++i)
    {
        if(isinf(AABS(idx(i))))
        {
            break;
        }

        int start = idx(i) / p;
        int end = idx(i) % p;
        if(DFS_cycle(A,start,end))
        {
            A(end,start) = 0.0;
            W(end,start) = 1;
        }
    }
/*
    #pragma omp parallel for
    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(fabs(A(i,j)) < tau)
            {
                W(i,j) = 1;
            }
            else 
            {
                W(i,j) = 0;
            }
        }
    }
*/
    // use OLS estimate 
    #pragma omp parallel for
    for(int i = 0; i < p; ++i)
    {
        uvec idx = find(A.row(i) != 0.0);
        int p_TEMP = (int) idx.size();

        if(p_TEMP > 0)
        {
            vec A_ROW_TEMP = XTX.col(i);
            A_ROW_TEMP = solve(XTX(idx,idx)+1e-3*eye(p_TEMP,p_TEMP), A_ROW_TEMP(idx));

            for(int j = 0; j < p_TEMP; ++j)
            {
                A(i,idx(j)) = A_ROW_TEMP(j);
                if(fabs(A(i,idx(j))) < tau && D(i,idx(j)) == 0)
                {
                    A(i,idx(j)) = 0.0;
                    W(i,idx(j)) = 1;
                }
            }
        }
    }

}


/*
DC_ADMM: solve the optimization problem by DC and ADMM algorithm
X: n by p real, data matrix 
A, Lambda: p by p real, initial estimate, A must be a DAG!!! Lambda must be compatible with A!!!
D: p by p integer, index of hypothesized edges, no penalty is imposed
tau: real, threshold par in TLP > 0
mu: real, sparsity par > 0
rho: real, ADMM par > 0
tol, rel_tol: real, absolute/relative tolerance > 0
dc_max_iter, admm_max_iter: integer, DC/ADMM max iteration 
auto_init: integer, indicate whether to use automatically initial estimate (A, Lambda), use with caution!!!
trace_obj: integer, indicate whether to print the obj function after each DC iter
*/
// [[Rcpp::export]]
List DC_ADMM(mat X, mat A, mat Lambda, umat D, double tau, double mu, double rho, 
            double tol_abs, double tol_rel, int dc_max_iter, int admm_max_iter, int auto_init, int trace_obj)
{   /*
    passing references to R objects???
    */

    int p = (int) X.n_cols;

    // cache X^TX, inv(X^T X + rho*I)
    mat XTX = X.t() * X;
    cube INV_XTX = zeros(p-1,p-1,p);

    #pragma omp parallel for
    for(int j = 0; j < p; ++j)
    {
        mat MAT_TEMP = XTX;
        MAT_TEMP.shed_row(j);
        MAT_TEMP.shed_col(j);
        INV_XTX.slice(j) = inv(MAT_TEMP + rho*eye(p-1, p-1));
    }

    // cache matrix M
    mat M = ones(p,p) / p;
    for(int j = 0; j < p; ++j)
    {
        M(0,j) = 0.0;
        M(j,j) = 2.0 / p;
    }
    for(int i = 0; i < p; ++i)
    {
        M(i,0) = 1;
    }
    M = M / tau;
    
    // auto initialization of A, Lambda
    if (auto_init == 1) 
    {
        #pragma omp parallel for
        for(int j = 0; j < p; ++j)
        {
            A.col(j) = XTX.col(j)/XTX(j,j);
            A(j,j) = 0.0;

            for(int i = 0; i < p; ++i)
            {
                if(fabs(A(i,j)) < 2*tau && D(i,j) == 0)
                {
                    A(i,j) = 0.0;
                }
            }
        }
    }

    // DC loop begins
    umat W(p,p);
    #pragma omp parallel for
    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(fabs(A(i,j)) < tau && D(i,j) == 0)
            {
                W(i,j) = 1;
            }
            else
            {
                W(i,j) = 0;
            }
        }
    }

    force_DAG(XTX, A, W, D, tau);
    mat A_curr = A;
    mat Lambda_curr = Lambda;

    // initial objective
    double obj_curr = objective(X,A_curr);

    for(int dc_iter = 0; dc_iter < dc_max_iter; ++dc_iter)
    {
        // trace_obj
        if(trace_obj == 1)
        {
            Rcout << "DC iteration: " << dc_iter << ", objective value: " << obj_curr << endl;
        }

        // ADMM loop begins: for details consult the reference or vignette
        mat B = A;
        cube Xi = zeros(p,p,p);
        cube Y = zeros(p,p,p);
        mat U = zeros(p,p);

        double tol_pri = 0;
        double tol_dual = 0;
        mat R_pri = zeros(p,p);
        mat R_dual = zeros(p,p);

        for(int admm_iter = 0; admm_iter < admm_max_iter; ++admm_iter)
        {
            // A direction
            #pragma omp parallel for
            for(int j = 0; j < p; ++j)
            {

                mat MAT_TEMP = B;
                MAT_TEMP.shed_col(j);
                vec VEC_TEMP = trans(MAT_TEMP.row(j));

                MAT_TEMP = U;
                MAT_TEMP.shed_col(j);
                VEC_TEMP -= trans(MAT_TEMP.row(j));

                MAT_TEMP = XTX;
                MAT_TEMP.shed_row(j);
                VEC_TEMP = INV_XTX.slice(j)*(VEC_TEMP*rho + MAT_TEMP.col(j));

                for(int i = 0, i_VEC_TEMP = 0; i < p; ++i)
                {
                    if(i != j)
                    {
                        A(j,i) = VEC_TEMP(i_VEC_TEMP);
                        i_VEC_TEMP++;
                    }
                }
            }

            // B direction
            R_dual = B;
            #pragma omp parallel for
            for(int i = 0; i < p; ++i)
            {
                for(int j = 0; j < p; ++j)
                {
                    if(W(i,j) == 0) 
                    {
                        B(i,j) = A(i,j) + U(i,j);
                    }
                    else
                    {
                        B(i,j) = 0;
                        if(i != j) 
                        {
                            for(int k = 0; k < p; ++k)
                            {
                                B(i,j) += Lambda(i,k) - Lambda(j,k) - Xi(i,j,k) - Y(i,j,k);
                            }
                            B(i,j) += tau * (p-1); 
                            B(i,j) = (fabs(A(i,j) + U(i,j)) + B(i,j))/(p+1);
                            if(D(i,j) == 0)
                            {
                                B(i,j) -= mu/(rho*(p+1));
                            }
                            B(i,j) = SIGN(A(i,j) + U(i,j))*MAX(0,B(i,j));
                        }
                    }
                }
            }
            R_dual = rho*(B - R_dual);

            // Lambda direction
            mat VV(p,p);
            VV.fill(tau);
            for(int j = 0; j < p; ++j)
            {
                VV(0,j) = 1;
            }

            #pragma omp parallel for
            for(int i = 1; i < p; ++i)
            {
                for(int k = 0; k < p; ++k)
                {
                    for(int j = 0; j < p; ++j)
                    {
                        if(j != i)
                        {
                            VV(i,k) += Xi(i,j,k) + Y(i,j,k) - Xi(j,i,k) - Y(j,i,k);
                            if(W(i,j)==1)
                            {
                                VV(i,k) += fabs(B(i,j));
                            }
                            else
                            {
                                VV(i,k) += tau;
                            }
                            if(W(j,i)==1)
                            {
                                VV(i,k) -= fabs(B(j,i));
                            }
                            else
                            {
                                VV(i,k) -= tau;
                            }
                        }
                    }
                    if(i != k)
                    {
                        VV(i,k) *= 0.5;
                    }
                    else
                    {
                        VV(i,k) = 0.5*(VV(i,k) - p*tau);
                    }
                }
            }
            Lambda = M * VV;

            // Xi and Y directions
            #pragma omp parallel for
            for(int i = 0; i < p; ++i)
            {
                for(int j = 0; j < p; ++j)
                {
                    if(i != j)
                    {
                        for(int k = 0; k < p; ++k)
                        {
                            if(W(i,j) == 1)
                            {
                                Xi(i,j,k) = tau*(Lambda(i,k) - Lambda(j,k)) - fabs(B(i,j)) - Y(i,j,k);
                            }
                            else
                            {
                                Xi(i,j,k) = tau*(Lambda(i,k) - Lambda(j,k) - 1) - Y(i,j,k);
                            }
                            if(j != k)
                            {
                                Xi(i,j,k) += tau;
                            }
                            Xi(i,j,k) = MAX(0,Xi(i,j,k));

                            if(W(i,j) == 1)
                            {
                                Y(i,j,k) += fabs(B(i,j)) + Xi(i,j,k) - tau*Lambda(i,k) + tau*Lambda(j,k);
                            }
                            else
                            {
                                Y(i,j,k) += Xi(i,j,k) - tau*Lambda(i,k) + tau*Lambda(j,k) + tau;
                            }
                            if(j != k)
                            {
                                Y(i,j,k) -= tau;
                            }
                        }
                    }
                }
            }

            // U direction
            R_pri = A - B;
            U += R_pri;

            // ADMM stopping criteria
            tol_pri = tol_abs + tol_rel*((A % A + B % B).max());
            tol_dual = tol_abs + tol_rel*((B % B).max());
            // may need to change the stopping criteria
            if((R_pri % R_pri).max() < tol_pri*tol_pri && (R_dual % R_dual).max() < tol_dual*tol_dual)
            {
                break;
            }

        } // ADMM loop ends

        // enforce DAG constraints. refer to step 3 of algorithm 1 in the paper
        force_DAG(XTX, A, W, D, tau);

        // DC stopping criteria
        double obj = objective(X,A);
        if(obj_curr - obj < tol_abs)
        {
            break;
        }

        obj_curr = obj;
        A_curr = A;
        Lambda_curr = Lambda;

    } // DC loop ends

    // return
    return List::create(Named("A") = A_curr, Named("Lambda") = Lambda_curr);
}


/*
test_DF: estimate degrees of freedom of testing
A: adjacency matrix
D: p by p integer, index of hypothesized edges, no penalty is imposed
*/
// // [[Rcpp::export]]
/*
int test_DF(mat A, umat D) 
{
    int p = (int) A.n_cols;
    int df = 0;

    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(D(i,j) == 1)
                df++;
        }
    }

    return df;
}
*/
