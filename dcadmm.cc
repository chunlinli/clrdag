/*
ref: Li, Shen, Pan "Likelihood inference for a large directed acyclic graph"
author: Chunlin Li (li000007@umn.edu), Ziyue Zhu (zhux0502@umn.edu)
version 2019.02: DFS-cycle detection in algorithm is implemented for cpp function
*/

#include <math.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

/*
obj: calculate the objective 
X: data matrix
A: adjacency matrix, A must be a DAG!!!
W: weight matrix
D: index of hypothesized edges, no penalty is imposed
mu: sparsity par
*/
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
A: adjacency matrix, A must be a DAG!!!
W: weight matrix
D: index of hypothesized edges, no penalty is imposed
mu: sparsity par
*/
inline double objective(const mat &X, const mat &A, const umat &W, const umat &D, double mu)
{
    mat Y = X - X * A.t();
    double res = 0.5*accu(Y % Y) / X.n_rows;

    int p = (int) X.n_cols;
    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(D(i,j) == 0 && W(i,j) == 1)
            {
                res += mu * fabs(A(i,j));
            }
        }
    }

    return  res; 
} 


/*



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




*/
void force_DAG(mat &A, umat &W, double tau)
{
    int p = (int) A.n_cols;
    mat AABS = abs(A);
    for(int i = 0; i < p; ++i)
    {
        for(int j = 0; j < p; ++j)
        {
            if(AABS(i,j) < tau)
            {
                A(i,j) = 0;
                W(i,j) = 1;
                AABS(i,j) = INFINITY;
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
            A(end,start) = 0;
        }
    }

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

}


/*
DC_ADMM: solve the optimization problem by DC and ADMM algorithm
X: data matrix
A, Lambda: initial estimate, A must be a DAG!!!
D: index of hypothesized edges, no penalty is imposed
tau: threshold par in TLP
mu: sparsity par
rho: ADMM par
tol, rel_tol: absolute/relative tolerance
dc_max_iter, admm_max_iter: DC/ADMM max iteration 
trace_obj: boolean value, indicate whether to print the obj function after each DC iter
*/
// [[Rcpp::export]]
List DC_ADMM(mat X, mat A, mat Lambda, umat D, double tau, double mu, double rho, 
            double tol_abs, double tol_rel, int dc_max_iter, int admm_max_iter, int trace_obj)
{   /*
    passing references to R objects???
    */

    int p = (int) X.n_cols;

    // cache X^TX, inv(X^T X + rho*I)
    mat XTX = X.t() * X;
    cube INV_XTX = zeros(p-1,p-1,p);
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
        M(0,j) = 0;
        M(j,j) = 2.0 / p;
    }
    for(int i = 0; i < p; ++i)
    {
        M(i,0) = 1;
    }
    M = M / tau;

    // DC loop begins
    umat W(p,p);
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

    // initial objective
    double obj_curr = objective(X,A,W,D,mu); 

    for(int dc_iter = 0; dc_iter < dc_max_iter; ++dc_iter)
    {

        // ADMM loop begins
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
        force_DAG(A, W, tau);









        // DC stopping criteria
        double obj_next = objective(X,A,W,D,mu);
        if(obj_curr - obj_next < tol_abs)
        {
            break;
        }
        obj_curr = obj_next;

        // trace_obj
        if(trace_obj == 1)
        {
            Rcout << "DC iteration: " << dc_iter << ", objective value: " << obj_curr << endl;
        }

    } // DC loop ends

    // return
    return List::create(Named("A") = A, Named("Lambda") = Lambda);
}

