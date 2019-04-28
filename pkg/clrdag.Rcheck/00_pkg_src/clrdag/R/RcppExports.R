DC_ADMM <- function(X, A, Lambda, D, tau, mu, rho, tol_abs, tol_rel, dc_max_iter, admm_max_iter, auto_init, trace_obj) {
    .Call(`_clrdag_DC_ADMM`, X, A, Lambda, D, tau, mu, rho, tol_abs, tol_rel, dc_max_iter, admm_max_iter, auto_init, trace_obj)
}

