LL_pln_linear = function(par,y1,y2,x1,x2,H=20,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)
    q = y2 - xb

    rule = gauss.quad(H, "hermite")
    S = rep(0, length(y2))
    for(k in 1:H){
        u = sqrt(2*(1-rho^2)) * rule$nodes[k] + rho/sigma*q
        w = rule$weights[k]
        S = S + w * dpois(y1, exp(wa + lambda*u))
    }
    S = pmax(S, 1e-16)
    LL = sum( log(S) - log(sqrt(2)*pi*sigma) - q^2/(2*sigma^2) )

    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

gradient_pln_linear = function(par,y1,y2,x1,x2,H=20,verbose=1,variance=FALSE){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)
    q = y2 - xb

    rule = gauss.quad(H, "hermite")
    S = rep(0, length(y2))
    dS = matrix(0, length(y2), length(par))
    colnames(dS) = names(par)
    for(k in 1:H){
        w = rule$weights[k]
        r = rule$nodes[k]

        u = sqrt(2*(1-rho^2)) * r + rho/sigma*q
        rate = exp(wa + lambda*u)
        pr_y1 = dpois(y1, rate)

        S = S + w * pr_y1

        drate = cbind(x1, -lambda*rho/sigma*x2, u, -lambda*rho*q/sigma^2,
                      lambda*(-2*rho*r/sqrt(2*(1-rho^2)) + q/sigma) )
        dS = dS + w * pr_y1 * (y1-rate) * drate
    }
    S = pmax(S, 1e-16)
    dL = 1/S * dS

    # adding additional terms for beta and sigma
    dL[, ncol(x1)+1:ncol(x2)] = dL[, ncol(x1)+1:ncol(x2)] + q/sigma^2*x2
    dL[, 'log_sigma'] = dL[, 'log_sigma'] - 1/sigma + q^2/sigma^3

    # accounting for transformation of parameters
    dL[, 'log_lambda'] = lambda * dL[, 'log_lambda']
    dL[, 'log_sigma'] = sigma * dL[, 'log_sigma']
    dL[, 'tau'] = (2 * exp(tau) / (exp(tau)+1)^2) * dL[, 'tau']

    gradient = colSums(dL)

    if(verbose>=2){
        cat("----Gradient:\n")
        print(gradient,digits=3)
    }
    if(any(is.na(gradient) | !is.finite(gradient))) gradient = rep(NA, length(gradient))
    if(variance){
        var = tryCatch( solve(crossprod(dL)), error = function(e){
            cat('BHHH cross-product not invertible: ', e$message, '\n')
            diag(length(par)) * NA
        } )
        return (list(g=gradient, var=var, I=crossprod(dL)))
    }
    return(gradient)
}

#' Recursive PLN-Linear Model
#' @description Estimate a Poisson Lognormal model (first-stage) and a linear model (second-stage) with bivariate normally distributed error terms. This command still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_pln Formula for the first-stage Poisson lognormal model
#' @param form_linear Formula for the second-stage linear model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param method Optimization algorithm.
#' @param init Initialization method
#' @param H Number of quadrature points
#' @param accu 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param verbose Level of output during estimation. Lowest is 0.
#' @return A list containing the results of the estimated model
#' @examples
#' library(MASS)
#' N = 1000
#' rho = -0.5
#' set.seed(1)
#'
#' x = rbinom(N, 1, 0.5)
#' z = rnorm(N)
#'
#' e = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
#' e1 = e[,1]
#' e2 = e[,2]
#'
#' y1 = rpois(N, exp(1 + x + z + e1))
#' y2 = 1 + x + y1 + e2
#'
#' est = pln_linear(y1~x+z, y2~x+y1)
#' est$estimates
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at SSRN: https://ssrn.com/abstract=3494856
pln_linear = function(form_pln, form_linear, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], H=20, verbose=0, accu=1e4){
    # 1.1 parse y1~x1
    mf1 = model.frame(form_pln, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y1 = model.response(mf1, "numeric")
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y2~x2
    mf2 = model.frame(form_linear, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y2 = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    begin = Sys.time()
    # PoissonLogNormal may take some time
    pln = pln(form_pln, data=data, verbose=-1, method=method)

    par_pln = pln$estimates[,1]
    names(par_pln) = paste0('pln.', names(par_pln))

    ols = lm(form_linear, data=data)
    par_ols = coef(summary(ols))[,1]

    par = c(par_pln[-length(par_pln)], par_ols, log_lambda=as.numeric(log(par_pln[length(par_pln)])), log_sigma=log(sigma(ols)), tau=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par

    # 2. Estimation
    resetIter()

    # optim
    # res = optim(par=par, fn=LL_pln_linear, gr=gradient_pln_linear, method=method, control=list(factr=accu,fnscale=-1), y1=y1, y2=y2, x1=x1, x2=x2, H=H, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_pln_linear, grad=gradient_pln_linear, start=par, method=method, y1=y1, y2=y2, x1=x1, x2=x2, H=H, verbose=verbose)
    # res$no_grad_est = maxLik(LL_pln_linear, start=par, method=method, y1=y1, y2=y2, x1=x1, x2=x2, H=H, verbose=verbose)

    res$par = res$estimate
    res$n_obs = length(y1)

    # 3. Compile results
    res = getVarSE(res, verbose=verbose)
    # gvar = gradient_pln_linear(res$par,y1,y2,x1,x2,H=H,verbose=0,variance=TRUE)
    # res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(lambda='log_lambda', sigma='log_sigma', rho='tau'), trans_types=c('exp', 'exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - pln$LL - logLik(ols) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$pln = pln
    res$ols = summary(ols)
    res$iter = endogeneity.env$iter

    if(verbose>0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
