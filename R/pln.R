LL_pln = function(par, y, x, H=20, verbose=1){
    if(length(par) != ncol(x)+1) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x)]
    lambda = exp(par[length(par)])
    wa = as.vector(x %*% alpha)

    rule = gauss.quad(H, "hermite")
    Li = rep(0, length(y))
    for(k in 1:H){
        Li = Li + rule$weights[k] * dpois(y, exp(wa + lambda*sqrt(2)*rule$nodes[k]))
    }
    Li = pmax(Li, 1e-16)
    LL = sum(log(Li/sqrt(pi)))

    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

gradient_pln = function(par, y, x, H=20, verbose=1, variance=FALSE){
    if(length(par) != ncol(x)+1) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x)]
    lambda = exp(par[length(par)])
    wa = as.vector(x %*% alpha)

    rule = gauss.quad(H, "hermite")
    S = rep(0, length(y))
    dS = matrix(0, length(y), length(par))
    colnames(dS) = names(par)
    for(k in 1:H){
        w = rule$weights[k]
        u = sqrt(2) * rule$nodes[k]
        rate = exp(wa + lambda*u)

        S = S + w * dpois(y, rate)
        dS = dS + w * dpois(y, rate) * (y-rate) * cbind(x, u)
    }
    S = pmax(S, 1e-16)
    dL = 1/S * dS

    # accounting for transformation of parameters
    dL[, 'log_lambda'] = lambda * dL[, 'log_lambda']
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

predict.pln = function(model, data=NULL, type="response"){
    mf = model.frame(model$formula, data=data, na.action=NULL, drop.unused.levels=TRUE)
    x = model.matrix(attr(mf, "terms"), data=mf)

    par = model$par
    alpha = par[1:ncol(x)]
    lambda = exp(par[length(par)])
    wa = as.vector(x %*% alpha)

    exp(wa + 0.5*lambda^2)
}


#' Poisson Lognormal Model
#' @description Estimate a Poisson model with a log-normally distributed heterogeneity term, which is also referred to as the Poisson-Normal model.\cr\cr
#' \deqn{E[y_i|x_i,u_i]=exp(\boldsymbol{\alpha}'\mathbf{x_i}+\lambda u_i)}{E[y_i | x_i, u_i] = exp(\alpha' * x_i + \lambda * u_i)}
#' The estimates of this model are often similar to those of a negative binomial model.
#' @param form Formula
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param method Optimization algorithm.
#' @param init Initialization method
#' @param H Number of quadrature points
#' @param verbose A integer indicating how much output to display during the estimation process.
#' * <0 - No ouput
#' * 0 - Basic output (model estimates)
#' * 1 - Moderate output, basic ouput + parameter and likelihood in each iteration
#' * 2 - Extensive output, moderate output + gradient values on each call
#' @return A list containing the results of the estimated model, some of which are inherited from the return of maxLik
#' * estimates: Model estimates with 95% confidence intervals
#' * estimate or par: Point estimates
#' * variance_type: covariance matrix used to calculate standard errors. Either BHHH or Hessian.
#' * var: covariance matrix
#' * se: standard errors
#' * gradient: Gradient function at maximum
#' * hessian: Hessian matrix at maximum
#' * gtHg: \eqn{g'H^-1g}, where H^-1 is simply the covariance matrix. A value close to zero (e.g., <1e-3 or 1e-6) indicates good convergence.
#' * LL or maximum: Likelihood
#' * AIC: AIC
#' * BIC: BIC
#' * n_obs: Number of observations
#' * n_par: Number of parameters
#' * LR_stat: Likelihood ratio test statistic for the heterogeneity term \eqn{\lambda=0}
#' * LR_p: p-value of likelihood ratio test
#' * iterations: number of iterations taken to converge
#' * message: Message regarding convergence status.
#'
#' Note that the list inherits all the components in the output of maxLik. See the documentation of maxLik for more details.
#' @md
#' @examples
#' library(MASS)
#' N = 2000
#' set.seed(1)
#'
#' # Works well when the variance of the normal term is not overly large
#' # When the variance is very large, it tends to be underestimated
#' x = rbinom(N, 1, 0.5)
#' z = rnorm(N)
#' y = rpois(N, exp(-1 + x + z + 0.5 * rnorm(N)))
#' est = pln(y~x+z)
#' print(est$estimates, digits=3)
#' @export
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at https://doi.org/10.1287/isre.2022.1113
pln = function(form, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], H=20, verbose=0){
    # 1.1 parse y~x
    mf = model.frame(form, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)

    # 1.2 Initialize parameters
    psn = glm(form, data=data, family='poisson')
    par_psn = coef(summary(psn))[,1]
    par = c(par_psn, log_lambda=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par

    # 2. Estimation
    resetIter()
    begin = Sys.time()

    # optim
    # res = optim(par=par, fn=LL_pln, gr=NULL, method=method, control=list(fnscale=-1), y=y, x=x, H=H, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_pln, grad=gradient_pln, start=par, method=method, y=y, x=x, H=H, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y)

    # 3. Compile results
    res = getVarSE(res, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(lambda='log_lambda'), trans_types=c('exp'))
    res$LR_stat = 2 * (res$LL - logLik(psn))
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$psn = summary(psn)
    res$iter = endogeneity.env$iter
    res$formula = form

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
