LL_pln_probit = function(par,y1,y2,x1,x2,H=20,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par['log_lambda'])
    tau = par['tau']
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)

    rule = gauss.quad(H, "hermite")
    Li = rep(0, length(y2))
    for(k in 1:H){
        Phi = pnorm((2*y2-1) / sqrt(1-rho^2) * (xb + sqrt(2)*rho*rule$nodes[k]))
        Pm = dpois(y1, exp(wa + sqrt(2)*lambda*rule$nodes[k]))
        Li = Li + rule$weights[k] * Pm * Phi
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

gradient_pln_probit = function(par,y1,y2,x1,x2,H=20,verbose=1,variance=FALSE){
    if(length(par) != ncol(x1)+ncol(x2)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par['log_lambda'])
    tau = par['tau']
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)
    y_rho = (2*y2-1) / sqrt(1-rho^2)

    rule = gauss.quad(H, "hermite")
    S = rep(0, length(y2))
    dS = matrix(0, length(y2), length(par))
    colnames(dS) = names(par)
    for(k in 1:H){
        w = rule$weights[k]
        u = sqrt(2)*rule$nodes[k]
        q = y_rho * (xb + rho*u)
        rate = exp(wa + lambda*u)

        Phi = pnorm(q)
        phi = dnorm(q)
        Pm = dpois(y1, rate)
        S = S + w * Pm * Phi

        ix = c(1:ncol(x1), length(par)-1)
        dS[, ix] = dS[, ix] + w*Phi*Pm*(y1-rate)*cbind(x1, u)
        dS[, ncol(x1)+1:ncol(x2)] = dS[, ncol(x1)+1:ncol(x2)] + w*Pm*phi*y_rho*x2
        dS[, 'tau'] = dS[, 'tau'] + w*(y_rho*u + rho*q/(1-rho^2))*Pm*phi
    }
    S = pmax(S, 1e-16)
    dL = 1/S * dS

    # accounting for transformation of parameters
    dL[, 'log_lambda'] = lambda * dL[, 'log_lambda']
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

#' Recursive PLN-Probit Model
#' @description Estimate a Poisson Lognormal model (first-stage) and a probit model (second-stage) whose error terms are bivariate normally distributed. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_pln Formula for the first-stage Poisson lognormal model
#' @param form_probit Formula for the second-stage probit model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param H Number of quadrature points
#' @param method Optimization algorithm. Without gradient, NM is much faster than BFGS
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
#' y1 = rpois(N, exp(-1 + x + z + e1))
#' y2 = as.numeric(1 + x + z + log(1+y1) + e2 > 0)
#' est = pln_probit(y1~x+z, y2~x+z+log(1+y1))
#' est$estimates
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at SSRN: https://ssrn.com/abstract=3494856
pln_probit = function(form_pln, form_probit, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], H=20, verbose=0, accu=1e4){
    # 1.1 parse y1~x1
    mf1 = model.frame(form_pln, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y1 = model.response(mf1, "numeric")
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y2~x2
    mf2 = model.frame(form_probit, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y2 = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    begin = Sys.time()
    pln = pln(form_pln, data=data, verbose=-1, method=method)
    probit = glm(form_probit, data=data, family=binomial(link='probit'))
    par_pln = pln$estimates[,1]
    par_probit = coef(summary(probit))[,1]
    names(par_pln) = paste0('pln.', names(par_pln))
    par = c(par_pln[-length(par_pln)], par_probit, log_lambda=as.numeric(log(par_pln[length(par_pln)])), tau=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par
    # print(par)

    # 2. Estimation
    resetIter()

    # optim
    # res = optim(par=par, fn=LL_pln_probit, gr=NULL, method=method, control=list(factr=accu,fnscale=-1), y1=y1, y2=y2, x1=x1, x2=x2, H=H, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_pln_probit, grad=gradient_pln_probit, start=par, method=method, y1=y1, y2=y2, x1=x1, x2=x2, H=H, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y1)

    # 3. Compile results
    res = getVarSE(res, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(lambda='log_lambda', rho='tau'), trans_types=c('exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - pln$LL - logLik(probit) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$pln = pln
    res$probit = summary(probit)
    res$iter = endogeneity.env$iter

    if(verbose>0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
