LL_probit_linear = function(par,y,z,x,w,verbose=1){
    if(length(par) != ncol(x)+ncol(z)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(z)]
    beta = par[ncol(z)+1:ncol(x)]
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)


    za = as.vector(z %*% alpha)
    xb = as.vector(x %*% beta)

    P = (w - za) / sigma
    Q = (2*y-1) * (xb  + rho*P) / sqrt(1-rho^2)

    LL = sum( pnorm(Q, log.p=TRUE) - log(sqrt(2*pi) * sigma) - P^2/2 )


    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}


Gradient_probit_linear = function(par,y,z,x,w,verbose=1,variance=FALSE){
    if(length(par) != ncol(x)+ncol(z)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(z)]
    beta = par[ncol(z)+1:ncol(x)]
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    za = as.vector(z %*% alpha)
    xb = as.vector(x %*% beta)

    P = (w - za) / sigma
    Q = (2*y-1) * (xb  + rho*P) / sqrt(1-rho^2)

    M = dnorm(Q) / pnorm(Q) * (2*y-1) / sqrt(1-rho^2)

    g_alpha = - ( (M * rho - P) * z ) / sigma
    g_beta = ( M * x )
    g_sigma = ( - rho * M * P - 1 + P^2 ) / sigma
    g_rho = ( M*P ) / (1-rho^2)

    g_log_sigma = g_sigma * sigma
    g_tau = g_rho * 2 * exp(tau) / (exp(tau)+1)^2


    dL = cbind(g_alpha, g_beta, g_log_sigma, g_tau)
    gradient = colSums(dL)
    names(gradient) = names(par)

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
        return (list(g=gradient, var=var, I=crossprod(dL), probit_prob = pnorm((2*y-1) * xb)))
    }
    return(gradient)
}

#' Recursive Probit-Linear Model
#' @description Estimate probit and linear models with bivariate normally distributed error terms. This command supports two models with opposite first and second stages.
#'
#' 1) Recursive Probit-Linear: the endogenous treatment effect model\cr
#' 2) Recursive Linear-Probit: the ivprobit model. The identification of this model requires an instrument.
#'
#' This command still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_probit Formula for the probit model
#' @param form_linear Formula for the linear model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param method Optimization algorithm. Default is BFGS
#' @param accu 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param verbose Level of output during estimation. Lowest is 0.
#' @return A list containing the results of the estimated model
#' @examples
#' library(MASS)
#' N = 2000
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
#' y1 = as.numeric(1 + x + z + e1 > 0)
#' y2 = 1 + x + z + y1 + e2
#'
#' est = probit_linear(y1~x+z, y2~x+z+y1)
#' est$estimates
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at SSRN: https://ssrn.com/abstract=3494856
probit_linear = function(form_probit, form_linear, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], verbose=0, accu=1e4){
    # 1.1 parse w~z from linear
    mf = model.frame(form_linear, data=data, na.action=NULL, drop.unused.levels=TRUE)
    w = model.response(mf, "numeric")
    z = model.matrix(attr(mf, "terms"), data=mf)
    # 1.2 parse y~x from probit
    mf2 = model.frame(form_probit, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf2, "numeric")
    x = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    est_linear = lm(form_linear, data=data)
    par_linear = coef(summary(est_linear))[,1]
    est_probit = glm(form_probit, data=data, family=binomial(link="probit"))
    par_probit = coef(summary(est_probit))[,1]
    names(par_linear) = paste0('linear.', names(par_linear))
    names(par_probit) = paste0('probit.', names(par_probit))
    par_linear[is.na(par_linear)] = 0
    par_probit[is.na(par_probit)] = 0
    par = c(par_probit, par_linear, log_sigma=0, tau=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par
    # print(par)

    # 2. Estimation
    resetIter()
    begin = Sys.time()

    # # use optim (hessian is forced to be symmetric)
    # res = optim(par=par, fn=LL_probit_linear, gr=Gradient_probit_linear, method="BFGS", control=list(factr=accu,fnscale=-1), y=y, z=z, x=x, w=w, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_probit_linear, grad=Gradient_probit_linear, start=par, method=method, y=y, z=z, x=x, w=w, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y)

    # 3. Compile results
    gvar = Gradient_probit_linear(res$par,y,z,x,w,verbose=verbose-1,variance=TRUE)
    if(any(is.na(gvar$var))){
        warning('NA in Hessian matrix, will reestimate with unif initialization!')
        return(probit_linear(form_linear, form_probit, data, par, init='unif', verbose=verbose, accu=accu))
    }

    res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(sigma='log_sigma', rho='tau'), trans_types=c('exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(est_linear) - logLik(est_probit) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$probit_prob = gvar$probit_prob
    res$iter = endogeneity.env$iter

    if(verbose>0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}

