LL_bilinear = function(par,y1,y2,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)

    Li = dnorm(y1, mean=wa, sd=lambda, log=TRUE) + dnorm(y2, mean=xb+sigma*rho/lambda*(y1-wa), sd=sigma*sqrt(1-rho^2), log=TRUE)
    LL = sum(Li)

    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}


Gradient_bilinear = function(par,y1,y2,x1,x2,verbose=1,variance=FALSE){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(x1 %*% alpha)
    xb = as.vector(x2 %*% beta)
    x3 = sigma*rho/lambda*(y1-wa)

    mu1 = y1-wa
    mu2 = y2-xb-x3
    sigma22 = sigma^2 * (1-rho^2)

    # rowwise multiplication on matrix x1
    dL_alpha = mu1/lambda^2 * x1 - (sigma*rho/lambda) * mu2/sigma22 * x1
    dL_beta = mu2/sigma22 * x2
    dL_lam = -1/lambda + mu1^2/lambda^3 - mu2*x3/(lambda*sigma22)
    dL_sigma = -1/sigma + (mu2*x3+mu2^2)/(sigma*sigma22)
    dL_rho = rho/(1-rho^2) + (mu2*(sigma*mu1/lambda)*(1-rho^2) - rho*mu2^2)/((1-rho^2)*sigma22)

    dL_loglam = dL_lam * lambda
    dL_logsigma = dL_sigma * sigma
    dL_tau = dL_rho * (2 * exp(tau) / (exp(tau)+1)^2)

    dL = cbind(dL_alpha, dL_beta, dL_loglam, dL_logsigma, dL_tau)
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
        return (list(g=gradient, var=var, I=crossprod(dL)))
    }
    return(gradient)
}

#' Recusrive Bivariate Linear Model
#' @description Estimate two linear models with bivariate normally distributed error terms.\cr\cr
#' First stage (Linear):
#' \deqn{m_i=\boldsymbol{\alpha}'\mathbf{w_i}+\lambda u_i}{m_i = \alpha' * w_i + \lambda * u_i}
#' Second stage (Linear):
#' \deqn{y_i = \boldsymbol{\beta}'\mathbf{x_i} + {\gamma}m_i + \sigma v_i}{y_i = \beta' * x_i + \gamma * m_i + \sigma * v_i}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' The identification of this model requires an instrumental variable that appears in w but not x. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form1 Formula for the first linear model
#' @param form2 Formula for the second linear model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param method Optimization algorithm. Default is BFGS
#' @param verbose A integer indicating how much output to display during the estimation process.
#' * <0 - No ouput
#' * 0 - Basic output (model estimates)
#' * 1 - Moderate output, basic ouput + parameter and likelihood in each iteration
#' * 2 - Extensive output, moderate output + gradient values on each call
#' @return A list containing the results of the estimated model, some of which are inherited from the return of maxLik
#' * estimates: Model estimates with 95% confidence intervals. Prefix "1" means first stage variables.
#' * estimate or par: Point estimates
#' * variance_type: covariance matrix used to calculate standard errors. Either BHHH or Hessian.
#' * var: covariance matrix
#' * se: standard errors
#' * var_bhhh: BHHH covariance matrix, inverse of the outer product of gradient at the maximum
#' * se_bhhh: BHHH standard errors
#' * gradient: Gradient function at maximum
#' * hessian: Hessian matrix at maximum
#' * gtHg: \eqn{g'H^-1g}, where H^-1 is simply the covariance matrix. A value close to zero (e.g., <1e-3 or 1e-6) indicates good convergence.
#' * LL or maximum: Likelihood
#' * AIC: AIC
#' * BIC: BIC
#' * n_obs: Number of observations
#' * n_par: Number of parameters
#' * LR_stat: Likelihood ratio test statistic for \eqn{\rho=0}
#' * LR_p: p-value of likelihood ratio test
#' * iterations: number of iterations taken to converge
#' * message: Message regarding convergence status.
#'
#' Note that the list inherits all the components in the output of maxLik. See the documentation of maxLik for more details.
#' @md
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
#' m = -1 + x + z + e1
#' y = -1 + x + m + e2
#'
#' est = bilinear(m~x+z, y~x+m)
#' print(est$estimates, digits=3)
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at https://doi.org/10.1287/isre.2022.1113
bilinear = function(form1, form2, data=NULL, par=NULL, method='BFGS', verbose=0){
    # 1.1 parse y1~x1
    mf1 = model.frame(form1, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y1 = model.response(mf1, "numeric")
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y2~x2
    mf2 = model.frame(form2, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y2 = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    ols1 = lm(form1, data=data)
    ols2 = lm(form2, data=data)
    par_ols1 = coef(summary(ols1))[,1]
    par_ols2 = coef(summary(ols2))[,1]
    names(par_ols1) = paste0('1', names(par_ols1))
    par = c(par_ols1, par_ols2, log_lambda=log(sigma(ols1)), log_sigma=log(sigma(ols2)), tau=0)
    # print(par)

    # 2. Estimation
    resetIter()
    begin = Sys.time()

    # res = optim(par=par, fn=LL_bilinear, gr=Gradient_bilinear, method="BFGS", control=list(fnscale=-1), y1=y1, y2=y2, x1=x1, x2=x2, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_bilinear, grad=Gradient_bilinear, start=par, method=method, y1=y1, y2=y2, x1=x1, x2=x2, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y1)

    # 3. Compile results
    gvar = Gradient_bilinear(res$par,y1,y2,x1,x2,verbose=verbose-1,variance=TRUE)
    res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(lambda='log_lambda', sigma='log_sigma', rho='tau'), trans_types=c('exp', 'exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(ols1) - logLik(ols2) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$iter = endogeneity.env$iter

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
