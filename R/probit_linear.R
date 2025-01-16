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
#' @description Estimate probit and linear models with bivariate normally distributed error terms.\cr\cr
#' First stage (Probit):
#' \deqn{m_i=1(\boldsymbol{\alpha}'\mathbf{w_i}+u_i>0)}{m_i = 1(\alpha' * w_i + u_i > 0)}
#' Second stage (Linear):
#' \deqn{y_i = \boldsymbol{\beta}'\mathbf{x_i} + {\gamma}m_i + \sigma v_i}{y_i = \beta' * x_i + \gamma * m_i + \sigma * v_i}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' w and x can be the same set of variables. Identification can be weak if w are not good predictors of m. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_probit Formula for the probit model
#' @param form_linear Formula for the linear model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param method Optimization algorithm. Default is BFGS
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
#' m = as.numeric(1 + x + z + e1 > 0)
#' y = 1 + x + z + m + e2
#'
#' est = probit_linear(m~x+z, y~x+z+m)
#' print(est$estimates, digits=3)
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2023) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research, 34(1):67-84. Available at https://doi.org/10.1287/isre.2022.1113
probit_linear = function(form_probit, form_linear, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], verbose=0){
    # Note: Be aware that y~x represents the first stage probit model
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
    par = c(par_linear, par_probit, log_sigma=0, tau=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par
    # print(par)

    # 2. Estimation
    resetIter()
    begin = Sys.time()

    # # use optim (hessian is forced to be symmetric)
    # res = optim(par=par, fn=LL_probit_linear, gr=Gradient_probit_linear, method="BFGS", control=list(fnscale=-1), y=y, z=z, x=x, w=w, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_probit_linear, grad=Gradient_probit_linear, start=par, method=method, y=y, z=z, x=x, w=w, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y)

    # 3. Compile results
    gvar = Gradient_probit_linear(res$par,y,z,x,w,verbose=verbose-1,variance=TRUE)
    if(any(is.na(gvar$var))){
        warning('NA in Hessian matrix, will reestimate with unif initialization!')
        return(probit_linear(form_probit, form_linear, data, par, init='unif', verbose=verbose))
    }

    res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(sigma='log_sigma', rho='tau'), trans_types=c('exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(est_linear) - logLik(est_probit) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$probit_prob = gvar$probit_prob
    res$iter = endogeneity.env$iter

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}

#' Recursive Linear-Probit Model
#' @description Estimate linear and probit models with bivariate normally distributed error terms.\cr\cr
#' First stage (Linear):
#' \deqn{m_i=\boldsymbol{\alpha}'\mathbf{w_i}+\sigma u_i}{m_i = \alpha' * w_i + \sigma*u_i}
#' Second stage (Probit):
#' \deqn{y_i = 1(\boldsymbol{\beta}'\mathbf{x_i} + {\gamma}m_i + v_i>0)}{y_i = 1(\beta' * x_i + \gamma * m_i + v_i > 0)}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' The identification of this model requires an instrumental variable that appears in w but not x. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_linear Formula for the linear model
#' @param form_probit Formula for the probit model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param method Optimization algorithm. Default is BFGS
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
#' m = 1 + x + z + e1
#' y = as.numeric(1 + x + m + e2 > 0)
#'
#' est = linear_probit(m~x+z, y~x+m)
#' print(est$estimates, digits=3)
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2023) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research, 34(1):67-84. Available at https://doi.org/10.1287/isre.2022.1113
linear_probit = function(form_linear, form_probit, data=NULL, par=NULL, method='BFGS', init=c('zero', 'unif', 'norm', 'default')[4], verbose=0){
    probit_linear(form_probit, form_linear, data=data, par=par, method=method, init=init, verbose=verbose)
}
