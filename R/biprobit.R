LL_biprobit = function(par,y1,y2,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+1) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    q1 = 2*y1-1
    q2 = 2*y2-1
    w1 = as.vector(x1 %*% alpha) * q1
    w2 = as.vector(x2 %*% beta) * q2
    rho_star = q1*q2*rho

    Li = pbivnorm(w1,w2,rho_star)
    Li = pmax(Li, 0) # make sure Li is nonnegative
    LL = sum(log(Li))

    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}


Gradient_biprobit = function(par,y1,y2,x1,x2,verbose=1,variance=FALSE){
    if(length(par) != ncol(x1)+ncol(x2)+1) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    q1 = 2*y1-1
    q2 = 2*y2-1
    w1 = as.vector(x1 %*% alpha) * q1
    w2 = as.vector(x2 %*% beta) * q2
    rho_star = q1*q2*rho

    Phi2 = pbivnorm(w1,w2,rho_star)
    phi2 = dbinorm(w1,w2,rho_star)

    g1 = dnorm(w1) * pnorm((w2-rho_star*w1)/sqrt(1-rho_star^2))
    g2 = dnorm(w2) * pnorm((w1-rho_star*w2)/sqrt(1-rho_star^2))


    # rowwise multiplication on matrix x1
    dL1 = (q1*g1/Phi2) * x1
    dL2 = (q2*g2/Phi2) * x2
    dL_rho = q1*q2*phi2/Phi2
    dL_tau = dL_rho * (2 * exp(tau) / (exp(tau)+1)^2)

    dL = cbind(dL1, dL2, dL_tau)
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

#' Recusrive Bivariate Probit Model
#' @description Estimate two probit models with bivariate normally distributed error terms.\cr\cr
#' First stage (Probit):
#' \deqn{m_i=1(\boldsymbol{\alpha}'\mathbf{w_i}+u_i>0)}{m_i = 1(\alpha' * w_i + u_i > 0)}
#' Second stage (Probit):
#' \deqn{y_i = 1(\boldsymbol{\beta}'\mathbf{x_i} + {\gamma}m_i + \sigma v_i>0)}{y_i = 1(\beta' * x_i + \gamma * m_i + \sigma * v_i > 0)}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' w and x can be the same set of variables. Identification can be weak if w are not good predictors of m. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form1 Formula for the first probit model
#' @param form2 Formula for the second probit model
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
#' m = as.numeric(1 + x + z + e1 > 0)
#' y = as.numeric(1 + x + z + m + e2 > 0)
#'
#' est = biprobit(m~x+z, y~x+z+m)
#' print(est$estimates, digits=3)
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at https://doi.org/10.1287/isre.2022.1113
biprobit = function(form1, form2, data=NULL, par=NULL, method='BFGS', verbose=0){
    # 1.1 parse y1~x1
    mf1 = model.frame(form1, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y1 = model.response(mf1, "numeric")
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y2~x2
    mf2 = model.frame(form2, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y2 = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    probit1 = glm(form1, data=data, family=binomial(link="probit"))
    probit2 = glm(form2, data=data, family=binomial(link="probit"))
    par_probit1 = coef(summary(probit1))[,1]
    par_probit2 = coef(summary(probit2))[,1]
    names(par_probit1) = paste0('1', names(par_probit1))
    par = c(par_probit1, par_probit2, tau=0)
    # print(par)

    # 2. Estimation
    resetIter()
    begin = Sys.time()

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_biprobit, grad=Gradient_biprobit, start=par, method=method, y1=y1, y2=y2, x1=x1, x2=x2, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y1)

    # use optim
    # res$optim = optim(par=par, fn=LL_biprobit, gr=Gradient_biprobit, method=method, control=list(fnscale=-1), y1=y1, y2=y2, x1=x1, x2=x2, verbose=verbose, hessian = TRUE)

    # 3. Compile results
    gvar = Gradient_biprobit(res$par,y1,y2,x1,x2,verbose=verbose-1,variance=TRUE)
    res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(rho='tau'), trans_types=c('correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(probit1) - logLik(probit2) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$iter = endogeneity.env$iter

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
