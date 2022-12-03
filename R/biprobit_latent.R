LL_biprobit_latent = function(par,y,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    gamma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    q = 2*y-1
    w1 = as.vector(x1 %*% alpha)
    w2 = as.vector(x2 %*% beta)

    Li = pbivnorm(w1,q*(w2+gamma),q*rho) + pbivnorm(-w1,q*w2,-q*rho)

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

Q_biprobit_latent = function(par,P1,y,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+2) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    gamma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    q = 2*y-1
    w1 = as.vector(x1 %*% alpha)
    w2 = as.vector(x2 %*% beta)

    L1 = pmax(1e-100, pbivnorm(w1,q*(w2+gamma),q*rho))
    L0 = pmax(1e-100, pbivnorm(-w1,q*w2,-q*rho))

    Li = P1 * log(L1) + (1-P1) * log(L0)
    LL = sum(Li)

    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

EM_biprobit = function(par,y,x1,x2,maxIter=200,tol=1e-6,tol_LL=1e-8,verbose=1){
    par_last = par
    LL_last = -Inf
    for(i in 1:maxIter){
        if(verbose>=1) {
            cat('-------- EM Iteration',i,', LL=', LL_last,'-----\n')
            print(par)
        }
        # 1. E step
        if(length(par) != ncol(x1)+ncol(x2)+2) stop("Number of parameters incorrect")
        alpha = par[1:ncol(x1)]
        beta = par[ncol(x1)+1:ncol(x2)]
        gamma = exp(par[length(par)-1])
        tau = par[length(par)]
        rho = 1 - 2/(exp(tau)+1)

        q = 2*y-1
        w1 = as.vector(x1 %*% alpha)
        w2 = as.vector(x2 %*% beta)

        L1 = pmax(1e-100, pbivnorm(w1,q*(w2+gamma),q*rho))
        L0 = pmax(1e-100, pbivnorm(-w1,q*w2,-q*rho))
        P1 = L1 / (L1+L0)


        # 2. M step
        est = optim(par=par, fn=Q_biprobit_latent, gr=NULL, method="BFGS", control=list(fnscale=-1), P1=P1, y=y, x1=x1, x2=x2, verbose=0)
        par = est$par
        LL = LL_biprobit_latent(par,y,x1,x2,verbose=0)

        if(all(abs(par_last-par)<abs(par_last)*tol) || abs(LL-LL_last)<abs(LL_last)*tol_LL) break
        par_last = par
        LL_last = LL
    }
    est$hessian = numericHessian(LL_biprobit_latent, t0=par, y=y, x1=x1, x2=x2, verbose=0)
    est$gradient = numericGradient(LL_biprobit_latent, t0=par, y=y, x1=x1, x2=x2, verbose=0)
    rownames(est$hessian) = colnames(est$hessian) = names(est$gradient) = names(par)
    est$maximum = LL
    est$iter = i
    est
}


#' Recursive Bivariate Probit Model with Latent First Stage
#' @description Estimate two probit models with bivariate normally distributed error terms, in which the dependent variable of the first stage model is unobserved.\cr\cr
#' First stage (Probit, \eqn{m_i^*} is unobserved):
#' \deqn{m_i^*=1(\boldsymbol{\alpha}'\mathbf{w_i}+u_i>0)}{m_i^* = 1(\alpha' * w_i + u_i > 0)}
#' Second stage (Probit):
#' \deqn{y_i = 1(\boldsymbol{\beta}'\mathbf{x_i} + {\gamma}m_i^* + \sigma v_i>0)}{y_i = 1(\beta' * x_i + \gamma * m_i^* + \sigma * v_i > 0)}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' w and x can be the same set of variables. The identification of this model is generally weak, especially if w are not good predictors of m. \eqn{\gamma} is assumed to be positive to ensure that the model estimates are unique.
#' @param form1 Formula for the first probit model, in which the dependent variable is unobserved. Use a formula like ~w to avoid specifying the dependent variable.
#' @param form2 Formula for the second probit model, the latent dependent variable of the first stage is automatically added as a regressor in this model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param EM Whether to maximize likelihood use the Expectation-Maximization (EM) algorithm, which is slower but more robust. Defaults to FLASE, but should change to TRUE is the model has convergence issues.
#' @param method Optimization algorithm. Default is BFGS
#' @param maxIter max iterations for EM algorithm
#' @param tol tolerance for convergence of EM algorithm
#' @param tol_LL tolerance for convergence of likelihood
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
#' * gradient: Gradient function at maximum
#' * hessian: Hessian matrix at maximum
#' * gtHg: \eqn{g'H^-1g}, where H^-1 is simply the covariance matrix. A value close to zero (e.g., <1e-3 or 1e-6) indicates good convergence.
#' * LL or maximum: Likelihood
#' * AIC: AIC
#' * BIC: BIC
#' * n_obs: Number of observations
#' * n_par: Number of parameters
#' * iterations: number of iterations taken to converge
#' * message: Message regarding convergence status.
#'
#' Note that the list inherits all the components in the output of maxLik. See the documentation of maxLik for more details.
#' @md
#' @examples
#' \donttest{
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
#'
#' est_latent = biprobit_latent(~x+z, y~x+z)
#' print(est_latent$estimates, digits=3)
#' }
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at https://doi.org/10.1287/isre.2022.1113
biprobit_latent = function(form1, form2, data=NULL, EM=FALSE, par=NULL, method='BFGS', verbose=0, maxIter=500, tol=1e-5, tol_LL=1e-6){
    # 1.1 parse ~x1
    mf1 = model.frame(form1, data=data, na.action=NULL, drop.unused.levels=TRUE)
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y~x2
    mf2 = model.frame(form2, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    par = rep(0, ncol(x1) + ncol(x2) + 2)
    names(par) = c(paste0('1', colnames(x1)), colnames(x2), 'log_gamma', 'tau')
    # 2. Estimation
    resetIter()
    begin = Sys.time()

    if(EM==TRUE){
        res = EM_biprobit(par,y,x1,x2,maxIter=maxIter,tol=tol,tol_LL=tol_LL,verbose=verbose)
    } else {
        # for biprobit_latent, maximize LL directly seems to work very well
        # The gtHg often converges to zero and the likelihood do not increase during EM
        # res = optim(par=par, fn=LL_biprobit_latent, gr=NULL, method=method, control=list(fnscale=-1), y=y, x1=x1, x2=x2, verbose=verbose, hessian = TRUE)
        # res$g = numericGradient(LL_biprobit_latent, t0=res$par, y=y, x1=x1, x2=x2, verbose=0)

        # use maxLik (identical estimate with optim, but more reliable SE)
        res = maxLik(LL_biprobit_latent, start=par, method=method, y=y, x1=x1, x2=x2, verbose=verbose)
        res$par = res$estimate
    }

    # 3. Compile results
    res$n_obs = length(y)
    res = getVarSE(res, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(gamma='log_gamma', rho='tau'), trans_types=c('exp', 'correlation'))
    res$iter = endogeneity.env$iter

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
