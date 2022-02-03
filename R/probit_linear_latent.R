LL_probit_linear_latent = function(par,y,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    gamma = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    rho = 1 - 2/( exp(par[length(par)])+1 )

    w1 = as.vector(x1 %*% alpha)
    w2 = as.vector(x2 %*% beta)
    z1 = (y-w2-gamma)/sigma
    z0 = (y-w2)/sigma

    # errors can be amplified due to multiplication?
    Li = dnorm(y-w2-gamma, sd=sigma) * pnorm( (w1+z1*rho)/sqrt(1-rho^2) ) + dnorm(y-w2, sd=sigma) * pnorm( -(w1+z0*rho)/sqrt(1-rho^2) )

    LL = sum(log(Li))

    if(verbose>=1){
        cat(sprintf('==== Iteration %d: LL=%.5f ====\n', endogeneity.env$iter, LL))
        print(par,digits=3)
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

Q_probit_linear_latent = function(par,P1,y,x1,x2,verbose=1){
    if(length(par) != ncol(x1)+ncol(x2)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(x1)]
    beta = par[ncol(x1)+1:ncol(x2)]
    gamma = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    rho = 1 - 2/( exp(par[length(par)])+1 )

    w1 = as.vector(x1 %*% alpha)
    w2 = as.vector(x2 %*% beta)
    z1 = (y-w2-gamma)/sigma
    z0 = (y-w2)/sigma

    Li = P1 * ( dnorm(y-w2-gamma, sd=sigma, log=TRUE) + pnorm( (w1+z1*rho)/sqrt(1-rho^2), log.p=TRUE) ) + (1-P1) * ( dnorm(y-w2, sd=sigma, log=TRUE) + pnorm( -(w1+z0*rho)/sqrt(1-rho^2), log.p=TRUE) )

    LL = sum(Li)

    if(is.na(LL) || !is.finite(LL)) LL = NA
    return (LL)
}

EM_probit_linear = function(par,y,x1,x2,maxIter=200,tol=1e-8,tol_LL=1e-8,verbose=1,accu=1e4){
    par_last = par
    LL_last = -Inf
    for(i in 1:maxIter){
        if(verbose>=1) {
            cat('-------- EM Iteration',i,', LL=', LL_last,'-----\n')
            print(par)
        }
        # 1. E step
        alpha = par[1:ncol(x1)]
        beta = par[ncol(x1)+1:ncol(x2)]
        gamma = exp(par[length(par)-2])
        sigma = exp(par[length(par)-1])
        rho = 1 - 2/( exp(par[length(par)])+1 )

        w1 = as.vector(x1 %*% alpha)
        w2 = as.vector(x2 %*% beta)
        z1 = (y-w2-gamma)/sigma
        z0 = (y-w2)/sigma

        L1 = dnorm(y-w2-gamma, sd=sigma) * pnorm( (w1+z1*rho)/sqrt(1-rho^2) )
        L0 =  dnorm(y-w2, sd=sigma) * pnorm( -(w1+z0*rho)/sqrt(1-rho^2) )
        P1 = L1 / (L1+L0)


        # 2. M step
        est = optim(par=par, fn=Q_probit_linear_latent, gr=NULL, method="BFGS", control=list(factr=accu,fnscale=-1), P1=P1, y=y, x1=x1, x2=x2, verbose=0)
        par = est$par
        LL = LL_probit_linear_latent(par,y,x1,x2,verbose=0)

        # unlike EM_biprobit, the LL is highly sensitive to parameters and often increase during EM iterations
        if(all(abs(par_last-par)<abs(par_last)*tol) || abs(LL-LL_last)<abs(LL_last)*tol_LL) break
        par_last = par
        LL_last = LL
    }
    # Standard errors estimated based on the original likelihood function
    est$hessian = numericHessian(LL_probit_linear_latent, t0=par, y=y, x1=x1, x2=x2, verbose=0)
    # est$hessian2 = numericHessian(Q_probit_linear_latent, t0=par, P1=P1, y=y, x1=x1, x2=x2, verbose=0)
    est$gradient = numericGradient(LL_probit_linear_latent, t0=par, y=y, x1=x1, x2=x2, verbose=0)
    # est$g2 = numericGradient(Q_probit_linear_latent, t0=par, P1=P1, y=y, x1=x1, x2=x2, verbose=0)
    rownames(est$hessian) = colnames(est$hessian) = names(est$gradient) = names(par)
    est$maximum = LL
    est$iter = i
    est
}

#' Recursive Probit-Linear Model with Latent First Stage
#' @description The first stage is a probit model with unobserved dependent variable, the second stage is a linear model that includes the first-stage dependent variable as a regressor.
#' @param form_probit Formula for the first-stage probit model, in which the dependent variable is latent
#' @param form_linear Formula for the second stage linear model. The latent dependent variable of the first stage is automatically added as a regressor in this model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
#' @param EM Whether to maximize likelihood use the Expectation-Maximization algorithm. EM is slower but more robust
#' @param method Optimization algorithm. Default is BFGS
#' @param accu 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim
#' @param maxIter max iterations for EM algorithm
#' @param tol tolerance for convergence of EM algorithm
#' @param tol_LL tolerance for convergence of likelihood
#' @param verbose Level of output during estimation. Lowest is 0.
#' @return A list containing the results of the estimated model
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
#' y1 = as.numeric(1 + x + z + e1 > 0)
#' y2 = 1 + x + z + y1 + e2
#' est = probit_linear(y1~x+z, y2~x+z+y1)
#' est$estimates
#'
#' est_latent = probit_linear_latent(~x+z, y2~x+z)
#' est_latent$estimates
#' }
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at SSRN: https://ssrn.com/abstract=3494856
probit_linear_latent = function(form_probit, form_linear, data=NULL, EM=TRUE, par=NULL, method='BFGS', verbose=0, accu=1e4, maxIter=500, tol=1e-6, tol_LL=1e-8){
    # 1.1 parse ~x1
    mf1 = model.frame(form_probit, data=data, na.action=NULL, drop.unused.levels=TRUE)
    x1 = model.matrix(attr(mf1, "terms"), data=mf1)
    # 1.2 parse y~x2
    mf2 = model.frame(form_linear, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf2, "numeric")
    x2 = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    par = rep(0, ncol(x1) + ncol(x2) + 3)
    names(par) = c(paste0('probit.', colnames(x1)), paste0('linear.', colnames(x2)), 'log_gamma', 'log_sigma', 'tau')
    # print(par)
    # 2. Estimation
    resetIter()
    begin = Sys.time()

    if(EM==TRUE){
        res = EM_probit_linear(par,y,x1,x2,maxIter=maxIter,tol=tol,tol_LL=tol_LL,verbose=verbose,accu=accu)
    } else {
        # res = optim(par=par, fn=LL_probit_linear_latent, gr=NULL, method=method, control=list(factr=accu,fnscale=-1), y=y, x1=x1, x2=x2, verbose=verbose, hessian = TRUE)
        # res$g = numericGradient(LL_probit_linear_latent, t0=res$par, y=y, x1=x1, x2=x2, verbose=0)

        # use maxLik (identical estimate with optim, but more reliable SE)
        res = maxLik(LL_probit_linear_latent, start=par, method=method, y=y, x1=x1, x2=x2, verbose=verbose)
        res$par = res$estimate
    }

    # 3. Compile results
    res$n_obs = length(y)
    res = getVarSE(res, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(gamma='log_gamma', sigma='log_sigma', rho='tau'), trans_types=c('exp', 'exp', 'correlation'))
    res$iter = endogeneity.env$iter

    if(verbose>0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}

