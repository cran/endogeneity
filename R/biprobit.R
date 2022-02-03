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
#' @description Estimate two probit models with bivariate normally distributed error terms. This command still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form1 Formula for the first probit model
#' @param form2 Formula for the second probit model
#' @param data Input data, a data frame
#' @param par Starting values for estimates
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
#' y2 = as.numeric(1 + x + z + y1 + e2 > 0)
#'
#' est = biprobit(y1~x+z, y2~x+z+y1)
#' est$estimates
#' @export
#' @family endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at SSRN: https://ssrn.com/abstract=3494856
biprobit = function(form1, form2, data=NULL, par=NULL, method='BFGS', verbose=0, accu=1e4){
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
    # res$optim = optim(par=par, fn=LL_biprobit, gr=Gradient_biprobit, method=method, control=list(factr=accu,fnscale=-1), y1=y1, y2=y2, x1=x1, x2=x2, verbose=verbose, hessian = TRUE)

    # 3. Compile results
    gvar = Gradient_biprobit(res$par,y1,y2,x1,x2,verbose=verbose-1,variance=TRUE)
    res = getVarSE(res, gvar=gvar, verbose=verbose)
    res$estimates = transCompile(res, trans_vars=c(rho='tau'), trans_types=c('correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(probit1) - logLik(probit2) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$iter = endogeneity.env$iter

    if(verbose>0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f **** \n', res$iterations, res$LL, res$gtHg))
        cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}
