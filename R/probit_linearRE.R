# library(endogeneity)
# library(maxLik)
# library(statmod)
# source('D:/Dropbox/Code/endogeneity/R/utilities.R')
#
#
# matVecProd <- function(m, v) {
#     .Call('_endogeneity_matVecProd', PACKAGE = 'endogeneity', m, v)
# }
#
# groupProd <- function(v, group) {
#     .Call('_endogeneity_groupProd', PACKAGE = 'endogeneity', v, group)
# }
#
# groupSum <- function(v, group) {
#     .Call('_endogeneity_groupSum', PACKAGE = 'endogeneity', v, group)
# }
#
# groupSumMat <- function(m, group) {
#     .Call('_endogeneity_groupSumMat', PACKAGE = 'endogeneity', m, group)
# }


LL_probit_linearRE = function(par,y,d,x,w,group,H=20,verbose=1){
    if(length(par) != ncol(x)+ncol(w)+3) stop("Number of parameters incorrect")
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    rule = gauss.quad(H, "hermite")

    Li = rep(0, length(d))
    for(k in 1:H){
        v = sqrt(2)*rule$nodes[k]
        omega = rule$weights[k] / sqrt(pi)

        s = (2*d-1)*(wa+rho*v)/sqrt(1-rho^2)
        u = xb + lambda*v

        log_phi = dnorm(y, u, sigma, log=T)
        phi_prod = exp(groupSum(log_phi, group))

        Li = Li + omega * pnorm(s) * phi_prod
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    LL = sum(log(Li))
    if(verbose>=1){
        writeLines(paste("==== Iteration ", endogeneity.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
        print(round(par,digits=3))
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)){
        if(verbose>=2) writeLines("NA or infinite likelihood, will try others")
        LL = -1e300
    }
    return (LL)
}


integrate_LL_probit_linearRE = function(y, d, xb, wa, group, lambda, sigma, rho, mu, scale, H){
    rule = gauss.quad(H, "hermite")
    Li = rep(0, length(d))
    for(k in 1:H){
        v = mu + sqrt(2)*rule$nodes[k]*scale
        omega = sqrt(2)*rule$weights[k]*exp(rule$nodes[k]^2) * dnorm(v) * scale

        s = (2*d-1)*(wa+rho*v)/sqrt(1-rho^2)
        u = xb + lambda*rep(v, times=diff(group))

        log_phi = dnorm(y, u, sigma, log=T)
        phi_prod = exp(groupSum(log_phi, group))

        Li = Li + omega * pnorm(s) * phi_prod
    }
    Li = pmax(Li, 1e-100) # in case that some Li=0
    sum(log(Li))
}


LL_probit_linearRE_AGQ = function(par,y,d,x,w,group,H=20,verbose=1){
    if(length(par) != ncol(x)+ncol(w)+3){
        print(names(par))
        print(colnames(x))
        print(colnames(w))
        stop(sprintf("Number of parameters incorrect, length(par)=%d, ncol(x)=%d, ncol(w)=%d", length(par), ncol(x), ncol(w)))
    }
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    rule = gauss.quad(H, "hermite")

    mu = endogeneity.env$mu
    scale = endogeneity.env$scale

    n_count = 1
    # Usually converge in a few iterations
    while(endogeneity.env$stopUpdate==F){
        mu_new = rep(0, length(d))
        scale_new = rep(0, length(d))
        Li = rep(0, length(d))
        for(k in 1:H){
            v = mu + sqrt(2)*rule$nodes[k]*scale
            omega = sqrt(2)*rule$weights[k]*exp(rule$nodes[k]^2) * dnorm(v) * scale

            s = (2*d-1)*(wa+rho*v)/sqrt(1-rho^2)
            u = xb + lambda*rep(v, times=diff(group))

            log_phi = dnorm(y, u, sigma, log=T)
            phi_prod = exp(groupSum(log_phi, group))

            Lik = omega * pnorm(s) * phi_prod
            Li = Li + Lik

            mu_new = mu_new + v*Lik
            scale_new = scale_new + v^2*Lik
        }
        Li = pmax(Li, 1e-100) # in case that some Li=0
        LL = sum(log(Li))
        mu_new = mu_new / Li
        scale_new = sqrt(pmax(scale_new / Li - mu_new^2, 1e-6))

        # Update only when LL improves
        if(!is.na(LL) && LL < endogeneity.env$LL) break

        # Stop on error
        if(any(is.na(scale_new)) || any(is.na(mu_new)) || any(scale_new<0) || n_count >= 100){
            if(verbose>=0){
                print(sprintf('--- Adaptive parameters not converging after %d inner and %d outer iterations ----', n_count, endogeneity.env$iter))
                print(par)
                print(mu_new)
                diff_mu = abs(mu_new - mu)/abs(mu)
                diff_scale = abs(scale_new - scale)/abs(scale)
                ix = which.max(diff_mu)
                print(c(mu[ix], mu_new[ix], scale[ix], scale_new[ix]))
                ix = which.max(diff_scale)
                print(c(mu[ix], mu_new[ix], scale[ix], scale_new[ix]))
                print(summary(diff_mu))
                print(summary(diff_scale))
            }
            break
        }

        # Check if mu and scale 99% converged
        if(mean(abs(mu_new - mu) < pmax(1e-4*abs(mu), 1e-6))>0.99
           && mean(abs(scale_new - scale) < pmax(1e-4*abs(scale), 1e-6))>0.99){
        # if(all(abs(mu_new - mu) < pmax(1e-4*abs(mu), 1e-6))
        #    && all(abs(scale_new - scale) < pmax(1e-4*abs(scale), 1e-6))){
            endogeneity.env$mu = mu_new
            endogeneity.env$scale = scale_new
            if(verbose>=2) print(sprintf('Adaptive parameters converged after %d iterations', n_count))
            break
        }
        n_count = n_count + 1
        mu = mu_new
        scale = scale_new
    }
    LL = integrate_LL_probit_linearRE(y, d, xb, wa, group, lambda, sigma, rho, endogeneity.env$mu, endogeneity.env$scale, H)

    if(verbose>=1){
        writeLines(paste("==== Iteration ", endogeneity.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
        print(round(par,digits=3))
    }
    addIter()
    if(is.na(LL) || !is.finite(LL)) LL = NA
    if(!is.na(LL) && LL > endogeneity.env$LL) {
        # stop updating when likelihood almost converges
        if(endogeneity.env$stopUpdate==F && (LL - endogeneity.env$LL < 1e-6 * abs(endogeneity.env$LL)) ) {
            endogeneity.env$stopUpdate=T
            if(verbose>=2) print('~~ Stop updating mu and scale now')
        }
        endogeneity.env$LL = LL
    }
    return (LL)
}

Gradient_probit_linearRE = function(par,y,d,x,w,group,H=20,verbose=1,variance=FALSE){
    if(length(par) != ncol(x)+ncol(w)+3){
        print(names(par))
        print(colnames(x))
        print(colnames(w))
        stop(sprintf("Number of parameters incorrect, length(par)=%d, ncol(x)=%d, ncol(w)=%d", length(par), ncol(x), ncol(w)))
    }
    alpha = par[1:ncol(w)]
    beta = par[ncol(w)+1:ncol(x)]
    lambda = exp(par[length(par)-2])
    sigma = exp(par[length(par)-1])
    tau = par[length(par)]
    rho = 1 - 2/(exp(tau)+1)

    wa = as.vector(w %*% alpha)
    xb = as.vector(x %*% beta)
    rule = gauss.quad(H, "hermite")

    mu = endogeneity.env$mu
    scale = endogeneity.env$scale

    Li = rep(0, length(d))
    dL_alpha = matrix(0, length(d), length(alpha) + 1) # also includes rho
    dL_beta = matrix(0, length(d), length(beta) + 2) # also includes lambda and sigma
    dw = cbind((2*d-1) / sqrt(1-rho^2) * w, rho=0)
    dx = cbind(x, lambda=0, sigma=0)
    for(k in 1:H){
        # GQ
        # v = sqrt(2)*rule$nodes[k]
        # omega = rule$weights[k] / sqrt(pi)

        # AGQ (same as GQ when mu=0 and scale=1)
        v = mu + sqrt(2)*rule$nodes[k]*scale
        omega = sqrt(2)*rule$weights[k]*exp(rule$nodes[k]^2) * dnorm(v) * scale

        s = (2*d-1)*(wa+rho*v)/sqrt(1-rho^2)
        u = xb + lambda*rep(v, times=diff(group))

        log_phi = dnorm(y, u, sigma, log=T)
        phi_prod = exp(groupSum(log_phi, group))

        Lik = omega * pnorm(s) * phi_prod
        Li = Li + Lik

        dx[, 1:length(beta)] = 1/sigma^2 * (y-u) * x
        dx[, 'lambda'] = 1/sigma^2 * (y-u) * rep(v, times=diff(group))
        dx[, 'sigma'] = ((y-u)^2 - sigma^2) / sigma^3
        dw[, 'rho'] = (2*d-1) * (v + wa * rho) / (1-rho^2)^1.5

        dL_alpha = dL_alpha + matVecProd(dw, omega * dnorm(s) * phi_prod)
        dL_beta = dL_beta + matVecProd(groupSumMat(dx, group), Lik)
    }
    dL = cbind(dL_alpha[, 1:length(alpha)], dL_beta, dL_alpha[, length(alpha)+1])
    Li = pmax(Li, 1e-100) # in case that some Li=0
    dL = matVecProd(dL, 1/Li)
    colnames(dL) = names(par)

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


#' Recursive Probit-LinearRE Model
#' @description A panel extension of the probit_linear model. The first stage is a probit model at the individual level. The second stage is a panel linear model at the individual-time level with individual-level random effects. The random effect is correlated with the error term in the first stage.\cr\cr
#' First stage (Probit):
#' \deqn{m_i=1(\boldsymbol{\alpha}'\mathbf{w_i}+u_i>0)}{m_i = 1(\alpha' * w_i + u_i > 0)}
#' Second stage (Panel linear model with individual-level random effects):
#' \deqn{y_{it} = \boldsymbol{\beta}'\mathbf{x_{it}} + {\gamma}m_i + \lambda v_i +\sigma \epsilon_{it}}{y_it = \beta' * x_it + \gamma * m_i + \lambda * v_i + \sigma * \epsilon_it}
#' Endogeneity structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}. \cr\cr
#' This model uses Adaptive Gaussian Quadrature to overcome numerical challenges with long panels. w and x can be the same set of variables. Identification can be weak if w are not good predictors of m. This model still works if the first-stage dependent variable is not a regressor in the second stage.
#' @param form_probit Formula for the probit model at the individual level
#' @param form_linear Formula for the linear model at the individual-time level
#' @param id group id, character if data  supplied or numerical vector if data not supplied
#' @param data Input data, must be a data.table object
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param method Optimization algorithm. Default is BFGS
#' @param stopUpdate Adaptive Gaussian Quadrature disabled if TRUE
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
#' * time: Time takes to estimate the model
#' * LR_stat: Likelihood ratio test statistic for \eqn{\rho=0}
#' * LR_p: p-value of likelihood ratio test
#' * iterations: number of iterations taken to converge
#' * message: Message regarding convergence status.
#'
#' Note that the list inherits all the components in the output of maxLik. See the documentation of maxLik for more details.
#' @md
#' @examples
#' library(MASS)
#' library(data.table)
#' N = 500
#' period = 5
#' obs = N*period
#' rho = -0.5
#' set.seed(100)
#'
#' e = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
#' e1 = e[,1]
#' e2 = e[,2]
#'
#' t = rep(1:period, N)
#' id = rep(1:N, each=period)
#' w = rnorm(N)
#' m = as.numeric(1+w+e1>0)
#' m_long = rep(m, each=period)
#'
#' x = rnorm(obs)
#' y = 1 + x + m_long + rep(e2, each=period) + rnorm(obs)
#'
#' dt = data.table(y, x, id, t, m=rep(m, each=period), w=rep(w, each=period))
#'
#' est = probit_linearRE(m~w, y~x+m, 'id', dt)
#' print(est$estimates, digits=3)
#' @export
#' @family endogeneity
#' @references Chen, H., Peng, J., Li, H., & Shankar, R. (2022). Impact of Refund Policy on Sales of Paid Information Services: The Moderating Role of Product Characteristics. Available at SSRN: https://ssrn.com/abstract=4114972.
probit_linearRE = function(form_probit, form_linear, id, data=NULL, par=NULL, method='BFGS', H=20, stopUpdate=F, init=c('zero', 'unif', 'norm', 'default')[4], verbose=0){
    # 1.1 Sort data by id
    if(is.null(data) && is.numeric(id)){
        ord = order(id)
        id = id[ord]
        group = c(0,cumsum(table(as.integer(factor(id)))))
        ind_data = NULL
    } else if(!is.null(data) && is.character(id)){
        id_var = id
        id = data[, id]
        ord = order(id)
        id = id[ord]
        group = c(0,cumsum(table(as.integer(factor(id)))))
        data = data[ord, ]
        # retain first row of each id
        ind_data = unique(data, by=id_var)
    } else {
        stop('data and id type mismatch. id should be a character if data  supplied, or a numerical vector if data not supplied')
    }


    if(!is.null(data)) data = data[ord,]

    # 1.1 parse linear formula
    mf = model.frame(form_linear, data=data, na.action=NULL, drop.unused.levels=TRUE)
    y = model.response(mf, "numeric")
    x = model.matrix(attr(mf, "terms"), data=mf)

    # 1.2 parse probit formula
    mf2 = model.frame(form_probit, data=ind_data, na.action=NULL, drop.unused.levels=TRUE)
    d = model.response(mf2, "numeric")
    w = model.matrix(attr(mf2, "terms"), data=mf2)
    # 1.3 Initialize parameters
    est_linear = lm(form_linear, data=data)
    par_linear = coef(summary(est_linear))[,1]
    est_probit = glm(form_probit, data=ind_data, family=binomial(link="probit"))
    par_probit = coef(summary(est_probit))[,1]
    names(par_linear) = paste0('linear.', names(par_linear))
    names(par_probit) = paste0('probit.', names(par_probit))
    par_linear[is.na(par_linear)] = 0
    par_probit[is.na(par_probit)] = 0
    par = c(par_probit, par_linear, log_lambda=0, log_sigma=0, tau=0)
    if(init=='unif') par = par - par + runif(length(par))
    if(init=='norm') par = par - par + rnorm(length(par))
    if(init=='zero') par = par - par
    # print(par)

    # 2. Estimation
    endogeneity.env$LL = -.Machine$double.xmax
    endogeneity.env$mu = rep(0, length(group)-1)
    endogeneity.env$scale = rep(1, length(group)-1)
    endogeneity.env$stopUpdate = stopUpdate
    endogeneity.env$iter = 1
    begin = Sys.time()

    # # use optim (hessian is forced to be symmetric)
    # res = optim(par=par, fn=LL_probit_linearRE_AGQ, gr=Gradient_probit_linearRE, method="BFGS", control=list(fnscale=-1), y=y, d=d, x=x, w=w, group=group, H=H, verbose=verbose, hessian = TRUE)

    # use maxLik (identical estimate with optim, but more reliable SE)
    res = maxLik(LL_probit_linearRE_AGQ, grad=Gradient_probit_linearRE, start=par, y=y, d=d, x=x, w=w, group=group, H=H, method=method, verbose=verbose)
    res$par = res$estimate
    res$n_obs = length(y)

    # 3. Compile results
    # res = getVarSE(res, verbose=verbose)
    gvar = Gradient_probit_linearRE(res$par,y,d,x,w,group,H,verbose=verbose-1,variance=TRUE)
    res = getVarSE(res, gvar=gvar, verbose=verbose)

    # res$num_g = numericGradient(LL_probit_linearRE,res$par,y=y, d=d, x=x, w=w, group=group, H=H)
    # cat('-------Gradient difference------\n')
    # print(res$num_g - gvar$g)

    res$estimates = transCompile(res, trans_vars=c(lambda='log_lambda', sigma='log_sigma', rho='tau'), trans_types=c('exp', 'exp', 'correlation'))
    res$LR_stat = 2 * ( res$LL - logLik(est_linear) - logLik(est_probit) )
    res$LR_p = 1 - pchisq(res$LR_stat, 1)
    res$ord = ord
    res$iter = endogeneity.env$iter
    res$mu = endogeneity.env$mu
    res$scale = endogeneity.env$scale

    if(verbose>=0){
        cat(sprintf('==== Converged after %d iterations, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
        # LR test still not correct as the linear model does not have RE
        # cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
        print(res$time <- Sys.time() - begin)
    }
    return (res)
}


# library(MASS)
# library(data.table)
# N = 1000
# period = 5
# obs = N*period
# rho = -0.5
# set.seed(100)
#
# e = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
# e1 = e[,1]
# e2 = e[,2]
#
# t = rep(1:period, N)
# id = rep(1:N, each=period)
# w = rnorm(N)
# d = as.numeric(1+w+e1>0)
# d_long = rep(d, each=period)
#
# x = rnorm(obs)
# y = 1 + x + d_long + rep(e2, each=period) + rnorm(obs)
#
# est = probit_linearRE(d~w, y~x+d_long, id, verbose=0, stopUpdate=F)
# # est2 = probit_linearRE(d~w, y~x+d_long, id, verbose=0, stopUpdate=T)
# print(est$estimates, digits=3)
# # est2$estimates
#
# # alternative way to call the function
# dt = data.table(y,x,id,t,d=rep(d, each=period),w=rep(w, each=period))
# est3 = probit_linearRE(d~w, y~x+d, 'id', dt, verbose=0, stopUpdate=F)
# est3$estimates

