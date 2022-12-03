#' Recursive two-stage models to address endogeneity
#' @description  This package supports various recursive two-stage models to address the endogeneity issue. The details of the implemented models are discussed in Peng (2022). In a recursive two-stage model, the dependent variable of the first stage is also the endogenous variable of interest in the second stage. The endogeneity is captured by the correlation in the error terms of the two stages. \cr\cr
#' Recursive two-stage models can be used to address the endogeneity of treatment variables in observational study and the endogeneity of mediators in experiments. \cr\cr
#' The first-stage supports linear model, probit model, and Poisson lognormal model. The second-stage supports linear and probit models. These models can be used to address the endogeneity of continuous, binary, and count variables. When the endogenous variable is binary, it can be unobserved or partially unobserved, but the identification can be weak. \cr\cr
#' @section Functions:
#' bilinear: recursive bivariate linear model \cr \cr
#' biprobit: recursive bivariate probit model \cr \cr
#' biprobit_latent: recursive bivariate probit model with latent first stage \cr \cr
#' biprobit_partial: recursive bivariate probit model with partially observed first stage \cr \cr
#' linear-probit: recursive linear-probit model \cr \cr
#' probit_linear: recursive probit-linear model \cr \cr
#' probit_linear_latent: recursive probit-linear model with latent first stage \cr \cr
#' probit_linear_partial: recursive probit-linear model with partially observed first stage \cr \cr
#' probit_linearRE: recursive probit-linearRE model in which the second stage is a panel linear model with random effects \cr \cr
#' pln: Poisson lognormal (PLN) model \cr \cr
#' pln_linear: recursive PLN-linear model \cr \cr
#' pln_probit: recursive PLN-probit model \cr \cr
#' @docType package
#' @name endogeneity
#' @importFrom statmod gauss.quad
#' @importFrom stats binomial rnorm dnorm pnorm qnorm dpois dlogis plogis model.frame model.matrix model.response optim pchisq poisson runif lm glm coef sigma logLik
#' @importFrom maxLik maxLik numericGradient numericHessian
#' @importFrom pbivnorm pbivnorm
#' @importFrom MASS mvrnorm
#' @rawNamespace import(data.table)
#' @importFrom Rcpp sourceCpp
#' @useDynLib endogeneity
#' @references Peng, Jing. (2022) Identification of Causal Mechanisms from Randomized Experiments: A Framework for Endogenous Mediation Analysis. Information Systems Research (Forthcoming), Available at https://doi.org/10.1287/isre.2022.1113
NULL
#> NULL