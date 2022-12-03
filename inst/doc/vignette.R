## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Table, echo=FALSE--------------------------------------------------------
tb = data.frame(
    Model = c('biprobit', 'biprobit_latent', 'biprobit_partial', 'probit_linear', 'probit_linear_latent', 'probit_linear_partial', 'probit_linearRE', 'pln_linear', 'pln_probit'),
    `First Stage` = c(rep('probit', 7), rep('pln', 2)),
    `Second Stage` = c(rep('probit',3), rep('linear', 3), 'linearRE', 'linear', 'probit'),
    `Endogenous Variable` = c(rep(c('binary', 'binary (unobserved)', 'binary (partially observed)'), 2), 'binary', rep('count', 2)),
    `Outcome Variable` = c(rep('binary',3), rep('continuous', 5), 'binary')
)
knitr::kable(tb, col.names = gsub("[.]", " ", names(tb)), caption='**Table 1. Recursive Two-Stage Models Supported by the Endogeneity Package**')

## -----------------------------------------------------------------------------
library(endogeneity)
example(probit_linear, prompt.prefix=NULL)

