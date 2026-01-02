#' MLE of Contaminated Negative Binomial Regression Model
#' Extended to Include Covariate Estimation of all Distribution Parameters.
#'
#'
#' MLE of the contaminated negative binomial regression model via an Expectation-Maximization algorithm.
#'
#' @param formula An object of class 'formula': a symbolic description of the model to be fitted.
#' @param alpha.formula A formula for the dispersion parameter (sigma). Defaults to ~1.
#' @param delta.formula A formula for the proportion of extreme values parameter (delta). Defaults to ~1.
#' @param eta.formula A formula for the inflation parameter (eta). Defaults to ~1.
#' @param data A mandatory data frame containing the variables in the model.
#' @param start Vector of initial values. If NULL, then values produced by glm.nb are used as initial values
#' @param method Optimization method to be used. Default is "Nelder-Mead". Other options include "CG" and "BFGS".
#' @param reltol Relative convergence tolerance in each optimization process. Defaults to 1e-15.
#' @param maxit The number of inner iterations in each optimization process. Defaults to 10000.
#' @param em.tol  The EM convergence tolerance. Defaults to 1e-10.
#' @param em.maxit  The number of EM iterations. Defaults to 1000.
#'
#' @details The \code{ml.cnb} function fits the contaminated negative binomial regression model see (https://github.com/arnootto/cNB/tree/master).
#'
#' @return A list of elements:
#'    \item{results}{A data frame with parameter estimates and standard errors. }
#'    \item{alpha}{Maximum likelihood estimate of alpha.}
#'    \item{delta}{Maximum likelihood estimate of delta.}
#'    \item{beta}{Maximum likelihood estimates of the regression coefficients.}
#'    \item{mu}{Maximum likelihood estimates of the mean parameter.}
#'    \item{X}{The design matrix for mu.}
#'    \item{U}{The design matrix for alpha.}
#'    \item{V}{The design matrix for delta.}
#'    \item{Z}{The design matrix for eta.}
#'    \item{y}{The response variable.}
#'    \item{loglike}{The log-likelihood value at convergence.}
#'    \item{AIC}{Akaike Information Criterion (AIC) for the fitted model.}
#'    \item{BIC}{Bayesian Information Criterion (BIC) for the fitted model.}
#'    \item{HQIC}{Hannan-Quinn Information Criterion (HQIC) for the fitted model.}
#'
#' @import MASS
#' @import cNB
#'
#' @examples
#' library(AER)
#'
#' # Load NMES1988 dataset
#' data("NMES1988")
#' head(NMES1988)
#'
#' # Construct design matrix
#' x <- model.matrix(~ visits + nvisits + ovisits + novisits + emergency + hospital +
#'                     health + chronic + adl + region + afam + gender + married +
#'                     school + income + employed + insurance + medicaid,
#'                   data = NMES1988)
#' x <- as.data.frame(x)
#'
#' options(contrasts = c("contr.treatment", "contr.poly"))
#'
#'
#'
#'
#' # Regression on mean, dispersion, delta, and eta
#' cnb_ex_trips_mu_alpha_delta_eta <- ml.ex.cnb(formula = visits ~ hospital + insuranceyes +
#'                                                           healthpoor + chronic + adllimited +
#'                                                           gendermale + school + afamyes,
#'                                             alpha.formula = ~ adllimited,
#'                                             delta.formula = ~ nvisits,
#'                                             eta.formula = ~ chronic,
#'                                             data = x,
#'                                             method = "BFGS")
#'
#' # Parameter estimates
#' round(cnb_ex_trips_mu_alpha_delta_eta$results, 4)
#'
#' # Information criteria and log-likelihood
#' round(cnb_ex_trips_mu_alpha_delta_eta$AIC, 4)
#' round(cnb_ex_trips_mu_alpha_delta_eta$BIC, 4)
#' round(cnb_ex_trips_mu_alpha_delta_eta$HQIC, 4)
#' round(cnb_ex_trips_mu_alpha_delta_eta$loglike, 4)
#'
#' @export

ml.ex.cnb <- function(formula, alpha.formula=~1, delta.formula=~1, eta.formula=~1, data, start = NULL, method = "Nelder-Mead", reltol=1e-15, maxit=1000, hessian=T) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  betaX <- model.matrix(formula, data = data)  # design matrix for mean
  alphaU <- model.matrix(alpha.formula, data = data)  # design matrix for sigma
  deltaV <- model.matrix(delta.formula, data = data)  # design matrix for delta
  etaZ <- model.matrix(eta.formula, data = data)
  nb2.reg.ml <- function(b.hat, X, U, V, Z, y) {
    beta.hat <- b.hat[1:ncol(X)] #coefficients for the mean (mu)
    theta.hat <- b.hat[(ncol(X)+1):(ncol(X)+ncol(U))] #coefficients for alpha
    gamma.hat <-b.hat[(ncol(X)+ncol(U)+1):(ncol(X)+ncol(U)+ncol(V))]
    lambda.hat <- b.hat[(ncol(X)+ncol(U)+ncol(V)+1):(ncol(X)+ncol(U)+ncol(V)+ncol(Z))]

    xb.hat <- X %*% beta.hat  # mean regression
    mu.hat <- exp(xb.hat)

    ua.hat <- U %*% theta.hat  # alpha regression
    alpha.hat <- exp(ua.hat)

    vg.hat <- V %*% gamma.hat  # delta regression
    delta.hat <- exp(vg.hat) / (1 + exp(vg.hat))

    zl.hat <- Z %*% lambda.hat  # eta regression
    eta.hat <- exp(zl.hat) 

    ll <- sum(dcnbinom	(x = y, mu = mu.hat, alpha = alpha.hat, delta = delta.hat, eta=eta.hat, log=T))
    return(ll)
  }

  if (is.null(start)) {#initial param
    nb <- MASS::glm.nb(formula = formula, data = na.omit(data), control = glm.control(maxit = 2500))
    nb$loglike
    beta.hat <- nb$coefficients
    if (anyNA(beta.hat)) {
      warning("NA detected in glm.nb coefficients. Replacing with small fallback values.")
      beta.hat[is.na(beta.hat)] <- runif(sum(is.na(beta.hat)))
    }
    theta.hat <- rep(log(0.1), ncol(alphaU))
    gamma.hat = rep(log(0.0055/(1-0.0055)),ncol(deltaV))
    lambda.hat = rep(log(0.006), ncol(etaZ))


    start <- c(beta.hat, theta.hat, gamma.hat,lambda.hat)

  }

  print(start)
  summary(start)
  any(!is.finite(start))

  fit <- optim(par = start,
               fn = nb2.reg.ml,
               X = betaX,
               U = alphaU,
               V = deltaV,
               Z = etaZ,
               y = y,
               method = method,
               control = list(fnscale = -1, maxit = maxit, reltol = reltol),
               hessian = hessian
  )
  if (fit$convergence != 0) warning("Optimization may not have converged.")

  beta.hat <- fit$par[1:ncol(betaX)] #coefficients for the mean (mu)
  theta.hat <- fit$par[(ncol(betaX)+1):(ncol(betaX)+ncol(alphaU))] #coefficients for alpha
  gamma.hat <-fit$par[(ncol(betaX)+ncol(alphaU)+1):(ncol(betaX)+ncol(alphaU)+ncol(deltaV))]
  lambda.hat <- fit$par[(ncol(betaX)+ncol(alphaU)+ncol(deltaV)+1):(ncol(betaX)+ncol(alphaU)+ncol(deltaV)+ncol(etaZ))]

  xb.hat <- betaX %*% beta.hat  #' mean regression
  mu.hat <- exp(xb.hat)

  ua.hat <- alphaU %*% theta.hat  #' alpha regression
  alpha.hat <- exp(ua.hat)

  vg.hat <- deltaV %*% gamma.hat  #' delta regression
  delta.hat <- exp(vg.hat) / (1 + exp(vg.hat))

  zl.hat <- etaZ %*% lambda.hat  #' eta regression
  eta.hat <- exp(zl.hat) 


  lc=sum(dcnbinom(y, mu=mu.hat, alpha = alpha.hat, delta=delta.hat, eta = eta.hat, log = T))


  if (hessian==T)
  {cov.mat <- tryCatch(solve(-fit$hessian), error = function(e) MASS::ginv(-fit$hessian))
  std.errors <- sqrt(diag(cov.mat))
  tvalue <- fit$par/std.errors
  pval <- 2*pt(abs(tvalue),df=nrow(data)-length(fit$par),lower.tail = F)
  results <- data.frame(
    Estimate = c(fit$par),
    `Std.Error` = std.errors,
    `t.value` = tvalue,
    `P.t` = round(pval,5)
  )}
  if(hessian==F){
    results <- data.frame(
      Estimate = c(fit$par))
  }
  rownames(results) <- c(paste0("beta_", colnames(betaX)), paste0("theta_", colnames(alphaU)), paste0("gamma_", colnames(deltaV)), paste0("lambda_", colnames(etaZ)))


  AIC <- -2*lc+nrow(results)*2
  BIC <- -2*lc+nrow(results)*log(nrow(data))
  HQIC <- -2*lc +nrow(results)*2*log(log(nrow(data)))

  return(list(
    results = results,
    beta = beta.hat,
    theta = theta.hat,
    gamma = gamma.hat,
    lambda = lambda.hat,
    mu = mu.hat,
    alpha = alpha.hat,
    delta = delta.hat,
    eta = eta.hat,
    X = betaX,
    U = alphaU,
    V = deltaV,
    Z = etaZ,
    y = y,
    loglike = lc,
    AIC = AIC,
    BIC = BIC,
    HQIC = HQIC
  ))
}
