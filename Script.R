#----TCC----

# Pacotes:
library(geoR)
library(geoRglm)
library(MASS)
library(sp)
library(lme4)
library(bbmle)
library(parallel)
library(doParallel)
library(xtable)

# Funções:

monta.sigma <- function(cov.pars, cov.model, nugget = 0, kappa, mat.dist){
    Sigma <- varcov.spatial(dists.lowertri = mat.dist,
                            cov.model = cov.model,kappa = kappa,
                            nugget = nugget, cov.pars = cov.pars)
    return(Sigma)
}

gauss.mult <- function(b, det.Sigma, inv.Sigma){
    n <- length(b)
    dens <- (-n/2) * log(2 * pi) - 0.5 * det.Sigma -
        0.5 * t(b) %*% inv.Sigma %*% b
    return(dens)
}

newton.raphson_old <- function(initial, escore, hessiano,
                               tol=0.0001, max.iter, n.dim, ...){
    solucao <- matrix(NA, max.iter, n.dim)
    solucao[1,] <- initial
    for(i in 2:max.iter){
        HSS <- hessiano(initial, ...)
        ESC <- t(escore(initial, ...))
        solucao[i,] <- initial - solve(HSS,ESC)
        initial <- solucao[i,]
        tolera <- abs(solucao[i,] - solucao[i-1,])
        if(all(tolera < tol) == TRUE)break
    }
    saida <- list()
    saida[1][[1]] <- HSS
    saida[2][[1]] <- initial
    saida <<- initial
    return(saida)
}

newton.raphson <- function(initial, escore, hessiano, tol=0.0001, max.iter, n.dim, ...) {
  solucao <- matrix(NA, max.iter, n.dim)
  solucao[1,] <- initial
  for (i in 2:max.iter) {
    HSS <- hessiano(initial, ...)
    ESC <- t(escore(initial, ...))
    solucao[i,] <- initial - solve(HSS, ESC)
    initial <- solucao[i,]
    tolera <- abs(solucao[i,] - solucao[i-1,])
    if (all(tolera < tol)) break
  }
  saida <- list(HSS = HSS, solution = initial)
  return(saida)
}

# Função para o integrando
Q.b <- function(b, Xbeta, Y, det.Sigma, inv.Sigma, family,
                prec = NULL, ntrial = NULL){
  eta <- Xbeta + b
  dens <- switch (
      family,
      "poisson" = sum(dpois(Y, lambda = exp(eta),log=TRUE)) +
          gauss.mult(b, det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma),
      "binomial" = sum(dbinom(Y, size = ntrial,
                              prob = (1/(1 + exp(- eta))), log = TRUE)) +
          gauss.mult(b, det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma),
      "negative.binomial" = sum(dnbinom(Y, size = prec, mu = exp(eta),
                                        log = TRUE)) +
          gauss.mult(b,det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma),
      "gamma" = sum(dgamma(Y, shape = prec, scale = exp(eta)/prec,
                           log=TRUE)) +
          gauss.mult(b,det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma),
      "beta" = sum(dbeta(Y, shape1 = (1/(1 + exp(- eta))) * prec,
                         shape2 = (1 - (1/(1 + exp(- eta)))) * prec,
                         log = TRUE)) +
          gauss.mult(b, det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma),
      "geometric" = sum(dnbinom(Y, size = 1, prob = (1/(1 + exp(eta))), log = TRUE)) +
          gauss.mult(b, det.Sigma = det.Sigma,
                     inv.Sigma = inv.Sigma)
  )
  return(as.numeric(dens))
}

# Função para o gradiente:
Q.b.grad <- function(b, Xbeta, Y, det.Sigma, inv.Sigma, family,
                     prec = NULL, ntrial = NULL){
    grad <- switch(
        family,
        "poisson" = {
            t((Y - exp(Xbeta + b))) - t(b) %*% inv.Sigma
        },
        "binomial" = {
            b1 <- ntrial * (1/(1 + exp(-(Xbeta + b))))
            t(Y - b1) - t(b) %*% inv.Sigma
        },
        "negative.binomial" = {
            et1 <- exp(Xbeta)
            et2 <- exp(b)
            et12 <- et1 * et2
            t((Y - (prec + Y) * (et12 * (et12 + prec)^-1))) -
                t(b) %*% inv.Sigma
        },
        "gamma" = {
            p1 <- -prec + prec * (exp(-Xbeta - b)) * Y
            t(p1) - t(b) %*% inv.Sigma
        },
        "beta" = {
            eta <- Xbeta + b
            mu <- 1/(1 + exp(- eta))
            D.mu.b <- exp(-eta) / ((1 + exp(-eta))^2)
            dg1 <- digamma((1 - mu) * prec)
            dg2 <- digamma(mu * prec)
            logY <- log(Y / (1 - Y))
            part1 <- D.mu.b * prec * (dg1 - dg2 + logY)
            t(part1) - t(b) %*% inv.Sigma
        },
        "geometric" = {
            eta <- exp(Xbeta + b)
            t((Y - (1 + Y) * (eta * (eta + 1)^-1))) -
                t(b) %*% inv.Sigma
        }
    )
}

# Função para o hessiano:
Q.b.hess <- function(b, Xbeta, Y, det.Sigma, inv.Sigma,
                     family, prec = NULL, ntrial = NULL) {
    Hess <- switch(
        family,
        "poisson" = {
            eta <- exp(Xbeta + b)
            diag(inv.Sigma) <- eta + diag(inv.Sigma)
            -inv.Sigma
        },
        "binomial" = {
            b1 <- 1/(1 + exp(-(Xbeta + b)))
            b2 <- exp(2 * Xbeta + 2 * b) / ((1 + exp(Xbeta + b))^2)
            D <- b1 - b2
            diag(inv.Sigma) <- ntrial * D + diag(inv.Sigma)
            -inv.Sigma
        },
        "negative.binomial" = {
            et1 <- exp(Xbeta)
            et2 <- exp(b)
            et12 <- et1 * et2
            p1 <- et12 * ((et12 + prec)^-1)
            p2 <- p1^2
            D <- p1 - p2
            diag(inv.Sigma) <- (prec + Y) * D + diag(inv.Sigma)
            -inv.Sigma
        },
        "gamma" = {
            p2 <- prec * Y * exp(-Xbeta - b)
            diag(inv.Sigma) <- p2 + diag(inv.Sigma)
            -inv.Sigma
        },
        "beta" = {
            eta <- Xbeta + b
            mu <- 1/(1 + exp(-eta))
            p1 <- mu * prec
            p2 <- (1 - mu) * prec
            med <- mu * (1 - mu)
            logY <- log(Y / (1 - Y))
            d2 <- (1 - mu)^2 - mu^2
            part1 <- -prec^2 * (trigamma(p1) + trigamma(p2)) * med
            part2 <- prec * (digamma(p2) - digamma(p1) + logY) * d2
            part3 <- (part1 + part2) * med
            diag(inv.Sigma) <- -part3 + diag(inv.Sigma)
            -inv.Sigma
        },
        "geometric" = {
            et12 <- exp(Xbeta + b)
            p1 <- et12 * ((et12 + 1)^-1)
            p2 <- p1^2
            D <- p1 - p2
            diag(inv.Sigma) <- (1 + Y) * D + diag(inv.Sigma)
            -inv.Sigma
        }
    )
    return(Hess)
}

laplace_old <- function(Q.b, gr, hess, method.integrate, otimizador, n.dim, ...){
    log.integral <- -sqrt(.Machine$double.xmax)
    inicial <- rep(0,n.dim)
    if(method.integrate == "BFGS"){
        temp <- try(optim(inicial,Q.b, gr = gr, ...,
                          method = otimizador, hessian = TRUE,
                          control = list(fnscale = -1)), silent=TRUE)
    }
    if(method.integrate == "NR"){
        temp <- try(newton.raphson(initial = inicial, escore = gr,
                                   hessiano = hess, n.dim = n.dim,
                                   max.iter = 100, ...), silent=TRUE)
    }
    if(method.integrate == "QNR"){
        temp <- try(qq.newton.raphson(initial = inicial, escore = gr,
                                      hessiano = hess, n.dim = n.dim,
                                      max.iter = 100, ...), silent=TRUE)
    }
    if(class(temp) != "try-error" & method.integrate == "BFGS"){
        log.integral <- temp$value + ((n.dim/2)*log(2 * pi) -
                                      0.5*determinant(-temp$hessian)$modulus)
    }
    if(class(temp) != "try-error" & method.integrate == "NR"){
        value <- Q.b(b = temp[2][[1]], ...)
        log.integral <- value + ((n.dim/2) * log(2 * pi) -
                                 0.5*determinant(-temp[1][[1]])$modulus)
    }
    .preditos <<- new.env(parent = .GlobalEnv)
    if(method.integrate == "NR"){
        assign("pred", value = temp[2][[1]], envir = .preditos)}
    if(method.integrate == "BFGS"){
        assign("pred", value = temp$par, envir = .preditos)}
    return(log.integral)
}

laplace <- function(Q.b, gr, hess, otimizador, n.dim, method.integrate, ...) {
    log.integral <- -sqrt(.Machine$double.xmax)
    inicial <- rep(0, n.dim)
    pred <- NULL
    if (method.integrate == "BFGS") {
        temp <- try(optim(inicial, Q.b, gr = gr, ..., method = otimizador, hessian = TRUE,
                          control = list(fnscale = -1)), silent = TRUE)
    } else if (method.integrate == "NR") {
        temp <- try(newton.raphson(initial = inicial, escore = gr, hessiano = hess,
                                   n.dim = n.dim, max.iter = 100, ...), silent = TRUE)
    } else if (method.integrate == "QNR") {
        temp <- try(qq.newton.raphson(initial = inicial, escore = gr, hessiano = hess,
                                      n.dim = n.dim, max.iter = 100, ...), silent = TRUE)
    }
    if (class(temp) != "try-error" && method.integrate == "BFGS") {
        log.integral <- temp$value + ((n.dim / 2) * log(2 * pi) - 0.5 * determinant(-temp$hessian)$modulus)
        pred <- temp$par
    } else if (class(temp) != "try-error" && method.integrate == "NR") {
        value <- Q.b(b = temp$solution, ...)
        log.integral <- value + ((n.dim / 2) * log(2 * pi) - 0.5 * determinant(-temp[1][[1]])$modulus)
        pred <- temp$solution
    }
    return(list(log.integral = log.integral, pred = pred))
}

loglik.sglmm_old <- function(par, Y, X, kappa, nugget, mat.dist, cov.model,
                             family, method.integrate = "NR", ntrial = 1, offset = NA){
    I = -sqrt(.Machine$double.xmax)
    n <- length(Y)
    n.beta <- dim(X)[2]
    beta <- as.numeric(par[1:n.beta])
    Xbeta <- X %*% beta
    if(is.na(offset)[1] != TRUE){
        Xbeta <- cbind(X, log(offset)) %*% c(beta, 1)
    }
    sigma <- exp(as.numeric(par[c(n.beta + 1)]))
    phi <- exp(as.numeric(par[c(n.beta + 2)]))
    if(nugget == TRUE){
        tau2 <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(nugget == FALSE){
        tau2 <- 0
    }
    if(family == "negative.binomial"   & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "gamma" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "beta" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "negative.binomial" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "gamma" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "beta" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if (!is.null(kappa)) {
        kappa = exp(as.numeric(kappa))
    }
    Sigma <- as.matrix(forceSymmetric(
        monta.sigma(cov.pars = c(sigma,phi),
                    cov.model = cov.model, nugget = tau2,
                    kappa = kappa, mat.dist = mat.dist)$varcov)
        )
    chol.Sigma <- try(chol(Sigma),silent=TRUE)
    det.Sigma <- try(sum(log(diag(chol.Sigma)))*2,silent=TRUE)
    inv.Sigma <- try(chol2inv(chol.Sigma), silent=TRUE)
    if(class(chol.Sigma)[1] != "try-error"){
        if(class(inv.Sigma)[1] != "try-error"){
            I <- laplace(Q.b, gr = Q.b.grad, hess = Q.b.hess,
                         method.integrate = method.integrate,
                         otimizador="BFGS", n.dim = n, Xbeta = Xbeta,
                         Y = Y, det.Sigma = det.Sigma,
                         inv.Sigma = inv.Sigma, family = family)
        }
    }
    return(-I)
}

loglik.sglmm <- function(par, Y, X, kappa, nugget, mat.dist, cov.model,
                         family, method.integrate = "NR", ntrial = 1, offset = NA){
    I = -sqrt(.Machine$double.xmax)
    n <- length(Y)
    n.beta <- dim(X)[2]
    beta <- as.numeric(par[1:n.beta])
    Xbeta <- X %*% beta
    if(is.na(offset)[1] != TRUE){
        Xbeta <- cbind(X, log(offset)) %*% c(beta, 1)
    }
    sigma <- exp(as.numeric(par[c(n.beta + 1)]))
    phi <- exp(as.numeric(par[c(n.beta + 2)]))
    if(nugget == TRUE){
        tau2 <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(nugget == FALSE){
        tau2 <- 0
    }
    if(family == "negative.binomial"   & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "gamma" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "beta" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "negative.binomial" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "gamma" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "beta" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if (!is.null(kappa)) {
        kappa = as.numeric(kappa)
    }
    Sigma <- as.matrix(forceSymmetric(
        monta.sigma(cov.pars = c(sigma,phi),
                    cov.model = cov.model, nugget = tau2,
                    kappa = kappa, mat.dist = mat.dist)$varcov)
        )
    chol.Sigma <- try(chol(Sigma),silent=TRUE)
    det.Sigma <- try(sum(log(diag(chol.Sigma)))*2,silent=TRUE)
    inv.Sigma <- try(chol2inv(chol.Sigma), silent=TRUE)
    if(class(chol.Sigma)[1] != "try-error"){
        if(class(inv.Sigma)[1] != "try-error"){
            if (any(family == c("poisson", "binomial", "geometric"))) {
                I <- laplace(Q.b, gr = Q.b.grad, hess = Q.b.hess,
                             method.integrate = method.integrate,
                             otimizador="BFGS", n.dim = n, Xbeta = Xbeta,
                             Y = Y, det.Sigma = det.Sigma,
                             inv.Sigma = inv.Sigma, family = family)
            }
        }
    }
    if(class(chol.Sigma)[1] != "try-error"){
        if(class(inv.Sigma)[1] != "try-error"){
            if (any(family == c("negative.binomial", "gamma", "beta"))) {
                I <- laplace(Q.b, gr = Q.b.grad, hess = Q.b.hess,
                             method.integrate = method.integrate,
                             otimizador="BFGS", n.dim = n, Xbeta = Xbeta,
                             Y = Y, det.Sigma = det.Sigma,, prec = prec,
                             inv.Sigma = inv.Sigma, family = family)
            }
        }
    }
    return(-I[[1]])
}

preditos <- function(par, Y, X, kappa, nugget, mat.dist, cov.model,
                     family, method.integrate = "NR", ntrial = 1, offset = NA){
    I = -sqrt(.Machine$double.xmax)
    n <- length(Y)
    n.beta <- dim(X)[2]
    beta <- as.numeric(par[1:n.beta])
    Xbeta <- X %*% beta
    if(is.na(offset)[1] != TRUE){
        Xbeta <- cbind(X, log(offset)) %*% c(beta, 1)
    }
    sigma <- exp(as.numeric(par[c(n.beta + 1)]))
    phi <- exp(as.numeric(par[c(n.beta + 2)]))
    if(nugget == TRUE){
        tau2 <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(nugget == FALSE){
        tau2 <- 0
    }
    if(family == "negative.binomial"   & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "gamma" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "beta" & nugget == TRUE){
        prec <- exp(as.numeric(par[c(n.beta+4)]))
    }
    if(family == "negative.binomial" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "gamma" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if(family == "beta" & nugget == FALSE){
        prec <- exp(as.numeric(par[c(n.beta+3)]))
    }
    if (!is.null(kappa)) {
        kappa = as.numeric(kappa)
    }
    Sigma <- as.matrix(forceSymmetric(
        monta.sigma(cov.pars = c(sigma,phi),
                    cov.model = cov.model, nugget = tau2,
                    kappa = kappa, mat.dist = mat.dist)$varcov)
        )
    chol.Sigma <- try(chol(Sigma),silent=TRUE)
    det.Sigma <- try(sum(log(diag(chol.Sigma)))*2,silent=TRUE)
    inv.Sigma <- try(chol2inv(chol.Sigma), silent=TRUE)
    if(class(chol.Sigma)[1] != "try-error"){
        if(class(inv.Sigma)[1] != "try-error"){
            if (any(family == c("poisson", "binomial", "geometric"))) {
                I <- laplace(Q.b, gr = Q.b.grad, hess = Q.b.hess,
                             method.integrate = method.integrate,
                             otimizador="BFGS", n.dim = n, Xbeta = Xbeta,
                             Y = Y, det.Sigma = det.Sigma,
                             inv.Sigma = inv.Sigma, family = family)
            }
        }
    }
    if(class(chol.Sigma)[1] != "try-error"){
        if(class(inv.Sigma)[1] != "try-error"){
            if (any(family == c("negative.binomial", "gamma", "beta"))) {
                I <- laplace(Q.b, gr = Q.b.grad, hess = Q.b.hess,
                             method.integrate = method.integrate,
                             otimizador="BFGS", n.dim = n, Xbeta = Xbeta,
                             Y = Y, det.Sigma = det.Sigma,, prec = prec,
                             inv.Sigma = inv.Sigma, family = family)
            }
        }
    }
    return(I$pred)
}

rglgm <- function(n.sample, family, beta, extra.prec = NA, X ,
                  cov.pars, nugget, kappa, cov.model, ntrial = 1, offset = 1,
                  xlims = c(0, 1), ylims = c(0, 1)){
    s <- grf(n = n.sample, cov.pars = cov.pars, cov.model = cov.model,
             kappa = kappa, nugget = nugget, messages = FALSE,
             xlims = xlims, ylims = ylims)
    Xbeta <- X%*%beta
    if(family == "binomial"){
        p <- inv.logit(Xbeta + s$data)
        y <- rbinom(n.sample, size = ntrial, prob = p)}
    if(family == "poisson"){
        lambda <- offset*exp(Xbeta + s$data)
        y <- rpois(n.sample, lambda = lambda)}
    if(family == "negative.binomial"){
        mu <- offset*exp(Xbeta + s$data)
        y <- rnbinom(n.sample,mu=mu,size=extra.prec)}
    if(family == "gamma"){
        mu <- exp(Xbeta + s$data)
        y <- rgamma(n.sample, shape = extra.prec, scale = mu/extra.prec)}
    if(family == "beta"){
        mu <- inv.logit(Xbeta + s$data)
        y <- rbeta(n.sample, mu*extra.prec, (1-mu)*extra.prec)}
    saida <- data.frame("y"=y, "coord.X" = s$coords[,1], "coord.Y" = s$coords[,2], "efeito" = s$data)
    return(saida)
}

start.values.sglmm <- function(formula, data, coords, nugget,
                               family, ntrial = 1, offset = NULL){
    if (is.null(offset)) {
        offset <- rep(1, nrow(data))
    }
    mf <- model.frame(formula,data)
    Y <- model.response(mf)
    X <- model.matrix(formula ,data=data)
    if( family == "binomial"){
        response <- cbind(Y,ntrial -Y)
        fit <- glm(response ~ -1 + X, data = data, family = "binomial")
        print(logLik(fit))
        esp <- predict(fit, type = "response")
        res <<- Y/ntrial - esp
        sigma = sd(Y/ntrial - esp)
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit), log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida,"logtau" = log(nugget))
        }
    }
    if(family == "poisson"){
        fit <- glm(Y ~ -1 + X, family= "poisson",
                   data = data, offset = log(offset))
        print(logLik(fit))
        esp <- predict(fit)
        sigma <- var(esp - log(Y + 1))
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit), log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida, "logtau" = log(nugget))
        }
    }
    if(family == "negative.binomial"){
        fit <- glm.nb(Y ~ -1 + X + offset(log(offset)), data = data)
        print(logLik(fit))
        esp <- predict(fit)
        sigma <- var( esp - log(Y + 1))
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit), log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida,"logtau" = log(nugget))
        }
        saida <- c(saida, "logtheta" = log(fit$theta))
    }
    if(family == "gamma"){
        fit <- glm(Y ~ -1 + X, family = Gamma(link = "log"),
                   data = data)
        te <- summary(fit)
        print(logLik(fit))
        esp <- predict(fit)
        sigma <- var(esp - log(Y))
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit), log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida, "logtau" = log(nugget))
        }
        saida <- c(saida, "logtheta" = te$dispersion)
    }
    if(family == "beta"){
        fit <- betareg(Y ~ -1 + X, data = data)
        te <- summary(fit)
        n.beta <- dim(X)[2]
        print(logLik(fit))
        esp <- predict(fit)
        sigma <- var(esp - Y)
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit)[1:n.beta], log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida, "logtau" = log(nugget))
        }
        saida <- c(saida, "logtheta" = log(as.numeric(coef(fit)[n.beta+1])))
    }
    if(family == "geometric"){
        fit <- glm(Y ~ -1 + X + offset(log(offset)), data = data,
                   family = negative.binomial(theta = 1))
        print(logLik(fit))
        esp <- predict(fit)
        sigma <- var(esp - log(Y + 1))
        phi <- 0.1 * max(dist(coords))
        saida <- c(coef(fit), log(sigma), log(phi))
        names(saida) <- c(colnames(X), "logsigma", "logphi")
        if(nugget == TRUE){
            nugget = 0.1 * sigma
            saida <- c(saida,"logtau" = log(nugget))
        }
    }
    return(saida)
}

sglmm_old <- function(formula, cov.model, kappa, inits, data,
                      coords, nugget, family, ntrial = 1, offset = 1,
                      method.optim = "BFGS", method.integrate = "NR",
                      predict = TRUE){
    formula <- as.formula(formula)
    mf <- model.frame(formula,data = data)
    Y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    mat.dist <- dist(coords)
    names <- c(colnames(X), "logsigma2", "logphi")
    n.beta <- dim(X)[2]
    if(family == "negative.binomial" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "gamma" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "beta" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "negative.binomial" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(family == "gamma" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(family == "beta" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(nugget == TRUE & family == "poisson"){
        names <- c(names,"logtau2")
    }
    if(nugget == TRUE & family == "binomial"){
        names <- c(names,"logtau2")
    }
    parnames(loglik.sglmm) <- names
    if (is.null(inits)) {
        inits <- start.values.glgm(formula, data, coords, nugget,
                                   family, ntrial, offset = offset)
    }
    names(inits) <- names
    estimativas <- mle2(loglik.sglmm, start = inits,
                        vecpar = TRUE,
                        method = method.optim,
                        control = list(maxit=1000),
                        skip.hessian = FALSE,
                        data = list(Y = Y, X = X, mat.dist = mat.dist,
                                    cov.model = cov.model,
                                    nugget = nugget, ntrial = ntrial,
                                    family = family, kappa =  kappa,
                                    method.integrate = method.integrate,
                                    offset = offset))
    n.pars <- length(coef(estimativas))
    summary.estimativas <- summary(estimativas)
    summary.estimativas@coef[,1][c(n.beta + 1):n.pars] <-
        exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])
    std.error = sqrt(exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])^2 *
                     (summary.estimativas@coef[,2][c(n.beta+1):n.pars]^2))
    summary.estimativas@coef[,2] <- c(summary.estimativas@coef[,2][1:n.beta],
                                      std.error)
    summary.estimativas@coef[,3] <- summary.estimativas@coef[,1]/
        summary.estimativas@coef[,2]
    summary.estimativas@coef[,4] <- NA
    if(predict == TRUE) {
        preditos <- .preditos$pred
        saida <- list()
        saida[1][[1]] <- summary.estimativas
        saida[7][[1]] <- logLik(estimativas)
        saida[2][[1]] <- preditos
        saida[3][[1]] <- coords
        saida[4][[1]] <- cov.model
        saida[5][[1]] <- family
        saida[6][[1]] <- exp(coef(estimativas)[c(n.beta + 1):n.pars])
        saida[8][[1]] <- ifelse(cov.model == "matern", exp(kappa), "NULL")
        saida[9][[1]] <- estimativas
        return(saida)}
    if(predict == FALSE){
        return(estimativas)
    }
}

sglmm <- function(formula, cov.model, kappa, inits, data,
                  coords, nugget, family, ntrial = 1, offset = 1,
                  method.optim = "BFGS", method.integrate = "NR",
                  predict = TRUE){
    formula <- as.formula(formula)
    mf <- model.frame(formula,data = data)
    Y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    mat.dist <- dist(coords)
    names <- c(colnames(X), "logsigma2", "logphi")
    n.beta <- dim(X)[2]
    if(family == "negative.binomial" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "gamma" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "beta" & nugget == TRUE){
        names <- c(names,"logtau2","logprec")
    }
    if(family == "negative.binomial" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(family == "gamma" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(family == "beta" & nugget == FALSE){
        names <- c(names,"logprec")
    }
    if(nugget == TRUE & family == "poisson"){
        names <- c(names,"logtau2")
    }
    if(nugget == TRUE & family == "binomial"){
        names <- c(names,"logtau2")
    }
    if(nugget == TRUE & family == "geometric"){
        names <- c(names,"logtau2")
    }
    parnames(loglik.sglmm) <- names
    if (is.null(inits)) {
        inits <- start.values.glgm(formula, data, coords, nugget,
                                   family, ntrial, offset = offset)
    }
    names(inits) <- names
    estimativas <- mle2(loglik.sglmm, start = inits,
                        vecpar = TRUE,
                        method = method.optim,
                        control = list(maxit=1000),
                        skip.hessian = FALSE,
                        data = list(Y = Y, X = X, mat.dist = mat.dist,
                                    cov.model = cov.model,
                                    nugget = nugget, ntrial = ntrial,
                                    family = family, kappa =  kappa,
                                    method.integrate = method.integrate,
                                    offset = offset))
    preditos <- preditos(par = coef(estimativas), Y = Y,
                         X = X, kappa = kappa, nugget = nugget,
                         mat.dist = mat.dist, cov.model = cov.model,
                         family = family, method.integrate = method.integrate,
                         offset = offset)
    n.pars <- length(coef(estimativas))
    summary.estimativas <- summary(estimativas)
    summary.estimativas@coef[,1][c(n.beta + 1):n.pars] <-
        exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])
    std.error = sqrt(exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])^2 *
                     (summary.estimativas@coef[,2][c(n.beta+1):n.pars]^2))
    summary.estimativas@coef[,2] <- c(summary.estimativas@coef[,2][1:n.beta],
                                      std.error)
    summary.estimativas@coef[,3] <- summary.estimativas@coef[,1]/
        summary.estimativas@coef[,2]
    summary.estimativas@coef[,4] <- NA
    if(predict == TRUE) {
        saida <- list()
        saida[1][[1]] <- summary.estimativas
        saida[7][[1]] <- logLik(estimativas)
        saida[2][[1]] <- preditos
        saida[3][[1]] <- coords
        saida[4][[1]] <- cov.model
        saida[5][[1]] <- family
        saida[6][[1]] <- exp(coef(estimativas)[c(n.beta + 1):n.pars])
        saida[8][[1]] <- ifelse(cov.model == "matern", kappa, "NULL")
        saida[9][[1]] <- estimativas
        return(saida)}
    if(predict == FALSE){
        return(estimativas)
    }
}

# Simulando de uma Poisson:

set.seed(5050)
sim.g <- grf(n = 200, grid = "irreg", cov.pars = c(0.5, 50),
             mean = 2, nugget = 0.05, xlims = c(0, 200),
             ylims = c(0, 200), cov.model = "matern",
             kappa = 2)

sim <- list(coords=sim.g$coords)
attr(sim,"class") <- "geodata"
sim$data <- rpois(200, lambda = exp(sim.g$data))

# Gráfico dos dados:

plot(sim)

# Transformando em data.frame:
sim_df <- data.frame(x1 = sim$coords[, 1],
                     x2 = sim$coords[, 2],
                     y = sim$data)

# Valores iniciais para theta:
inicial_simul <- start.values.sglmm(y ~ 1, family="poisson",
                                    data = sim_df, coords = sim$coords,
                                    nugget = T, offset = rep(1, length(sim$data)))

# Ajuste do modelo:
fit_sglmm <- sglmm(y ~ 1, cov.model = "matern", kappa = 2,
                   inits = inicial_simul, data = sim_df, coords = sim$coords,
                   nugget = T, family = "poisson", offset = rep(1, nrow(sim_df)),
                   method.optim = "BFGS", method.integrate = "NR")

# Estimativas pontuais:
fit_sglmm[[9]]

# Perfil de verossimilhança:
perfil_beta0 <- profile(fit_sglmm[[9]], which = 1)
perfil_sigma2 <- profile(fit_sglmm[[9]], which = 2)
perfil_phi <- profile(fit_sglmm[[9]], which = 3)
perfil_tau2 <- profile(fit_sglmm[[9]], which = 4)

par(mfrow = c(2, 2))

plot(perfil_beta0)
plot(perfil_sigma2)
plot(perfil_phi)
plot(perfil_tau2)

#----Ajuste para Weed----

load("~/Documentos/Projetos/LabB/Dados/geoCount/data/Weed.RData")

names(Weed) <- c("x1", "x2", "y", "Estimado")

png(filename = "~/Documentos/Projetos/LabB/Imagens/Plot_Geo_Weed.png")
plot(as.geodata(Weed, coords.col = 1:2, data.col = 3))
dev.off()

weed_geo <- as.geodata(Weed, coords.col = 1:2,
                       data.col = 3)

bordas <- Weed[chull(Weed$x1, Weed$x2), 1:2]

weed_geo$borders <- bordas

plot(weed_geo)

#----Modelos somente com intercepto----

# Ajustes com efeito de pepita:

# Valores iniciais:
theta_ini <- start.values.sglmm(y ~ 1, data = Weed, family = "poisson",
                                coords = Weed[, 1:2], nugget = T,
                                offset = rep(1, nrow(Weed)))

# Modelo matern com k = 1
fit_weed_1 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                    family = "poisson", inits = theta_ini, nugget = T,
                    cov.model = "matern", kappa = 1,
                    offset = rep(1, nrow(Weed)))

# Modelo matern com k = 1.5
fit_weed_1.5 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                      family = "poisson", inits = theta_ini, nugget = T,
                      cov.model = "matern", kappa = 1.5,
                      offset = rep(1, nrow(Weed)))

# Modelo matern com k = 2
fit_weed_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                    family = "poisson", inits = theta_ini, nugget = T,
                    cov.model = "matern", kappa = 2,
                    offset = rep(1, nrow(Weed)))

# Modelo exponencial
fit_weed_exp <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                      family = "poisson", inits = theta_ini, nugget = T,
                      cov.model = "exponential", kappa = NULL,
                      offset = rep(1, nrow(Weed)), method.integrate = "BFGS",
                      method.optim = "Nelder-Mead")

# Modelo esférico:
fit_weed_esf <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                      family = "poisson", inits = theta_ini, nugget = T,
                      cov.model = "spherical", kappa = NULL,
                      offset = rep(1, nrow(Weed)), method.integrate = "BFGS")

# Ajustes sem efeito de pepita:

# Valores iniciais:
theta_ini_2 <- start.values.sglmm(y ~ 1, data = Weed, family = "poisson",
                                  coords = Weed[, 1:2], nugget = F,
                                  offset = rep(1, nrow(Weed)))

# Modelo matern com k = 1
fit_weed_1_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                      family = "poisson", inits = theta_ini_2, nugget = F,
                      cov.model = "matern", kappa = 1,
                      offset = rep(1, nrow(Weed)))

# Modelo matern com k = 1.5
fit_weed_1.5_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                        family = "poisson", inits = theta_ini_2, nugget = F,
                        cov.model = "matern", kappa = 1.5,
                        offset = rep(1, nrow(Weed)))

# Modelo matern com k = 2
fit_weed_2_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                      family = "poisson", inits = theta_ini_2, nugget = F,
                      cov.model = "matern", kappa = 2,
                      offset = rep(1, nrow(Weed)))

# Modelo exponencial
fit_weed_exp_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                        family = "poisson", inits = theta_ini_2, nugget = F,
                        cov.model = "exponential", kappa = NULL,
                        offset = rep(1, nrow(Weed)), method.integrate = "BFGS",
                        method.optim = "Nelder-Mead")

# Modelo esférico:
fit_weed_esf_2 <- sglmm(y ~ 1, data = Weed, coords = Weed[, 1:2],
                        family = "poisson", inits = theta_ini_2, nugget = F,
                        cov.model = "spherical", kappa = NULL,
                        offset = rep(1, nrow(Weed)), method.integrate = "BFGS")

# Verificando as estimativas dos parâmetros e a máxima
# verossimilhança do modelo:

parametros <- cbind(
    data.frame(Modelo = c("$\\text{Exponencial} + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1.5) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 2) + \\tau^{2}$",
                          "$\\text{Esférico} + \\tau^{2}$",
                          "$\\text{Exponencial}$",
                          "$\\text{Matèrn}(\\kappa = 1)$",
                          "$\\text{Matèrn}(\\kappa = 1.5)$",
                          "$\\text{Matèrn}(\\kappa = 2)$",
                          "$\\text{Esférico}$")),
    rbind(
        do.call(rbind,
                lapply(list(fit_weed_exp[[9]],
                            fit_weed_1[[9]],
                            fit_weed_1.5[[9]],
                            fit_weed_2[[9]],
                            fit_weed_esf[[9]]),
                       coef)
                ),
        cbind(
            do.call(rbind,
                    lapply(list(fit_weed_exp_2[[9]],
                                fit_weed_1_2[[9]],
                                fit_weed_1.5_2[[9]],
                                fit_weed_2_2[[9]],
                                fit_weed_esf_2[[9]]),
                           coef)
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    )[, 1],
    rbind(
        do.call(rbind,
                list(fit_weed_exp[[6]],
                     fit_weed_1[[6]],
                     fit_weed_1.5[[6]],
                     fit_weed_2[[6]],
                     fit_weed_esf[[6]])
                ),
        cbind(
            do.call(rbind,
                    list(fit_weed_exp_2[[6]],
                         fit_weed_1_2[[6]],
                         fit_weed_1.5_2[[6]],
                         fit_weed_2_2[[6]],
                         fit_weed_esf_2[[6]]),
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    ),
    data.frame(LogLik = as.numeric(c(fit_weed_exp[[7]],
                                     fit_weed_1[[7]],
                                     fit_weed_1.5[[7]],
                                     fit_weed_2[[7]],
                                     fit_weed_esf[[7]],
                                     fit_weed_exp_2[[7]],
                                     fit_weed_1_2[[7]],
                                     fit_weed_1.5_2[[7]],
                                     fit_weed_2_2[[7]],
                                     fit_weed_esf_2[[7]])))
)

colnames(parametros) <- c("Modelo", "$\\hat{\\beta_{0}}$", "$\\hat{\\sigma^{2}}$",
                          "$\\hat{\\phi}$", "$\\hat{\\tau^{2}}$",
                          "logLik")

parametros[, -1] <- round(parametros[, -1], 4)

parametros[parametros == 0] <- "-"

mat <- xtable(parametros, auto = T, caption = "Estimativas dos parâmetros",
              digits = 4, label = "tab:parmodelos", align = "l|l|rrrrr|")

print(mat, sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      include.rownames = F)

anova(fit_weed_1[[9]], fit_weed_1_2[[9]])

a_1 <- anova(fit_weed_1[[9]], fit_weed_1_2[[9]])

print(xtable(matrix(a_1, nrow = 2, dimnames = dimnames(a_1)),
             auto = T, digits = 2))

perfil_beta0_1 <- profile(fit_weed_1_2[[9]], which = 1)
perfil_sigma_1 <- profile(fit_weed_1_2[[9]], which = 2)
perfil_phi_1 <- profile(fit_weed_1_2[[9]], which = 3)

par(mfrow = c(2, 2))

plot(perfil_beta0_1)
plot(perfil_sigma_1)
plot(perfil_phi_1)

## Modelos Binomial Negativa

theta_ini_bn <- start.values.sglmm(y ~ 1, data = Weed, family = "negative.binomial",
                                   coords = Weed[, 1:2], nugget = T)

theta_ini_2_bn <- start.values.sglmm(y ~ 1, data = Weed, family = "negative.binomial",
                                     coords = Weed[, 1:2], nugget = F)

fit_weed_exp_nb <- sglmm(y ~ 1, data = Weed, cov.model = "exponential", kappa = NULL,
                         coords = Weed[, 1:2], inits = theta_ini_bn, nugget = T,
                         family = "negative.binomial")

fit_weed_1_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 1,
                       coords = Weed[, 1:2], inits = theta_ini_bn, nugget = T,
                       family = "negative.binomial")

fit_weed_1.5_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 1.5,
                         coords = Weed[, 1:2], inits = theta_ini_bn, nugget = T,
                         family = "negative.binomial")

fit_weed_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 2,
                       coords = Weed[, 1:2], inits = theta_ini_bn, nugget = T,
                       family = "negative.binomial")

fit_weed_esf_nb <- sglmm(y ~ 1, data = Weed, cov.model = "spherical", kappa = NULL,
                         coords = Weed[, 1:2], inits = theta_ini_bn, nugget = T,
                         family = "negative.binomial")

fit_weed_exp_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "exponential", kappa = NULL,
                         coords = Weed[, 1:2], inits = theta_ini_2_bn, nugget = F,
                         family = "negative.binomial")

fit_weed_1_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 1,
                         coords = Weed[, 1:2], inits = theta_ini_2_bn, nugget = F,
                         family = "negative.binomial")

fit_weed_1.5_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 1.5,
                           coords = Weed[, 1:2], inits = theta_ini_2_bn, nugget = F,
                           family = "negative.binomial")

fit_weed_2_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "matern", kappa = 2,
                         coords = Weed[, 1:2], inits = theta_ini_2_bn, nugget = F,
                         family = "negative.binomial")

fit_weed_esf_2_nb <- sglmm(y ~ 1, data = Weed, cov.model = "spherical", kappa = NULL,
                           coords = Weed[, 1:2], inits = theta_ini_2_bn, nugget = F,
                           family = "negative.binomial")

parametros_nb <- cbind(
    data.frame(Modelo = c("$\\text{Exponencial} + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1.5) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 2) + \\tau^{2}$",
                          "$\\text{Esférico} + \\tau^{2}$",
                          "$\\text{Exponencial}$",
                          "$\\text{Matèrn}(\\kappa = 1)$",
                          "$\\text{Matèrn}(\\kappa = 1.5)$",
                          "$\\text{Matèrn}(\\kappa = 2)$",
                          "$\\text{Esférico}$")),
    rbind(
        do.call(rbind,
                lapply(list(fit_weed_exp_nb[[9]],
                            fit_weed_1_nb[[9]],
                            fit_weed_1.5_nb[[9]],
                            fit_weed_2_nb[[9]],
                            fit_weed_esf_nb[[9]]),
                       coef)
                ),
        cbind(
            do.call(rbind,
                    lapply(list(fit_weed_exp_2_nb[[9]],
                                fit_weed_1_2_nb[[9]],
                                fit_weed_1.5_2_nb[[9]],
                                fit_weed_2_2_nb[[9]],
                                fit_weed_esf_2_nb[[9]]),
                           coef)
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    )[, 1],
    rbind(
        do.call(rbind,
                list(fit_weed_exp_nb[[6]],
                     fit_weed_1_nb[[6]],
                     fit_weed_1.5_nb[[6]],
                     fit_weed_2_nb[[6]],
                     fit_weed_esf_nb[[6]])
                ),
        cbind(
            do.call(rbind,
                    list(fit_weed_exp_2_nb[[6]],
                         fit_weed_1_2_nb[[6]],
                         fit_weed_1.5_2_nb[[6]],
                         fit_weed_2_2_nb[[6]],
                         fit_weed_esf_2_nb[[6]]),
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    ),
    data.frame(LogLik = as.numeric(c(fit_weed_exp_nb[[7]],
                                     fit_weed_1_nb[[7]],
                                     fit_weed_1.5_nb[[7]],
                                     fit_weed_2_nb[[7]],
                                     fit_weed_esf_nb[[7]],
                                     fit_weed_exp_2_nb[[7]],
                                     fit_weed_1_2_nb[[7]],
                                     fit_weed_1.5_2_nb[[7]],
                                     fit_weed_2_2_nb[[7]],
                                     fit_weed_esf_2_nb[[7]])))
)

colnames(parametros_nb) <- c("Modelo", "$\\hat{\\beta_{0}}$", "$\\hat{\\sigma^{2}}$",
                             "$\\hat{\\phi}$", "$\\hat{\\tau^{2}}$", "$\\hat{\\psi}$",
                             "logLik")

parametros_nb[, -1] <- round(parametros_nb[, -1], 4)

parametros_nb[parametros_nb == 0] <- "-"

mat_nb <- xtable(parametros_nb, auto = T, caption = "Estimativas dos parâmetros",
                 digits = 4, label = "tab:parmodelosnb", align = "l|l|rrrrrr|")

print(mat_nb, sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      include.rownames = F)

p_b_0 <- profile(fit_weed_esf_2_nb[[9]], which = 1)
p_s_0 <- profile(fit_weed_esf_2_nb[[9]], which = 2)
p_ph_0 <- profile(fit_weed_esf_2_nb[[9]], which = 3)
p_pr_0 <- profile(fit_weed_esf_2_nb[[9]], which = 4)

p_b_0@profile$`(Intercept)` <- p_b_0@profile$`(Intercept)`[-1,]

par(mfrow = c(2, 2))

plot(p_b_0, xlab = "Intercepto", main = expression(paste("Perfil de Verossimilhança : ", beta[0])))
plot(p_s_0, xlab = expression(log(sigma^2)), main = expression(paste("Perfil de Verossimilhança : ", log(sigma^2))))
plot(p_ph_0, xlab = expression(log(phi)), main = expression(paste("Perfil de Verossimilhança : ", log(phi))))
plot(p_pr_0, xlab = expression(log(psi)), main = expression(paste("Perfil de Verossimilhança : ", log(psi))), xlim = c(0, 10))


# #----Modelos incluindo coordenadas----
#
# # Ajustes com efeito de pepita:
#
# # Valores iniciais:
# theta_ini_b <- start.values.sglmm(y ~ x1 + x2, data = Weed, family = "poisson",
#                                   coords = Weed[, 1:2], nugget = T,
#                                   offset = rep(1, nrow(Weed)))
#
# # Modelo matern com k = 1
# fit_weed_1_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                       family = "poisson", inits = theta_ini_b, nugget = T,
#                       cov.model = "matern", kappa = 1,
#                       offset = rep(1, nrow(Weed)))
#
# # Modelo matern com k = 1.5
# fit_weed_1.5_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                         family = "poisson", inits = theta_ini_b, nugget = T,
#                         cov.model = "matern", kappa = 1.5,
#                         offset = rep(1, nrow(Weed)))
#
# # Modelo matern com k = 2
# fit_weed_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                       family = "poisson", inits = theta_ini_b, nugget = T,
#                       cov.model = "matern", kappa = 2,
#                       offset = rep(1, nrow(Weed)))
#
# # Modelo exponencial
# fit_weed_exp_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                         family = "poisson", inits = theta_ini_b, nugget = T,
#                         cov.model = "exponential", kappa = NULL,
#                         offset = rep(1, nrow(Weed)))
#
# # Modelo esférico:
# fit_weed_esf_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                         family = "poisson", inits = theta_ini_b, nugget = T,
#                         cov.model = "spherical", kappa = NULL,
#                         offset = rep(1, nrow(Weed)))
#
# # Ajustes sem efeito de pepita:
#
# # Valores iniciais:
# theta_ini_2_b <- start.values.sglmm(y ~ x1 + x2, data = Weed, family = "poisson",
#                                     coords = Weed[, 1:2], nugget = F,
#                                     offset = rep(1, nrow(Weed)))
#
#
# # Modelo matern com k = 1
# fit_weed_1_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                         family = "poisson", inits = theta_ini_2_b, nugget = F,
#                         cov.model = "matern", kappa = 1,
#                         offset = rep(1, nrow(Weed)))
#
# # Modelo matern com k = 1.5
# fit_weed_1.5_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                           family = "poisson", inits = theta_ini_2_b, nugget = F,
#                           cov.model = "matern", kappa = 1.5,
#                           offset = rep(1, nrow(Weed)))
#
# # Modelo matern com k = 2
# fit_weed_2_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                         family = "poisson", inits = theta_ini_2_b, nugget = F,
#                         cov.model = "matern", kappa = 2,
#                         offset = rep(1, nrow(Weed)))
#
# # Modelo exponencial
# fit_weed_exp_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                           family = "poisson", inits = theta_ini_2_b, nugget = F,
#                           cov.model = "exponential", kappa = NULL,
#                           offset = rep(1, nrow(Weed)))
#
# # Modelo esférico:
# fit_weed_esf_2_b <- sglmm(y ~ x1 + x2, data = Weed, coords = Weed[, 1:2],
#                           family = "poisson", inits = theta_ini_2_b, nugget = F,
#                           cov.model = "spherical", kappa = NULL,
#                           offset = rep(1, nrow(Weed)))
#
# # Verificando as estimativas dos parâmetros e a máxima
# # verossimilhança do modelo:
#
# parametros_b <- cbind(
#     data.frame(Modelo = c("$\\text{Exponencial} + \\tau^{2}$",
#                           "$\\text{Matèrn}(\\kappa = 1) + \\tau^{2}$",
#                           "$\\text{Matèrn}(\\kappa = 1.5) + \\tau^{2}$",
#                           "$\\text{Matèrn}(\\kappa = 2) + \\tau^{2}$",
#                           "$\\text{Esférico} + \\tau^{2}$",
#                           "$\\text{Exponencial}$",
#                           "$\\text{Matèrn}(\\kappa = 1)$",
#                           "$\\text{Matèrn}(\\kappa = 1.5)$",
#                           "$\\text{Matèrn}(\\kappa = 2)$",
#                           "$\\text{Esférico}$")),
#     rbind(
#         do.call(rbind,
#                 lapply(list(fit_weed_exp_b[[9]],
#                             fit_weed_1_b[[9]],
#                             fit_weed_1.5_b[[9]],
#                             fit_weed_2_b[[9]],
#                             fit_weed_esf_b[[9]]),
#                        coef)
#                 ),
#         cbind(
#             do.call(rbind,
#                     lapply(list(fit_weed_exp_2_b[[9]],
#                                 fit_weed_1_2_b[[9]],
#                                 fit_weed_1.5_2_b[[9]],
#                                 fit_weed_2_2_b[[9]],
#                                 fit_weed_esf_2_b[[9]]),
#                            coef)
#                     ),
#             data.frame(logtau2 = rep(0, 5))
#         )
#     ),
#     data.frame(LogLik = as.numeric(c(fit_weed_exp_b[[7]],
#                                      fit_weed_1_b[[7]],
#                                      fit_weed_1.5_b[[7]],
#                                      fit_weed_2_b[[7]],
#                                      fit_weed_esf_b[[7]],
#                                      fit_weed_exp_2_b[[7]],
#                                      fit_weed_1_2_b[[7]],
#                                      fit_weed_1.5_2_b[[7]],
#                                      fit_weed_2_2_b[[7]],
#                                      fit_weed_esf_2_b[[7]])))
# )
#
# colnames(parametros_b) <- c("Modelo", "$\\hat{\\beta_{0}}$",
#                             "$\\hat{\\beta_{1}}$", "$\\hat{\\beta_{2}}$",
#                             "$\\log(\\hat{\\sigma^{2}})$", "$\\log(\\hat{\\phi})$",
#                             "$\\log(\\hat{\\tau^{2}})$", "logLik")
#
#
# parametros_b[, -1] <- round(parametros_b[, -1], 4)
#
# parametros_b[parametros_b == 0] <- "-"
#
# mat_b <- xtable(parametros_b, auto = T, caption = "Estimativas dos parâmetros",
#                 digits = 4, label = "tab:parmodelosb")
#
# print(mat_b, sanitize.text.function = identity,
#       sanitize.colnames.function = identity,
#       include.rownames = F)
#
# parametros_b[which(parametros_b$logLik == max(parametros_b$logLik)), ]
#
# perfil_beta0_2 <- profile(fit_weed_1_2_b[[9]], which = 1)
# perfil_beta1_2 <- profile(fit_weed_1_2_b[[9]], which = 2)
# perfil_beta2_2 <- profile(fit_weed_1_2_b[[9]], which = 3)
# perfil_sigma_2 <- profile(fit_weed_1_2_b[[9]], which = 4)
# perfil_phi_2 <- profile(fit_weed_1_2_b[[9]], which = 5)
#
# perfil_beta0_2@profile$`(Intercept)` <- perfil_beta0_2@profile$`(Intercept)`[-1,]
# perfil_beta2_2@profile$x2 <- perfil_beta2_2@profile$x2[-1,]
#
# par(mfrow = c(3, 2))
# plot(perfil_beta0_2, xlim = c(3, 6.5))
# plot(perfil_beta1_2, xlim = c(-0.004, 0.005))
# plot(perfil_beta2_2, xlim = c(-0.006, 0.003))
# plot(perfil_sigma_2)
# plot(perfil_phi_2)

#----Ajuste usando MCMC----

beta <- theta_ini[1]
sigma2 <- exp(theta_ini[2])
phi <- exp(theta_ini[3])
tau2 <- exp(theta_ini[4])

mcmc.1 <- list(cov.pars = c(sigma2, phi, tau2), link = "log",
               beta = beta, family = "poisson", cov.model = "spherical")
S.prop <- mcmc.control(S.scale=0.1, thin=10)

tune.S <- glsm.mcmc(as.geodata(Weed, coords.col = 1:2,
                               data.col = 3),
                    model=mcmc.1, mcmc.input=S.prop)
S.control <- mcmc.control(S.scale=0.5, thin=50, burn.in=20000)

S.sims <- glsm.mcmc(sim, model = mcmc.1, mcmc.in = S.control)
lik.control <- prepare.likfit.glsm(S.sims)

mc.fit1 <- likfit.glsm(lik.control, ini.phi = phi, fix.nugget.rel = F)
mc.fit1

# Krigagem e predição espacial:

weed_pred <- data.frame(Pred = fit_weed_1_2[[2]])

weed_pred <- cbind(weed_pred, fit_weed_1_2[[3]])

weed_pred <- as.geodata(weed_pred, coords.col = 2:3, data.col = 1)

gr <- pred_grid(bordas, by = 7)
kc <- krige.control(beta = coef(fit_weed_1_2[[9]])[1],
                    cov.pars = fit_weed_1_2[[6]],
                    cov.model = "matern", kappa = 1)
oc <- output.control(n.predictive = 1000, simulations.predictive = T,
                     threshold = 250)
pred <- krige.conv(weed_pred, locations = gr, borders = bordas,
                   krige = kc, output = oc)

image(pred, col = hcl.colors(20, "blues", rev = T),
      x.leg = c(0, 200), y.leg = c(0, 30))

#----Ajuste para dados SPT----

spt <- readxl::read_xlsx("~/Downloads/SPT_PJ.xlsm")
bor <- read.delim("~/Downloads/bordersPJ.txt", sep = ",")
names(spt) <- c("Z", "x1", "x2", "y")

spt$y[spt$y == 99] <- NA
spt$Z <- -1 * as.numeric(spt$Z)
spt$x2 <- max(spt$x2) - spt$x2

spt_geo <- as.geodata(spt, coords.col = 2:3, data.col = 4)
spt_geo$borders <- bor

plot(spt_geo)

boxplot(y ~ Z, data = spt)

## Removendo NAs e observações truncadas em y == 50:
spt_4 <- subset(spt, Z == 4 & !is.na(y))# & y != 50)

spt_4_geo <- as.geodata(spt_4, coords.col = 2:3,
                        data.col = 4)

spt_13 <- subset(spt, Z == 13 & !is.na(y))# & y != 50)

spt_13_geo <- as.geodata(spt_13, coords.col = 2:3,
                         data.col = 4)

spt_4_geo$borders <- bor
spt_13_geo$borders <- bor

plot(spt_4_geo)
plot(spt_13_geo)

# Ajuste para profundidade 4 m:
inits_4 <- start.values.sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                              nugget = F, family = "poisson", offset = rep(1, length(spt_4$Z)))

inits_4_nug <- start.values.sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                                  nugget = T, family = "poisson", offset = rep(1, length(spt_4$Z)))

mod_exp_4 <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                   cov.model = "exponential", kappa = "NULL",
                   inits = inits_4, nugget = F, family = "poisson")

mod_exp_4_nug <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                       cov.model = "exponential", kappa = "NULL",
                       inits = inits_4_nug, nugget = T, family = "poisson")

mod_esf_4 <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                   cov.model = "spherical", kappa = NULL,
                   inits = inits_4, nugget = F, family = "poisson")

mod_esf_4_nug <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                       cov.model = "spherical", kappa = NULL,
                       inits = inits_4_nug, nugget = T, family = "poisson")

mod_1_4 <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                 cov.model = "matern", kappa = 1,
                 inits = inits_4, nugget = F, family = "poisson")

mod_1_4_nug <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                 cov.model = "matern", kappa = 1,
                 inits = inits_4_nug, nugget = T, family = "poisson")

mod_1.5_4 <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                   cov.model = "matern", kappa = 1.5,
                   inits = inits_4, nugget = F, family = "poisson")

mod_1.5_4_nug <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                       cov.model = "matern", kappa = 1.5,
                       inits = inits_4_nug, nugget = T, family = "poisson")

mod_2_4 <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                 cov.model = "matern", kappa = 2,
                 inits = inits_4, nugget = F, family = "poisson")

mod_2_4_nug <- sglmm(y ~ 1, data = spt_4, coords = spt_4[, 2:3],
                     cov.model = "matern", kappa = 2,
                     inits = inits_4_nug, nugget = T, family = "poisson")


list(mod_exp_4[[7]],
     mod_esf_4[[7]],
     mod_1_4[[7]],
     mod_1.5_4[[7]],
     mod_2_4[[7]])

parametros_spt_4 <- cbind(
    data.frame(Modelo = c("$\\text{Exponencial} + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1.5) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 2) + \\tau^{2}$",
                          "$\\text{Esférico} + \\tau^{2}$",
                          "$\\text{Exponencial}$",
                          "$\\text{Matèrn}(\\kappa = 1)$",
                          "$\\text{Matèrn}(\\kappa = 1.5)$",
                          "$\\text{Matèrn}(\\kappa = 2)$",
                          "$\\text{Esférico}$")),
    rbind(
        do.call(rbind,
                lapply(list(mod_exp_4_nug[[9]],
                            mod_1_4_nug[[9]],
                            mod_1.5_4_nug[[9]],
                            mod_2_4_nug[[9]],
                            mod_esf_4_nug[[9]]),
                       coef)
                ),
        cbind(
            do.call(rbind,
                    lapply(list(mod_exp_4[[9]],
                                mod_1_4[[9]],
                                mod_1.5_4[[9]],
                                mod_2_4[[9]],
                                mod_esf_4[[9]]),
                           coef)
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    )[, 1],
    rbind(
        do.call(rbind,
                list(mod_exp_4_nug[[6]],
                     mod_1_4_nug[[6]],
                     mod_1.5_4_nug[[6]],
                     mod_2_4_nug[[6]],
                     mod_esf_4_nug[[6]])
                ),
        cbind(
            do.call(rbind,
                    list(mod_exp_4[[6]],
                         mod_1_4[[6]],
                         mod_1.5_4[[6]],
                         mod_2_4[[6]],
                         mod_esf_4[[6]]),
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    ),
    data.frame(LogLik = as.numeric(c(mod_exp_4_nug[[7]],
                                     mod_1_4_nug[[7]],
                                     mod_1.5_4_nug[[7]],
                                     mod_2_4_nug[[7]],
                                     mod_esf_4_nug[[7]],
                                     mod_exp_4[[7]],
                                     mod_1_4[[7]],
                                     mod_1.5_4[[7]],
                                     mod_2_4[[7]],
                                     mod_esf_4[[7]])))
)

colnames(parametros_spt_4) <- c("Modelo", "$\\hat{\\beta_{0}}$",
                                "$\\hat{\\sigma^{2}}$", "$\\hat{\\phi}$",
                                "$\\hat{\\tau^{2}}$", "logLik")


mat_spt_4 <- xtable(parametros_spt_4, auto = T, caption = "Estimativas para profundidade 4 metros",
                    digits = 4, label = "tab:parmodelosspt4", align = "l|l|rrrrr|")

print(mat_spt_4, sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      include.rownames = F)

perf_b_1_4 <- profile(mod_2_4_nug[[9]], which = 1)
perf_s_1_4 <- profile(mod_2_4_nug[[9]], which = 2)
perf_p_1_4 <- profile(mod_2_4_nug[[9]], which = 3)
perf_t_1_4 <- profile(mod_2_4_nug[[9]], which = 4)

par(mfrow = c(2, 2))

plot(perf_b_1_4)
plot(perf_s_1_4, xlim = c(-10, 2))
plot(perf_p_1_4, xlim = c(-2, 5))
plot(perf_t_1_4, xlim = c(-3, -0.8))

# Ajuste para profundidade 13 m:

inits_13 <- start.values.sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                               nugget = F, family = "poisson", offset = rep(1, length(spt_13$Z)))

inits_13_nug <- start.values.sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                                   nugget = T, family = "poisson")

mod_exp <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                 cov.model = "exponential", kappa = "NULL",
                 inits = inits_13, nugget = F, family = "poisson")

mod_exp_nug <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                     cov.model = "exponential", kappa = "NULL",
                     inits = inits_13_nug, nugget = T, family = "poisson")

mod_esf <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                 cov.model = "spherical", kappa = NULL,
                 inits = inits_13, nugget = F, family = "poisson")

mod_esf_nug <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                     cov.model = "spherical", kappa = NULL,
                     inits = inits_13_nug, nugget = T, family = "poisson")

mod_1 <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
               cov.model = "matern", kappa = 1,
               inits = inits_13, nugget = F, family = "poisson")

mod_1_nug <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                   cov.model = "matern", kappa = 1,
                   inits = inits_13_nug, nugget = T, family = "poisson")

mod_1.5 <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                 cov.model = "matern", kappa = 1.5,
                 inits = inits_13, nugget = F, family = "poisson")

mod_1.5_nug <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                     cov.model = "matern", kappa = 1.5,
                     inits = inits_13_nug, nugget = T, family = "poisson")

mod_2 <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
               cov.model = "matern", kappa = 2,
               inits = inits_13, nugget = F, family = "poisson")

mod_2_nug <- sglmm(y ~ 1, data = spt_13, coords = spt_13[, 2:3],
                   cov.model = "matern", kappa = 2,
                   inits = inits_13_nug, nugget = T, family = "poisson")

list(mod_exp[[7]],
     mod_esf[[7]],
     mod_1[[7]],
     mod_1.5[[7]],
     mod_2[[7]])

parametros_spt_13 <- cbind(
    data.frame(Modelo = c("$\\text{Exponencial} + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 1.5) + \\tau^{2}$",
                          "$\\text{Matèrn}(\\kappa = 2) + \\tau^{2}$",
                          "$\\text{Esférico} + \\tau^{2}$",
                          "$\\text{Exponencial}$",
                          "$\\text{Matèrn}(\\kappa = 1)$",
                          "$\\text{Matèrn}(\\kappa = 1.5)$",
                          "$\\text{Matèrn}(\\kappa = 2)$",
                          "$\\text{Esférico}$")),
    rbind(
        do.call(rbind,
                lapply(list(mod_exp_nug[[9]],
                            mod_1_nug[[9]],
                            mod_1.5_nug[[9]],
                            mod_2_nug[[9]],
                            mod_esf_nug[[9]]),
                       coef)
                ),
        cbind(
            do.call(rbind,
                    lapply(list(mod_exp[[9]],
                                mod_1[[9]],
                                mod_1.5[[9]],
                                mod_2[[9]],
                                mod_esf[[9]]),
                           coef)
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    )[, 1],
    rbind(
        do.call(rbind,
                list(mod_exp_nug[[6]],
                     mod_1_nug[[6]],
                     mod_1.5_nug[[6]],
                     mod_2_nug[[6]],
                     mod_esf_nug[[6]])
                ),
        cbind(
            do.call(rbind,
                    list(mod_exp[[6]],
                         mod_1[[6]],
                         mod_1.5[[6]],
                         mod_2[[6]],
                         mod_esf[[6]]),
                    ),
            data.frame(logtau2 = rep(0, 5))
        )
    ),
    data.frame(LogLik = as.numeric(c(mod_exp_nug[[7]],
                                     mod_1_nug[[7]],
                                     mod_1.5_nug[[7]],
                                     mod_2_nug[[7]],
                                     mod_esf_nug[[7]],
                                     mod_exp[[7]],
                                     mod_1[[7]],
                                     mod_1.5[[7]],
                                     mod_2[[7]],
                                     mod_esf[[7]])))
)

colnames(parametros_spt_13) <- c("Modelo", "$\\hat{\\beta_{0}}$",
                                "$\\hat{\\sigma^{2}}$", "$\\hat{\\phi}$",
                                "$\\hat{\\tau^{2}}$", "logLik")

mat_spt_13 <- xtable(parametros_spt_13, auto = T, caption = "Estimativas para profundidade 13 metros",
                     digits = 4, label = "tab:parmodelosspt13", align = "l|l|rrrrr|")

print(mat_spt_13, sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      include.rownames = F)

perf_b_1 <- profile(mod_2[[9]], which = 1)
perf_s_1 <- profile(mod_2[[9]], which = 2)
perf_p_1 <- profile(mod_2[[9]], which = 3)

par(mfrow = c(2, 2))

plot(perf_b_1)
plot(perf_s_1)
plot(perf_p_1)

# Krigagem e predição espacial

spt_4_pred <- cbind(data.frame(Pred = mod_2_4_nug[[2]]),
                    mod_2_4_nug[[3]])

spt_4_pred <- as.geodata(spt_4_pred, data.col = 1, coords.col = 2:3)

gr_4 <- pred_grid(bor, by = 0.5)
kc_4 <- krige.control(beta = coef(mod_2_4_nug[[9]])[1],
                      cov.pars = mod_2_4_nug[[6]][1:2],
                      cov.model = "matern", kappa = 2,
                      nugget = mod_2_4_nug[[6]][3])
oc_4 <- output.control(n.predictive = 1000, simulations.predictive = T,
                       threshold = 250)
pred_4 <- krige.conv(spt_4_pred, locations = gr_4, borders = bor,
                     krige = kc_4, output = oc_4)

image(pred_4, col = hcl.colors(20, "blues", rev = T),
      x.leg = c(-15, 0), y.leg = c(-15, -10))

spt_13_pred <- cbind(data.frame(Pred = mod_2[[2]]),
                    mod_2[[3]])

spt_13_pred <- as.geodata(spt_13_pred, data.col = 1, coords.col = 2:3)

gr_13 <- pred_grid(bor, by = 0.5)
kc_13 <- krige.control(beta = coef(mod_2[[9]])[1],
                       cov.pars = mod_2[[6]],
                       cov.model = "matern", kappa = 2)
oc_13 <- output.control(n.predictive = 1000, simulations.predictive = T,
                        threshold = 250)
pred_13 <- krige.conv(spt_13_pred, locations = gr_13, borders = bor,
                      krige = kc_13, output = oc_13)

image(pred_13, col = hcl.colors(20, "blues", rev = T),
      x.leg = c(-15, 0), y.leg = c(-5, -3))

#----Estudo de simulação----

# Simulando de uma Poisson:

set.seed(5050)
sim.g <- grf(n = 200, grid = "irreg", cov.pars = c(0.5, 30),
             mean = 2, nugget = 0.05, xlims = c(0, 200),
             ylims = c(0, 200), cov.model = "matern",
             kappa = 2)

sim <- list(coords=sim.g$coords)
attr(sim,"class") <- "geodata"
sim$data <- rpois(200, lambda = exp(sim.g$data))

# Gráfico dos dados:

plot(sim)

# Transformando em data.frame:
sim_df <- data.frame(x1 = sim$coords[, 1],
                     x2 = sim$coords[, 2],
                     y = sim$data)

# Valores iniciais para theta:
inicial_simul <- start.values.sglmm(y ~ 1, family="poisson",
                                    data = sim_df, coords = sim$coords,
                                    nugget = T, offset = rep(1, length(sim$data)))

# Ajuste do modelo:
fit_sglmm <- sglmm(y ~ 1, cov.model = "matern", kappa = 2,
                   inits = inicial_simul, data = sim_df, coords = sim$coords,
                   nugget = T, family = "poisson", offset = rep(1, nrow(sim_df)),
                   method.optim = "BFGS", method.integrate = "NR")

# Estimativas pontuais:
fit_sglmm[[9]]

# Perfil de verossimilhança:
perfil_beta0 <- profile(fit_sglmm[[9]], which = 1)
perfil_sigma2 <- profile(fit_sglmm[[9]], which = 2)
perfil_phi <- profile(fit_sglmm[[9]], which = 3)
perfil_tau2 <- profile(fit_sglmm[[9]], which = 4)

par(mfrow = c(2, 2))

plot(perfil_beta0)
plot(perfil_sigma2)
plot(perfil_phi)
plot(perfil_tau2)

## Usando processamento em paralelo:

library(parallel)
library(doParallel)

# Definir funções de simulação:
simulation_function <- function(i, n = 100, cov.model = "exponential", kappa = "NULL", xlims = c(0, 200), ylims = c(0, 200)) {
  set.seed(i)
  sim.g <- grf(n = n, grid = "irreg", cov.pars = c(0.5, 30),
               mean = 2, nugget = 0.05, xlims = xlims,
               ylims = ylims, cov.model = "exponential")
  sim <- list(coords = sim.g$coords)
  attr(sim, "class") <- "geodata"
  sim$data <- rpois(n, lambda = exp(sim.g$data))
  sim_df <- data.frame(x1 = sim$coords[, 1],
                       x2 = sim$coords[, 2],
                       y = sim$data)
  val_ini_simu <- start.values.sglmm(y ~ 1, family = "poisson", data = sim_df,
                                     coords = sim_df[, 1:2], nugget = TRUE, offset = rep(1, nrow(sim_df)))
  ajuste_ini_simu <- sglmm(y ~ 1, cov.model = cov.model, kappa = kappa, inits = val_ini_simu, data = sim_df,
                           coords = sim_df[, 1:2], nugget = TRUE, family = "poisson", offset = rep(1, nrow(sim_df)),
                           method.optim = "BFGS", method.integrate = "NR")
  return(c(ajuste_ini_simu[[9]]@coef[1], ajuste_ini_simu[[6]]))
}

simulation_function_2 <- function(i, n = 100, cov.model = "exponential", kappa = "NULL", xlims = c(0, 200), ylims = c(0, 200)) {
  set.seed(i)
  sim.g <- grf(n = n, grid = "irreg", cov.pars = c(0.5, 30),
               mean = 2, nugget = 0.05, xlims = xlims,
               ylims = ylims, cov.model = "exponential")
  sim <- list(coords = sim.g$coords)
  attr(sim, "class") <- "geodata"
  sim$data <- rpois(n, lambda = exp(sim.g$data))
  sim_df <- data.frame(x1 = sim$coords[, 1],
                       x2 = sim$coords[, 2],
                       y = sim$data)
  val_ini_simu <- start.values.sglmm(y ~ 1, family = "poisson", data = sim_df,
                                     coords = sim_df[, 1:2], nugget = TRUE, offset = rep(1, nrow(sim_df)))
  mcmc <- list(cov.pars = exp(val_ini_simu[2:4]), link = "log",
               beta = val_ini_simu[1], family = "poisson", cov.model = cov.model, kappa = kappa)
  S.prop <- mcmc.control(S.scale=0.1, thin=10)
  tune.S <- glsm.mcmc(sim, model = mcmc,
                      mcmc.input = S.prop)
  S.control <- mcmc.control(S.scale = 0.5, thin = 50, burn.in = 10000)
  S.sims <- glsm.mcmc(sim, model = mcmc, mcmc.in = S.control)
  lik.control <- prepare.likfit.glsm(S.sims)
  mc.fit <- likfit.glsm(lik.control, ini.phi = exp(val_ini_simu[3]), fix.nugget.rel = F)
  return(c(mc.fit$beta, mc.fit$cov.pars, mc.fit$nugget.rel))
}

# Cria o cluster, com n-1 núcleos para processamento:
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

# Exporta as bibliotecas necessárias para o cluster:
clusterEvalQ(cl, {
    library(geoR)
    library(geoRglm)
    library(bbmle)
})

# Exporta as funções necessárias para o cluster:
clusterExport(cl, varlist = c("simulation_function",
                              "start.values.sglmm", "gauss.mult",
                              "sglmm","loglik.sglmm","parnames<-","mle2",
                              "forceSymmetric","monta.sigma", "laplace", "preditos",
                              "Q.b", "Q.b.grad", "Q.b.hess", "newton.raphson",
                              "simulation_function_2"))

# Ajustes para diferentes tamanhos de amostra:

resultados_50 <- parLapply(cl, 1:1000, simulation_function, n = 50)
resultados_100 <- parLapply(cl, 1:1000, simulation_function, n = 100)
resultados_200 <- parLapply(cl, 1:1000, simulation_function, n = 200)

resultados_50 <- do.call(rbind, resultados_50)
resultados_100 <- do.call(rbind, resultados_100)
resultados_200 <- do.call(rbind, resultados_200)

eqm_b0_50 <- mean((resultados_50[, 1] - 2)^2)
eqm_b0_100 <- mean((resultados_100[, 1] - 2)^2)
eqm_b0_200 <- mean((resultados_200[, 1] - 2)^2)

eqm_s_50 <- mean((resultados_50[, 2] - 0.5)^2)
eqm_s_100 <- mean((resultados_100[, 2] - 0.5)^2)
eqm_s_200 <- mean((resultados_200[, 2] - 0.5)^2)

eqm_p_50 <- mean((resultados_50[, 3] - 30)^2)
eqm_p_100 <- mean((resultados_100[, 3] - 30)^2)
eqm_p_200 <- mean((resultados_200[, 3] - 30)^2)

eqm_t_50 <- mean((resultados_50[, 4] - 0.05)^2)
eqm_t_100 <- mean((resultados_100[, 4] - 0.05)^2)
eqm_t_200 <- mean((resultados_200[, 4] - 0.05)^2)

media_n <- do.call(rbind, lapply(list(resultados_50,
                                      resultados_100,
                                      resultados_200),
                                 colMeans))


var_n <- do.call(rbind, lapply(list(resultados_50,
                                    resultados_100,
                                    resultados_200), FUN = function(x) {apply(x, 2, var)}))

eqm_n <- matrix(c(eqm_b0_50, eqm_b0_100, eqm_b0_200,
                   eqm_s_50, eqm_s_100, eqm_s_200,
                   eqm_p_50, eqm_p_100, eqm_p_200,
                   eqm_t_50, eqm_t_100, eqm_t_200), nrow = 3)

# Viés do estimador

vies_n <- eqm_n - var_n

rownames(media_n) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(var_n) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(eqm_n) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(vies_n) <- c("$n = 50$", "$n = 100$", "$n = 200$")

colnames(eqm_n) <- colnames(media_n) <- colnames(var_n) <- colnames(vies_n) <- c("$\\hat{\\beta_{0}}$","$\\hat{\\sigma^{2}}$",
                                                                                 "$\\hat{\\phi}$", "$\\hat{\\tau^{2}}$")

mat_media_n <- xtable(media_n, auto = T, caption = "Média dos estimadores",
                      digits = 4, label = "tab:media")
mat_var_n <- xtable(var_n, auto = T, caption = "Variância dos estimadores",
                    digits = 4, label = "tab:varestn")
mat_eqm_n <- xtable(eqm_n, auto = T, caption = "Erro Quadrático Médio dos estimadores",
                    digits = 4, label = "tab:eqmn")
mat_vies_n <- xtable(vies_n, auto = T, caption = "Viés dos estimadores",
                     digits = 4, label = "tab:viesparn")

print(
    xtable(
        data.frame(
            `Parâmetro` = c("$\\hat{\\beta_{0}}$", NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA,
                    "$\\hat{\\phi}$", NA, NA, "$\\hat{\\tau^{2}}$", NA, NA),
            n = c(rep(c("$n = 50$", "$n = 100$", "$n = 200$"), 4)),
            `Média` = c(media_n[, 1], media_n[, 2], media_n[, 3], media_n[, 4]),
            `Variância` = c(var_n[, 1], var_n[, 2], var_n[, 3], var_n[, 4]),
            `EQM` = c(eqm_n[, 1], eqm_n[, 2], eqm_n[, 3], eqm_n[, 4]),
            `Viés` = c(vies_n[, 1], vies_n[, 2], vies_n[, 3], vies_n[, 4])
        ),
        digits = 4,
        caption = "Estimadores para diferentes tamanhos de amostras",
        label = "tab:estamo"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F
)

# MCMC

resultados_50_mcmc <- parLapply(cl, 1:1000, simulation_function_2, n = 50)
resultados_100_mcmc <- parLapply(cl, 1:1000, simulation_function_2, n = 100)
resultados_200_mcmc <- parLapply(cl, 1:1000, simulation_function_2, n = 200)

resultados_50_mcmc <- do.call(rbind, resultados_50_mcmc)
resultados_100_mcmc <- do.call(rbind, resultados_100_mcmc)
resultados_200_mcmc <- do.call(rbind, resultados_200_mcmc)

eqm_b0_50_mcmc <- mean((resultados_50_mcmc[, 1] - 2)^2)
eqm_b0_100_mcmc <- mean((resultados_100_mcmc[, 1] - 2)^2)
eqm_b0_200_mcmc <- mean((resultados_200_mcmc[, 1] - 2)^2)

eqm_s_50_mcmc <- mean((resultados_50_mcmc[, 2] - 0.5)^2)
eqm_s_100_mcmc <- mean((resultados_100_mcmc[, 2] - 0.5)^2)
eqm_s_200_mcmc <- mean((resultados_200_mcmc[, 2] - 0.5)^2)

eqm_p_50_mcmc <- mean((resultados_50_mcmc[, 3] - 30)^2)
eqm_p_100_mcmc <- mean((resultados_100_mcmc[, 3] - 30)^2)
eqm_p_200_mcmc <- mean((resultados_200_mcmc[, 3] - 30)^2)

eqm_t_50_mcmc <- mean((resultados_50_mcmc[, 4] - 0.05)^2)
eqm_t_100_mcmc <- mean((resultados_100_mcmc[, 4] - 0.05)^2)
eqm_t_200_mcmc <- mean((resultados_200_mcmc[, 4] - 0.05)^2)

media_n_mcmc <- do.call(rbind, lapply(list(resultados_50_mcmc,
                                           resultados_100_mcmc,
                                           resultados_200_mcmc),
                                      colMeans))


var_n_mcmc <- do.call(rbind, lapply(list(resultados_50_mcmc,
                                         resultados_100_mcmc,
                                         resultados_200_mcmc), FUN = function(x) {apply(x, 2, var)}))

eqm_n_mcmc <- matrix(c(eqm_b0_50_mcmc, eqm_b0_100_mcmc, eqm_b0_200_mcmc,
                       eqm_s_50_mcmc, eqm_s_100_mcmc, eqm_s_200_mcmc,
                       eqm_p_50_mcmc, eqm_p_100_mcmc, eqm_p_200_mcmc,
                       eqm_t_50_mcmc, eqm_t_100_mcmc, eqm_t_200_mcmc), nrow = 3)

# Viés do estimador

vies_n_mcmc <- eqm_n_mcmc - var_n_mcmc

rownames(media_n_mcmc) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(var_n_mcmc) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(eqm_n_mcmc) <- c("$n = 50$", "$n = 100$", "$n = 200$")
rownames(vies_n_mcmc) <- c("$n = 50$", "$n = 100$", "$n = 200$")

colnames(eqm_n_mcmc) <- colnames(media_n_mcmc) <- colnames(var_n_mcmc) <- colnames(vies_n_mcmc) <- c("$\\hat{\\beta_{0}}$","$\\hat{\\sigma^{2}}$",
                                                                                                     "$\\hat{\\phi}$", "$\\hat{\\tau^{2}}$")
mat_media_n_mcmc <- xtable(media_n_mcmc, auto = T, caption = "Média dos estimadores por MCMC",
                           digits = 4, label = "tab:mediamcmc")
mat_var_n_mcmc <- xtable(var_n_mcmc, auto = T, caption = "Variância dos estimadores por MCMC",
                         digits = 4, label = "tab:varestnmcmc")
mat_eqm_n_mcmc <- xtable(eqm_n_mcmc, auto = T, caption = "Erro Quadrático Médio dos estimadores por MCMC",
                         digits = 4, label = "tab:eqmnmcmc")
mat_vies_n_mcmc <- xtable(vies_n_mcmc, auto = T, caption = "Viés dos estimadores por MCMC",
                          digits = 4, label = "tab:viesparnmcmc")

# Matrizes com os estimadores:

print(
    xtable(
        data.frame(
            Var = c("$\\hat{\\beta_{0}}$", NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA,
                    "$\\hat{\\phi}$", NA, NA, "$\\hat{\\tau^{2}}$", NA, NA),
            n = c(rep(c("$n = 50$", "$n = 100$", "$n = 200$"), 4)),
            `Média` = c(media_n[, 1], media_n[, 2], media_n[, 3], media_n[, 4]),
            `Variância` = c(var_n[, 1], var_n[, 2], var_n[, 3], var_n[, 4]),
            `EQM` = c(eqm_n[, 1], eqm_n[, 2], eqm_n[, 3], eqm_n[, 4]),
            `Viés` = c(vies_n[, 1], vies_n[, 2], vies_n[, 3], vies_n[, 4])
        ),
        digits = 4,
        caption = "Estimadores para diferentes tamanhos de amostras",
        label = "tab:estamo"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F
)

print(
    xtable(
        data.frame(
            `Parâmetro` = c("$\\hat{\\beta_{0}}$", NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA,
                    "$\\hat{\\phi}$", NA, NA, "$\\hat{\\tau^{2}}$", NA, NA),
            n = c(rep(c("$n = 50$", "$n = 100$", "$n = 200$"), 4)),
            `Média` = c(media_n_mcmc[, 1], media_n_mcmc[, 2], media_n_mcmc[, 3], media_n_mcmc[, 4]),
            `Variância` = c(var_n_mcmc[, 1], var_n_mcmc[, 2], var_n_mcmc[, 3], var_n_mcmc[, 4]),
            `EQM` = c(eqm_n_mcmc[, 1], eqm_n_mcmc[, 2], eqm_n_mcmc[, 3], eqm_n_mcmc[, 4]),
            `Viés` = c(vies_n_mcmc[, 1], vies_n_mcmc[, 2], vies_n_mcmc[, 3], vies_n_mcmc[, 4])
        ),
        digits = 4,
        caption = "Estimadores MCMC para diferentes tamanhos de amostras",
        label = "tab:estamomcmc"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F
)

# Ajustes para modelos de covariância incorretos:

resultados_esfe <- parLapply(cl, 1:1000, simulation_function, cov.model = "spherical")
resultados_mat_1 <- parLapply(cl, 1:1000, simulation_function, cov.model = "matern", kappa = 1)
resultados_mat_2 <- parLapply(cl, 1:1000, simulation_function, cov.model = "matern", kappa = 2)

resultados_esfe_mcmc <- parLapply(cl, 1:1000, simulation_function_2, cov.model = "spherical")
resultados_mat_1_mcmc <- parLapply(cl, 1:1000, simulation_function_2, cov.model = "matern", kappa = 1)
resultados_mat_2_mcmc <- parLapply(cl, 1:1000, simulation_function_2, cov.model = "matern", kappa = 2)

resultados_esfe <- do.call(rbind, resultados_esfe)
resultados_mat_1 <- do.call(rbind, resultados_mat_1)
resultados_mat_2 <- do.call(rbind, resultados_mat_2)

resultados_esfe_mcmc <- do.call(rbind, resultados_esfe_mcmc)
resultados_mat_1_mcmc <- do.call(rbind, resultados_mat_1_mcmc)
resultados_mat_2_mcmc <- do.call(rbind, resultados_mat_2_mcmc)


media_cov <- do.call(rbind, (lapply(list(resultados_100,
                                         resultados_mat_1,
                                         resultados_mat_2,
                                         resultados_esfe),
                                    colMeans)))

media_cov_mcmc <- do.call(rbind, (lapply(list(resultados_100_mcmc,
                                              resultados_mat_1_mcmc,
                                              resultados_mat_2_mcmc,
                                              resultados_esfe_mcmc),
                                         colMeans)))

var_cov <- do.call(rbind, lapply(list(resultados_100,
                                      resultados_mat_1,
                                      resultados_mat_2,
                                      resultados_esfe), FUN = function(x) {apply(x, 2, var)}))

var_cov_mcmc <- do.call(rbind, lapply(list(resultados_100_mcmc,
                                           resultados_mat_1_mcmc,
                                           resultados_mat_2_mcmc,
                                           resultados_esfe_mcmc), FUN = function(x) {apply(x, 2, var)}))

eqm_b0_esfe <- mean((resultados_esfe[, 1] - 2)^2)
eqm_b0_mat1 <- mean((resultados_mat_1[, 1] - 2)^2)
eqm_b0_mat2 <- mean((resultados_mat_2[, 1] - 2)^2)

eqm_s_esfe <- mean((resultados_esfe[, 2] - 0.5)^2)
eqm_s_mat1 <- mean((resultados_mat_1[, 2] - 0.5)^2)
eqm_s_mat2 <- mean((resultados_mat_2[, 2] - 0.5)^2)

eqm_p_esfe <- mean((resultados_esfe[, 3] - 30)^2)
eqm_p_mat1 <- mean((resultados_mat_1[, 3] - 30)^2)
eqm_p_mat2 <- mean((resultados_mat_2[, 3] - 30)^2)

eqm_t_esfe <- mean((resultados_esfe[, 4] - 0.05)^2)
eqm_t_mat1 <- mean((resultados_mat_1[, 4] - 0.05)^2)
eqm_t_mat2 <- mean((resultados_mat_2[, 4] - 0.05)^2)

eqm_cov <- matrix(c(eqm_b0_100, eqm_b0_mat1, eqm_b0_mat2, eqm_b0_esfe,
                    eqm_s_100, eqm_s_mat1, eqm_s_mat2, eqm_s_esfe,
                    eqm_p_100, eqm_p_mat1, eqm_p_mat2, eqm_p_esfe,
                    eqm_t_100, eqm_t_mat1, eqm_t_mat2, eqm_t_esfe), nrow = 4)

eqm_b0_esfe_mcmc <- mean((resultados_esfe_mcmc[, 1] - 2)^2)
eqm_b0_mat1_mcmc <- mean((resultados_mat_1_mcmc[, 1] - 2)^2)
eqm_b0_mat2_mcmc <- mean((resultados_mat_2_mcmc[, 1] - 2)^2)

eqm_s_esfe_mcmc <- mean((resultados_esfe_mcmc[, 2] - 0.5)^2)
eqm_s_mat1_mcmc <- mean((resultados_mat_1_mcmc[, 2] - 0.5)^2)
eqm_s_mat2_mcmc <- mean((resultados_mat_2_mcmc[, 2] - 0.5)^2)

eqm_p_esfe_mcmc <- mean((resultados_esfe_mcmc[, 3] - 30)^2)
eqm_p_mat1_mcmc <- mean((resultados_mat_1_mcmc[, 3] - 30)^2)
eqm_p_mat2_mcmc <- mean((resultados_mat_2_mcmc[, 3] - 30)^2)

eqm_t_esfe_mcmc <- mean((resultados_esfe_mcmc[, 4] - 0.05)^2)
eqm_t_mat1_mcmc <- mean((resultados_mat_1_mcmc[, 4] - 0.05)^2)
eqm_t_mat2_mcmc <- mean((resultados_mat_2_mcmc[, 4] - 0.05)^2)

eqm_cov_mcmc <- matrix(c(eqm_b0_100_mcmc, eqm_b0_mat1_mcmc, eqm_b0_mat2_mcmc, eqm_b0_esfe_mcmc,
                         eqm_s_100_mcmc, eqm_s_mat1_mcmc, eqm_s_mat2_mcmc, eqm_s_esfe_mcmc,
                         eqm_p_100_mcmc, eqm_p_mat1_mcmc, eqm_p_mat2_mcmc, eqm_p_esfe_mcmc,
                         eqm_t_100_mcmc, eqm_t_mat1_mcmc, eqm_t_mat2_mcmc, eqm_t_esfe_mcmc), nrow = 4)


# Viés do estimador

vies_cov <- eqm_cov - var_cov

vies_cov_mcmc <- eqm_cov_mcmc - var_cov_mcmc

rownames(media_cov) <- c("Exponencial", "Esférica", "$\\text{Matèrn,}\\kappa = 1$", "$\\text{Matèrn,}\\kappa = 2$")
rownames(var_cov) <- c("Exponencial", "Esférica", "$\\text{Matèrn,}\\kappa = 1$", "$\\text{Matèrn,}\\kappa = 2$")
rownames(eqm_cov) <- c("Exponencial", "Esférica", "$\\text{Matèrn,}\\kappa = 1$", "$\\text{Matèrn,}\\kappa = 2$")
rownames(vies_cov) <- c("Exponencial", "Esférica", "$\\text{Matèrn,}\\kappa = 1$", "$\\text{Matèrn,}\\kappa = 2$")

colnames(eqm_cov) <- colnames(media_cov)
colnames(vies_cov) <- colnames(media_cov)

mat_media_cov <- xtable(media_cov, auto = T, caption = "Média dos estimadores",
                        digits = 4, label = "tab:mediacov")
mat_var_cov <- xtable(var_cov, auto = T, caption = "Variância dos estimadores",
                      digits = 4, label = "tab:varestcov")
mat_eqm_cov <- xtable(eqm_cov, auto = T, caption = "Erro Quadrático Médio dos estimadores",
                      digits = 4, label = "tab:eqmcov")
mat_vies_cov <- xtable(vies_cov, auto = T, caption = "Viés dos estimadores",
                       digits = 4, label = "tab:viesparcov")

# Matrizes com os estimadores:

print(
    xtable(
        data.frame(
            Var = c("$\\hat{\\beta_{0}}$", NA, NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA, NA,
                    "$\\hat{\\phi}$", NA, NA, NA, "$\\hat{\\tau^{2}}$", NA, NA, NA),
            n = c(rep(c("Exponencial", "Esférica", "$\\text{Matèrn, }\\kappa = 1$", "$\\text{Matèrn, }\\kappa = 2$"), 4)),
            `Média` = c(media_cov[, 1], media_cov[, 2], media_cov[, 3], media_cov[, 4]),
            `Variância` = c(var_cov[, 1], var_cov[, 2], var_cov[, 3], var_cov[, 4]),
            `EQM` = c(eqm_cov[, 1], eqm_cov[, 2], eqm_cov[, 3], eqm_cov[, 4]),
            `Viés` = c(vies_cov[, 1], vies_cov[, 2], vies_cov[, 3], vies_cov[, 4])
        ),
        digits = 4,
        caption = "Estimadores para diferentes funções de covariância",
        label = "tab:estcov"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F
)

print(
    xtable(
        data.frame(
            Var = c("$\\hat{\\beta_{0}}$", NA, NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA, NA,
                    "$\\hat{\\phi}$", NA, NA, NA, "$\\hat{\\tau^{2}}$", NA, NA, NA),
            n = c(rep(c("Exponencial", "Esférica", "$\\text{Matèrn, }\\kappa = 1$", "$\\text{Matèrn, }\\kappa = 2$"), 4)),
            `Média` = c(media_cov_mcmc[, 1], media_cov_mcmc[, 2], media_cov_mcmc[, 3], media_cov_mcmc[, 4]),
            `Variância` = c(var_cov_mcmc[, 1], var_cov_mcmc[, 2], var_cov_mcmc[, 3], var_cov_mcmc[, 4]),
            `EQM` = c(eqm_cov_mcmc[, 1], eqm_cov_mcmc[, 2], eqm_cov_mcmc[, 3], eqm_cov_mcmc[, 4]),
            `Viés` = c(vies_cov_mcmc[, 1], vies_cov_mcmc[, 2], vies_cov_mcmc[, 3], vies_cov_mcmc[, 4])
        ),
        digits = 4,
        caption = "Estimadores MCMC para diferentes funções de covariância",
        label = "tab:estcovmcmc"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F
)

# Ajustes para espaço amostral maior/menor

resultados_menor <- parLapply(cl, 1:1000, simulation_function, xlims = c(0, 50), ylims = c(0, 50))
resultados_maior <- parLapply(cl, 1:1000, simulation_function, xlims = c(0, 500), ylims = c(0, 500))

resultados_menor_mcmc <- parLapply(cl, 1:1000, simulation_function_2, xlims = c(0, 50), ylims = c(0, 50))
resultados_maior_mcmc <- parLapply(cl, 1:1000, simulation_function_2, xlims = c(0, 500), ylims = c(0, 500))

resultados_menor <- do.call(rbind, resultados_menor)
resultados_maior <- do.call(rbind, resultados_maior)

resultados_menor_mcmc <- do.call(rbind, resultados_menor_mcmc)
resultados_maior_mcmc <- do.call(rbind, resultados_maior_mcmc)

media_grid <- do.call(rbind, (lapply(list(resultados_menor,
                                          resultados_100,
                                          resultados_maior),
                                     colMeans)))

media_grid_mcmc <- do.call(rbind, (lapply(list(resultados_menor_mcmc,
                                               resultados_100_mcmc,
                                               resultados_maior_mcmc),
                                     colMeans)))

var_grid <- do.call(rbind, lapply(list(resultados_menor,
                                       resultados_100,
                                       resultados_maior), FUN = function(x) {apply(x, 2, var)}))

var_grid_mcmc <- do.call(rbind, lapply(list(resultados_menor_mcmc,
                                            resultados_100_mcmc,
                                            resultados_maior_mcmc), FUN = function(x) {apply(x, 2, var)}))

eqm_b0_menor <- mean((resultados_menor[, 1] - 2)^2)
eqm_b0_maior <- mean((resultados_maior[, 1] - 2)^2)

eqm_s_menor <- mean((resultados_menor[, 2] - 0.5)^2)
eqm_s_maior <- mean((resultados_maior[, 2] - 0.5)^2)

eqm_p_menor <- mean((resultados_menor[, 3] - 30)^2)
eqm_p_maior <- mean((resultados_maior[, 3] - 30)^2)

eqm_t_menor <- mean((resultados_menor[, 4] - 0.05)^2)
eqm_t_maior <- mean((resultados_maior[, 4] - 0.05)^2)

eqm_b0_menor_mcmc <- mean((resultados_menor_mcmc[, 1] - 2)^2)
eqm_b0_maior_mcmc <- mean((resultados_maior_mcmc[, 1] - 2)^2)

eqm_s_menor_mcmc <- mean((resultados_menor_mcmc[, 2] - 0.5)^2)
eqm_s_maior_mcmc <- mean((resultados_maior_mcmc[, 2] - 0.5)^2)

eqm_p_menor_mcmc <- mean((resultados_menor_mcmc[, 3] - 30)^2)
eqm_p_maior_mcmc <- mean((resultados_maior_mcmc[, 3] - 30)^2)

eqm_t_menor_mcmc <- mean((resultados_menor_mcmc[, 4] - 0.05)^2)
eqm_t_maior_mcmc <- mean((resultados_maior_mcmc[, 4] - 0.05)^2)

eqm_grid <- matrix(c(eqm_b0_menor, eqm_b0_100, eqm_b0_maior,
                     eqm_s_menor, eqm_s_100, eqm_s_maior,
                     eqm_p_menor, eqm_p_100, eqm_p_maior,
                     eqm_t_menor, eqm_t_100, eqm_t_maior), nrow = 3)

eqm_grid_mcmc <- matrix(c(eqm_b0_menor_mcmc, eqm_b0_100_mcmc, eqm_b0_maior_mcmc,
                          eqm_s_menor_mcmc, eqm_s_100_mcmc, eqm_s_maior_mcmc,
                          eqm_p_menor_mcmc, eqm_p_100_mcmc, eqm_p_maior_mcmc,
                          eqm_t_menor_mcmc, eqm_t_100_mcmc, eqm_t_maior_mcmc), nrow = 3)

# Viés do estimador

vies_grid <- eqm_grid - var_grid

vies_grid_mcmc <- eqm_grid_mcmc - var_grid_mcmc

rownames(media_grid) <- c("$50 \\times 50$", "$100 \\times 100$", "$200 \\times 200$")
rownames(var_grid) <- c("$50 \\times 50$", "$100 \\times 100$", "$200 \\times 200$")
rownames(eqm_grid) <- c("$50 \\times 50$", "$100 \\times 100$", "$200 \\times 200$")
rownames(vies_grid) <- c("$50 \\times 50$", "$100 \\times 100$", "$200 \\times 200$")

colnames(eqm_grid) <- colnames(media_grid)
colnames(vies_grid) <- colnames(media_grid)

mat_media_grid <- xtable(media_grid, auto = T, caption = "Média dos estimadores",
                         digits = 4, label = "tab:mediagrid")
mat_var_grid <- xtable(var_grid, auto = T, caption = "Variância dos estimadores",
                       digits = 4, label = "tab:varestgrid")
mat_eqm_grid <- xtable(eqm_grid, auto = T, caption = "Erro Quadrático Médio dos estimadores",
                       digits = 4, label = "tab:eqmgrid")
mat_vies_grid <- xtable(vies_grid, auto = T, caption = "Viés dos estimadores",
                        digits = 4, label = "tab:viespargrid")

print(
    xtable(
        data.frame(
            `Parâmetro` = c("$\\hat{\\beta_{0}}$", NA, NA, "$\\hat{\\sigma^{2}}$", NA, NA,
                    "$\\hat{\\phi}$", NA, NA, "$\\hat{\\tau^{2}}$", NA, NA),
            `Correlação` = c(rep(c("$50 \\times 50$", "$100 \times 100$", "$500 \times 500$"), 4)),
            `Média` = c(media_grid[, 1], media_grid[, 2], media_grid[, 3], media_grid[, 4]),
            `Variância` = c(var_grid[, 1], var_grid[, 2], var_grid[, 3], var_grid[, 4]),
            `EQM` = c(eqm_grid[, 1], eqm_grid[, 2], eqm_grid[, 3], eqm_grid[, 4]),
            `Viés` = c(vies_grid[, 1], vies_grid[, 2], vies_grid[, 3], vies_grid[, 4])
        ),
        digits = 4,
        align = "l|l|lrrrr|",
        caption = "Estimadores para diferentes regiões amostrais",
        label = "tab:estgrid"
    ),
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    include.rownames = F,
    hline.after = c(-1, 0, 3, 6, 9, 12)
)

# Stop the cluster
stopCluster(cl)
