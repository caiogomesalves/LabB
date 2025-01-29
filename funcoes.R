library(geoR)
library(geoRglm)
library(MASS)
library(lme4)
library(bbmle)

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

newton.raphson <- function(initial, escore, hessiano, tol=0.0001,
                           max.iter, n.dim, ...) {
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
                     inv.Sigma = inv.Sigma)
  )
  return(as.numeric(dens))
}

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
        grad <- t(part1) - t(b) %*% inv.Sigma
    }
    )
}

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
    "negbin" = {
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
    }
    )
    return(Hess)
}

laplace <- function(Q.b, gr, hess, otimizador,
                    n.dim, method.integrate, ...){
    log.integral <- -sqrt(.Machine$double.xmax)
    inicial <- rep(0, n.dim)
    pred <- NULL
    if (method.integrate == "BFGS") {
        temp <- try(optim(inicial, Q.b, gr = gr, ...,
                          method = otimizador, hessian = TRUE,
                          control = list(fnscale = -1)), silent = TRUE)
    } else if (method.integrate == "NR") {
        temp <- try(newton.raphson(initial = inicial, escore = gr,
                                   hessiano = hess, n.dim = n.dim,
                                   max.iter = 100, ...), silent = TRUE)
    } else if (method.integrate == "QNR") {
        temp <- try(qq.newton.raphson(initial = inicial, escore = gr,
                                      hessiano = hess, n.dim = n.dim,
                                      max.iter = 100, ...), silent = TRUE)
    }
    if (class(temp) != "try-error" && method.integrate == "BFGS") {
        log.integral <- temp$value + ((n.dim / 2) *
                                      log(2 * pi) - 0.5 *
                                      determinant(-temp$hessian)$modulus)
        pred <- temp$par
    } else if (class(temp) != "try-error" && method.integrate == "NR") {
        value <- Q.b(b = temp$solution, ...)
        log.integral <- value + ((n.dim / 2) * log(2 * pi) - 0.5 *
                                 determinant(-temp[1][[1]])$modulus)
        pred <- temp$solution
    }
    return(list(log.integral = log.integral, pred = pred))
}

loglik.sglmm <- function(par, Y, X, kappa, nugget,
                         mat.dist, cov.model,
                         family, method.integrate = "NR",
                         ntrial = 1, offset = NA){
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
    return(-I[[1]])
}

preditos <- function(par, Y, X, kappa, nugget, mat.dist, cov.model,
                     family, method.integrate = "NR",
                     ntrial = 1, offset = NA){
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
    return(I$pred)
}

start.values.sglmm <- function(formula, data, coords, nugget,
                               family, ntrial = 1, offset = 1){
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
    return(saida)
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
                         family = family,
                         method.integrate = method.integrate,
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
        saida[8][[1]] <- ifelse(cov.model == "matern", exp(kappa), "NULL")
        saida[9][[1]] <- estimativas
        return(saida)}
    if(predict == FALSE){
        return(estimativas)
    }
}
