###### (1) Main function ----
fit701M.nb <- function(
  xv,
  isoyv,
  isowv,
  etx,
  dtx,
  int,
  constraints,
  eps = 1e-8
) {
  mtx <- dtx / etx # matrix of death rates
  mtx[mtx == 0] <- 10^(-8)

  n <- length(xv) # number of ages
  mY <- length(unique(isoyv))
  mW <- length(unique(isowv)) # number of years
  m <- nrow(etx)
  R <- dim(mtx)[3]

  # initialise parameter vectors
  beta1v <- matrix(rep(int, each = R), nrow = n, ncol = R)
  beta2v = beta3v = beta4v = beta5v = beta6v <- (1:n) * 0
  kappa2v <- matrix(0, nrow = mY, ncol = R)
  kappa3v <- matrix(0, nrow = mW, ncol = R)

  # Stage 0
  mx <- mean(xv)
  for (j in 1:n) {
    if (all(int == 0)) {
      beta1v[j, ] <- colMeans(log(mtx[, j, ]), na.rm = TRUE)
    }
    beta2v[j] <- 1
    beta3v[j] <- 1
  }

  for (r in 1:R) {
    kappa2v[, r] <- sapply(1:mY, function(t) {
      (sum(log(mtx[isoyv == t, , r]), na.rm = TRUE) - mW * sum(beta1v[, r])) /
        mW
    })
    kappa3v[, r] <- sapply(1:mW, function(w) {
      (sum(log(mtx[isowv == w, , r]), na.rm = TRUE) - mY * sum(beta1v[, r])) /
        mY
    })
  }

  # Dispersion values
  ds.x <- rep(log(5), n) - mean(rep(log(5), n))
  ds.r <- rep(log(100), R)
  exp.ds <- exp(outer(ds.x, ds.r, "+"))
  exp.ds <- replicate(m, exp.ds, simplify = "array") # gives an n x R x T array
  exp.ds <- aperm(exp.ds, c(3, 1, 2))

  # Stage 1: iterate
  l0 <- -100
  l1 <- -99
  iter <- 0

  while (abs((l1 - l0) / l0) > eps) {
    iter <- iter + 1
    l0 <- l1

    ### (1) Stage 1B2 optimise over the beta2(x) ----
    for (j in 1:(n)) {
      dv <- dtx[, j, ]
      ev <- etx[, j, ]
      phi <- exp.ds[, j, ]
      beta2v[j] <- llmaxM2B1.nb(
        b1 = beta1v[j, ],
        b2 = beta2v[j - 1 * (iter == 1 & j != 1)],
        b3 = beta3v[j],
        k2 = kappa2v,
        k3 = kappa3v,
        dv = dv,
        ev = ev,
        phi = phi
      )
    }

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (2) Stage 1K2 optimise over the kappa2(t,r) ----
    for (i in 1:mY) {
      dv <- dtx[isoyv == i, , ]
      ev <- etx[isoyv == i, , ]
      phi <- exp.ds[isoyv == i, , ]
      kappa2v[i, ] <- llmaxM2D1.nb(
        b1 = beta1v,
        b2 = beta2v,
        b3 = beta3v,
        k2 = kappa2v[i, ],
        k3 = kappa3v,
        dv = dv,
        ev = ev,
        phi = phi
      )
    }

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (3) Stage 1B3 optimize over the beta3(x) ----
    for (j in 1:n) {
      dv <- dtx[, j, ]
      ev <- etx[, j, ]
      phi <- exp.ds[, j, ]
      beta3v[j] <- llmaxM2B2.nb(
        b1 = beta1v[j, ],
        b2 = beta2v[j],
        b3 = beta3v[j],
        k2 = kappa2v,
        k3 = kappa3v,
        dv = dv,
        ev = ev,
        phi = phi
      )
    }

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (4) Stage 1K3 optimise over the kappa3(t,r) ----
    for (i in 1:mW) {
      dv <- dtx[isowv == i, , ] # actual deaths
      ev <- etx[isowv == i, , ] # exposure
      phi <- exp.ds[isowv == i, , ]
      kappa3v[i, ] <- llmaxM2D2.nb(
        b1 = beta1v,
        b2 = beta2v,
        b3 = beta3v,
        k2 = kappa2v,
        k3 = kappa3v[i, ],
        dv = dv,
        ev = ev,
        phi = phi
      )
    }

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (5) Stage 1A  optimise over the beta1(x) ----
    for (j in n:1) {
      dv <- dtx[, j, ]
      ev <- etx[, j, ]
      phi <- exp.ds[, j, ]
      beta1v[j, ] <- llmaxM2B0.nb(
        b1 = beta1v[j + ifelse(iter == 1 & j != n, 1, 0), ],
        b2 = beta2v[j],
        b3 = beta3v[j],
        k2 = kappa2v,
        k3 = kappa3v,
        dv = dv,
        ev = ev,
        phi = phi
      )
    }

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (6) Stage 1F  optimize over the phi(r) ----
    ds.r <- llmaxM2Phir.nb(
      mu = mhat,
      dv = dtx,
      ev = etx,
      phix = ds.x,
      phir = ds.r
    )

    exp.ds <- exp(outer(ds.x, ds.r, "+"))
    exp.ds <- replicate(m, exp.ds, simplify = "array") # gives an n x R x T array
    exp.ds <- aperm(exp.ds, c(3, 1, 2))

    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (7) Stage 1F  optimize over the phi(x) ----
    ds.x <- llmaxM2Phix.nb(
      mu = mhat,
      dv = dtx,
      ev = etx,
      phix = ds.x,
      phir = ds.r
    )

    exp.ds <- exp(outer(ds.x, ds.r, "+"))
    exp.ds <- replicate(m, exp.ds, simplify = "array") # gives an n x R x T array
    exp.ds <- aperm(exp.ds, c(3, 1, 2))

    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    ### (8) Constraints ----
    if (constraints == "ALL") {
      fac21 <- kappa2v[1, ]
      fac31 <- kappa3v[1, ]
      fac22 <- beta2v[1] # sum(beta2v)
      fac32 <- beta3v[1] # sum(beta3v)
      fac71 <- mean(ds.x)
      kappa2v <- fac22 * sweep(kappa2v, 2, fac21, '-')
      kappa3v <- fac32 * sweep(kappa3v, 2, fac31, '-')
      beta2v <- beta2v / fac22
      beta3v <- beta3v / fac32
      beta1v <- beta1v +
        (beta2v * fac22) %o% fac21 +
        (beta3v * fac32) %o% fac31
      ds.r <- ds.r + fac71
      ds.x <- ds.x - fac71
    }

    exp.ds <- exp(outer(ds.x, ds.r, "+"))
    exp.ds <- replicate(m, exp.ds, simplify = "array") # gives an n x R x T array
    exp.ds <- aperm(exp.ds, c(3, 1, 2))

    mhat <- mtx * NA
    for (j in 1:n) {
      mhat[, j, ] <- exp(sweep(
        beta2v[j] * kappa2v[isoyv, ] + beta3v[j] * kappa3v[isowv, ],
        2,
        beta1v[j, ],
        '+'
      ))
    }
    epsilon <- (dtx - etx * mhat) / sqrt(etx * mhat + (etx * mhat)^2 / exp.ds)
    l1 <- sum(
      lgamma(exp.ds + dtx) -
        lgamma(exp.ds) -
        lgamma(dtx + 1) +
        dtx * log(etx * mhat) -
        dtx * log(exp.ds + etx * mhat) +
        exp.ds * log(exp.ds) -
        exp.ds * log(exp.ds + etx * mhat),
      na.rm = TRUE
    )

    cat(l1, ' -> ')
  }

  ### Calculate hessians baseline model
  lcb <- list(
    beta1 = beta1v,
    beta2 = beta2v,
    beta3 = beta3v,
    kappa2 = kappa2v,
    kappa3 = kappa3v,
    phix = ds.x,
    phir = ds.r,
    mhat = mhat
  )
  H.beta1v <- hessian.beta1v(dv = dtx, ev = etx, lca = lcb)
  H.beta2v <- hessian.beta2v(dv = dtx, ev = etx, lca = lcb)
  H.kappa2v <- hessian.kappa2v(dv = dtx, ev = etx, lca = lcb)
  H.beta3v <- hessian.beta3v(dv = dtx, ev = etx, lca = lcb)
  H.kappa3v <- hessian.kappa3v(dv = dtx, ev = etx, lca = lcb)

  H.base <- cbind(H.beta1v, H.beta2v, H.kappa2v, H.beta3v, H.kappa3v)

  # set dimnames
  colnames(H.base) = rownames(H.base) <-
    c(
      outer(dimnames(dtx)[[2]], dimnames(dtx)[[3]], FUN = function(x, y) {
        paste0('Alpha_Age', x, "_Region", y)
      }),
      paste0('Beta2_Age', dimnames(dtx)[[2]]),
      outer(unique(isoyv), dimnames(dtx)[[3]], FUN = function(x, y) {
        paste0('Kappa2_Year', x, "_Region", y)
      }),
      paste0('Beta3_Age', dimnames(dtx)[[2]]),
      outer(unique(isowv), dimnames(dtx)[[3]], FUN = function(x, y) {
        paste0('Kappa3_Week', x, "_Region", y)
      })
    )

  # Remove some columns due to ID constraints
  ind1.rm <- c(
    length(beta1v) + 1,
    length(beta1v) + length(beta2v) + (0:(R - 1)) * max(isoyv) + 1,
    length(beta1v) + length(beta2v) + length(kappa2v) + 1,
    length(beta1v) +
      length(beta2v) +
      length(kappa2v) +
      length(beta3v) +
      (0:(R - 1)) * max(isowv) +
      1
  )
  H.base <- H.base[-ind1.rm, -ind1.rm]

  # Calculate number of parameters and deduct the number of constraints
  if (constraints == "ALL") {
    npar <- ncol(H.base)
  }

  # Calculate the BIC
  AIC <- -2 * l1 + 2 * npar
  BIC <- -2 * l1 + log(length(dtx)) * npar

  list(
    beta1 = beta1v,
    beta2 = beta2v,
    beta3 = beta3v,
    kappa2 = kappa2v,
    kappa3 = kappa3v,
    phix = ds.x,
    phir = ds.r,
    epsilon = epsilon,
    mhat = mhat,
    ll = l1,
    H.base = H.base,
    npar = npar,
    ind1.rm = ind1.rm,
    AIC = AIC,
    BIC = BIC
  )
}

###### (2) Newton Raphson iterations baseline model -----
llmaxM2B0.nb <- function(b1, b2, b3, k2, k3, dv, ev, phi) {
  b11 <- b1
  b10 <- b11 - 1
  thetat <- ev *
    exp(b2 * k2[isoyv, ] + b3 * k3[isowv, ])
  s1 <- colSums(dv, na.rm = TRUE)

  while (max(abs(b11 - b10)) > 0.0001) {
    b10 <- b11
    thetat.exp <- sweep(thetat, 2, exp(b10), '*')
    f0 <- colSums((dv + phi) * (1 - phi / (phi + thetat.exp)), na.rm = TRUE) -
      s1
    df0 <- colSums(
      (dv + phi) * (phi * thetat.exp) / ((phi + thetat.exp)^2),
      na.rm = TRUE
    )
    b11 <- b10 - f0 / df0
    b11
  }
  b11
}
llmaxM2B1.nb <- function(b1, b2, b3, k2, k3, dv, ev, phi) {
  b21 <- b2
  b20 <- b21 - 1
  thetat <- ev *
    exp(sweep(b3 * k3[isowv, ], 2, b1, '+'))
  s1 <- sum(dv * k2[isoyv, ], na.rm = TRUE)

  while (abs(b21 - b20) > 0.0001) {
    b20 <- b21
    thetat.exp <- thetat * exp(b20 * k2[isoyv, ])
    f0 <- sum(
      (dv + phi) * (1 - phi / (phi + thetat.exp)) * k2[isoyv, ],
      na.rm = TRUE
    ) -
      s1
    df0 <- sum(
      (dv + phi) * (phi * thetat.exp) / ((phi + thetat.exp)^2) * k2[isoyv, ]^2,
      na.rm = TRUE
    )
    b21 <- b20 - f0 / df0
    b21
  }
  b21
}
llmaxM2D1.nb <- function(b1, b2, b3, k2, k3, dv, ev, phi) {
  k21 <- k2
  k20 <- k21 - 1
  k3 <- k3
  thetat <- ev * exp(sweep(aperm(b3 %o% k3, c(2, 1, 3)), c(2, 3), b1, '+'))
  s1 <- apply(sweep(dv, 2, b2, '*'), 3, sum, na.rm = TRUE)

  while (max(abs(k21 - k20)) > 0.0001) {
    k20 <- k21
    thetat.exp <- sweep(thetat, c(2, 3), exp(b2 %o% k20), '*')
    f0 <- apply(
      sweep((dv + phi) * (1 - phi / (phi + thetat.exp)), 2, b2, '*'),
      3,
      sum,
      na.rm = TRUE
    ) -
      s1
    df0 <- apply(
      sweep(
        (dv + phi) * (phi * thetat.exp / ((phi + thetat.exp)^2)),
        2,
        b2^2,
        '*'
      ),
      3,
      sum,
      na.rm = TRUE
    )
    k21 <- k20 - f0 / df0
    k21
  }
  k21
}
llmaxM2B2.nb <- function(b1, b2, b3, k2, k3, dv, ev, phi) {
  b31 <- b3
  b30 <- b31 - 1
  thetat <- ev *
    exp(sweep(b2 * k2[isoyv, ], 2, b1, '+'))
  s1 <- sum(dv * k3[isowv, ], na.rm = TRUE)

  while (abs(b31 - b30) > 0.0001) {
    b30 <- b31
    thetat.exp <- thetat * exp(b30 * k3[isowv, ])
    f0 <- sum(
      (dv + phi) * (1 - phi / (phi + thetat.exp)) * k3[isowv, ],
      na.rm = TRUE
    ) -
      s1
    df0 <- sum(
      (dv + phi) * (phi * thetat.exp) / ((phi + thetat.exp)^2) * k3[isowv, ]^2,
      na.rm = TRUE
    )
    b31 <- b30 - f0 / df0
  }
  b31
}
llmaxM2D2.nb <- function(b1, b2, b3, k2, k3, dv, ev, phi) {
  k31 <- k3
  k30 <- k31 - 1
  thetat <- ev * exp(sweep(aperm(b2 %o% k2, c(2, 1, 3)), c(2, 3), b1, '+'))
  s1 <- apply(sweep(dv, 2, b3, '*'), 3, sum, na.rm = TRUE)

  while (max(abs(k31 - k30)) > 0.0001) {
    k30 <- k31
    thetat.exp <- sweep(thetat, c(2, 3), exp(b3 %o% k30), '*')
    f0 <- apply(
      sweep((dv + phi) * (1 - phi / (phi + thetat.exp)), 2, b3, '*'),
      3,
      sum,
      na.rm = TRUE
    ) -
      s1
    df0 <- apply(
      sweep(
        (dv + phi) * (phi * thetat.exp / ((phi + thetat.exp)^2)),
        2,
        b3^2,
        '*'
      ),
      3,
      sum,
      na.rm = TRUE
    )
    k31 <- k30 - f0 / df0
  }
  k31
}
llmaxM2Phir.nb <- function(mu, dv, ev, phix, phir) {
  phir1 <- phir
  phir0 <- phir - 0.1
  exp.phi0 <- exp(outer(phix, phir0, "+"))
  exp.phi0 <- replicate(nrow(dv), exp.phi0, simplify = "array") # gives an n x R x T array
  exp.phi0 <- aperm(exp.phi0, c(3, 1, 2))

  while (max(abs(phir1 - phir0)) > 0.001) {
    phir0 <- phir1
    exp.phi0 <- exp(outer(phix, phir0, "+"))
    exp.phi0 <- replicate(nrow(dv), exp.phi0, simplify = "array") # gives an n x R x T array
    exp.phi0 <- aperm(exp.phi0, c(3, 1, 2))

    dg1 <- digamma(dv + exp.phi0)
    tg1 <- trigamma(dv + exp.phi0)
    dg2 <- digamma(exp.phi0)
    tg2 <- trigamma(exp.phi0)
    c1 <- dv + exp.phi0
    c2 <- ev * mu + exp.phi0
    lc2 <- log(c2)

    f0 <- apply(
      (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
      3,
      sum,
      na.rm = TRUE
    )
    df0 <- apply(
      ((tg1 - tg2 + 1 / exp.phi0 - 2 / c2 + c1 * 1 / (c2^2)) * exp.phi0^2) +
        (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
      3,
      sum,
      na.rm = TRUE
    )

    phir1 <- phir0 - f0 / df0
  }
  phir1
}
llmaxM2Phix.nb <- function(mu, dv, ev, phix, phir) {
  phix[length(phix)] <- phix[length(phix) - 1]
  phix1 <- phix
  phix0 <- phix - 0.1
  exp.phi0 <- exp(outer(phix0, phir, "+"))
  exp.phi0 <- replicate(nrow(dv), exp.phi0, simplify = "array") # gives an n x R x T array
  exp.phi0 <- aperm(exp.phi0, c(3, 1, 2))

  while (max(abs(phix1 - phix0)) > 0.01) {
    phix[length(phix)] <- phix[length(phix) - 1]
    phix0 <- phix1
    exp.phi0 <- exp(outer(phix0, phir, "+"))
    exp.phi0 <- replicate(nrow(dv), exp.phi0, simplify = "array") # gives an n x R x T array
    exp.phi0 <- aperm(exp.phi0, c(3, 1, 2))

    dg1 <- digamma(dv + exp.phi0)
    tg1 <- trigamma(dv + exp.phi0)
    dg2 <- digamma(exp.phi0)
    tg2 <- trigamma(exp.phi0)
    c1 <- dv + exp.phi0
    c2 <- ev * mu + exp.phi0
    lc2 <- log(c2)

    f0 <- apply(
      (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
      2,
      sum,
      na.rm = TRUE
    )
    df0 <- apply(
      ((tg1 - tg2 + 1 / exp.phi0 - 2 / c2 + c1 * 1 / (c2^2)) * exp.phi0^2) +
        (dg1 - dg2 + log(exp.phi0) + 1 - lc2 - c1 * c2^(-1)) * exp.phi0,
      2,
      sum,
      na.rm = TRUE
    )

    phix1 <- phix0 - f0 / df0

    phix1
  }
  # fit <- lm(phix1[1:7] ~ I(1:(length(phix1) - 1)))
  # phix1[length(phix1)] <- fit$coefficients[1] + fit$coefficients[2]*length(phix1)
  phix1
}

###### (3) Hessian computations baseline model ----
hessian.beta1v <- function(dv, ev, lca) {
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array") # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat <- lca$mhat
  R <- dim(dv)[3]
  na <- dim(dv)[2]
  isoyv <- isoyear(as.Date(rownames(dv))) -
    min(isoyear(as.Date(rownames(dv)))) +
    1
  isowv <- isoweek(as.Date(rownames(dv)))
  ny <- length(unique(isoyv))
  nw <- length(unique(isowv))

  # Define grid
  grid <- expand.grid(Age = 1:na, Region = 1:R)

  # Calculation
  H.beta1v <- sapply(1:nrow(grid), function(i) {
    # Values
    x = grid$Age[i]
    r = grid$Region[i]

    # Second order derivatives
    base <- -(dv[, x, r] + exp.phi[, x, r]) *
      exp.phi[, x, r] *
      ev[, x, r] *
      mhat[, x, r] *
      (ev[, x, r] * mhat[, x, r] + exp.phi[, x, r])^(-2)

    v1 <- sapply(1:R, function(s) {
      c(rep(0, (x - 1)), sum(base, na.rm = TRUE) * (r == s), rep(0, (na - x)))
    })
    v2 <- c(
      rep(0, (x - 1)),
      sum(base * lca$kappa2[isoyv, r], na.rm = TRUE),
      rep(0, (na - x))
    )
    v3 <- sapply(1:R, function(s) {
      sapply(1:ny, function(t) {
        sum(base[isoyv == t] * lca$beta2[x], na.rm = TRUE) * (r == s)
      })
    })
    v4 <- unname(c(
      rep(0, (x - 1)),
      sum(base * lca$kappa3[isowv, r], na.rm = TRUE),
      rep(0, (na - x))
    ))
    v5 <- sapply(1:R, function(s) {
      sapply(1:nw, function(w) {
        sum(base[isowv == w] * lca$beta3[x], na.rm = TRUE) * (r == s)
      })
    })

    # Return

    c(v1, v2, v3, v4, v5)
  })

  H.beta1v
}
hessian.beta2v <- function(dv, ev, lca) {
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array") # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat <- lca$mhat
  R <- dim(dv)[3]
  na <- dim(dv)[2]
  isoyv <- isoyear(as.Date(rownames(dv))) -
    min(isoyear(as.Date(rownames(dv)))) +
    1
  isowv <- isoweek(as.Date(rownames(dv)))
  ny <- length(unique(isoyv))
  nw <- length(unique(isowv))

  # Calculation
  H.beta2v <- sapply(1:ncol(dv), function(x) {
    deriv <- (dv[, x, ] -
      (dv[, x, ] + exp.phi[, x, ]) *
        (1 - (exp.phi[, x, ] / (ev[, x, ] * mhat[, x, ] + exp.phi[, x, ]))))
    base <- -(dv[, x, ] + exp.phi[, x, ]) *
      exp.phi[, x, ] *
      ev[, x, ] *
      mhat[, x, ] *
      (ev[, x, ] * mhat[, x, ] + exp.phi[, x, ])^(-2)
    base <- sweep(base, c(1, 2), lca$kappa2[isoyv, ], '*')

    v1 <- sapply(1:R, function(r) {
      unname(c(rep(0, (x - 1)), sum(base[, r], na.rm = TRUE), rep(0, (na - x))))
    })
    v2 <- unname(c(
      rep(0, (x - 1)),
      sum(sweep(base, c(1, 2), lca$kappa2[isoyv, ], '*'), na.rm = TRUE),
      rep(0, (na - x))
    ))
    v3 <- t(sapply(1:ny, function(t) {
      colSums(base[isoyv == t, ] * lca$beta2[x], na.rm = TRUE) +
        colSums(deriv[isoyv == t, ], na.rm = TRUE)
    }))
    v4 <- unname(c(
      rep(0, (x - 1)),
      sum(sweep(base, c(1, 2), lca$kappa3[isowv, ], '*'), na.rm = TRUE),
      rep(0, (na - x))
    ))
    v5 <- t(sapply(1:nw, function(w) {
      colSums(base[isowv == w, ] * lca$beta3[x], na.rm = TRUE)
    }))

    c(v1, v2, v3, v4, v5)
  })

  H.beta2v
}
hessian.kappa2v <- function(dv, ev, lca) {
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array") # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat <- lca$mhat
  R <- dim(dv)[3]
  na <- dim(dv)[2]
  isoyv <- isoyear(as.Date(rownames(dv))) -
    min(isoyear(as.Date(rownames(dv)))) +
    1
  isowv <- isoweek(as.Date(rownames(dv)))
  ny <- length(unique(isoyv))
  nw <- length(unique(isowv))

  # Define grid
  grid <- expand.grid(Year = 1:length(unique(isoyv)), Region = 1:R)

  # Calculation
  H.kappa2v <- sapply(1:nrow(grid), function(j) {
    # Grid values
    t <- grid$Year[j]
    r <- grid$Region[j]

    # First order
    deriv <- (dv[isoyv == t, , r] -
      (dv[isoyv == t, , r] + exp.phi[isoyv == t, , r]) *
        (1 -
          exp.phi[isoyv == t, , r] /
            (ev[isoyv == t, , r] *
              mhat[isoyv == t, , r] +
              exp.phi[isoyv == t, , r])))

    # Second order derivatives
    base <- -(dv[isoyv == t, , r] + exp.phi[isoyv == t, , r]) *
      exp.phi[isoyv == t, , r] *
      ev[isoyv == t, , r] *
      mhat[isoyv == t, , r] *
      (ev[isoyv == t, , r] *
        mhat[isoyv == t, , r] +
        exp.phi[isoyv == t, , r])^(-2)
    base <- sweep(base, 2, lca$beta2, '*')

    v1 <- sapply(1:R, function(s) {
      colSums(base, na.rm = TRUE) * (r == s)
    })
    v2 <- colSums(base * lca$kappa2[t, r], na.rm = TRUE) +
      colSums(deriv, na.rm = TRUE)
    v3 <- sapply(1:R, function(s) {
      unname(c(
        rep(0, (t - 1)),
        sum(sweep(base, 2, lca$beta2, '*'), na.rm = TRUE) * (r == s),
        rep(0, (nrow(lca$kappa2) - t))
      ))
    })
    v4 <- colSums(
      sweep(base, 1, lca$kappa3[1:nrow(base), r], '*'),
      na.rm = TRUE
    )
    v5 <- sapply(1:R, function(s) {
      c(
        rowSums(sweep(base, 2, lca$beta3, '*'), na.rm = TRUE) * (r == s),
        rep(0, nw - nrow(base))
      )
    })

    # Return
    c(v1, v2, v3, v4, v5)
  })

  H.kappa2v
}
hessian.beta3v <- function(dv, ev, lca) {
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array") # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat <- lca$mhat
  R <- dim(dv)[3]
  na <- dim(dv)[2]
  isoyv <- isoyear(as.Date(rownames(dv))) -
    min(isoyear(as.Date(rownames(dv)))) +
    1
  isowv <- isoweek(as.Date(rownames(dv)))
  ny <- length(unique(isoyv))
  nw <- length(unique(isowv))

  # Calculation
  H.beta3v <- sapply(1:ncol(dv), function(x) {
    # First order
    deriv <- (dv[, x, ] -
      (dv[, x, ] + exp.phi[, x, ]) *
        (1 - exp.phi[, x, ] / (ev[, x, ] * mhat[, x, ] + exp.phi[, x, ])))

    # Second order derivatives
    base <- -(dv[, x, ] + exp.phi[, x, ]) *
      exp.phi[, x, ] *
      ev[, x, ] *
      mhat[, x, ] *
      (ev[, x, ] * mhat[, x, ] + exp.phi[, x, ])^(-2)
    base <- sweep(base, c(1, 2), lca$kappa3[isowv, ], '*')

    v1 <- sapply(1:R, function(r) {
      unname(c(rep(0, (x - 1)), sum(base[, r], na.rm = TRUE), rep(0, (na - x))))
    })
    v2 <- unname(c(
      rep(0, (x - 1)),
      sum(sweep(base, c(1, 2), lca$kappa2[isoyv, ], '*'), na.rm = TRUE),
      rep(0, (na - x))
    ))
    v3 <- t(sapply(1:nrow(lca$kappa2), function(t) {
      colSums(base[isoyv == t, ] * lca$beta2[x], na.rm = TRUE)
    }))
    v4 <- unname(c(
      rep(0, (x - 1)),
      sum(sweep(base, c(1, 2), lca$kappa3[isowv, ], '*'), na.rm = TRUE),
      rep(0, (na - x))
    ))
    v5 <- t(sapply(1:nrow(lca$kappa3), function(w) {
      colSums(base[isowv == w, ] * lca$beta3[x], na.rm = TRUE) +
        colSums(deriv[isowv == w, ], na.rm = TRUE)
    }))

    # Return
    c(v1, v2, v3, v4, v5)
  })

  H.beta3v
}
hessian.kappa3v <- function(dv, ev, lca) {
  # Presets
  exp.phi <- exp(outer(lca$phix, lca$phir, "+"))
  exp.phi <- replicate(nrow(dv), exp.phi, simplify = "array") # gives an n x R x T array
  exp.phi <- aperm(exp.phi, c(3, 1, 2))
  mhat <- lca$mhat
  R <- dim(dv)[3]
  na <- dim(dv)[2]
  isoyv <- isoyear(as.Date(rownames(dv))) -
    min(isoyear(as.Date(rownames(dv)))) +
    1
  isowv <- isoweek(as.Date(rownames(dv)))
  ny <- length(unique(isoyv))
  nw <- length(unique(isowv))

  # Define grid
  grid <- expand.grid(Week = 1:length(unique(isowv)), Region = 1:R)

  # Calculation
  H.kappa3v <- sapply(1:nrow(grid), function(j) {
    # Grid values
    w <- grid$Week[j]
    r <- grid$Region[j]

    # First order
    deriv <- (dv[isowv == w, , r] -
      (dv[isowv == w, , r] + exp.phi[isowv == w, , r]) *
        (1 -
          exp.phi[isowv == w, , r] /
            (ev[isowv == w, , r] *
              mhat[isowv == w, , r] +
              exp.phi[isowv == w, , r])))

    # Calculation
    base <- -(dv[isowv == w, , r] + exp.phi[isowv == w, , r]) *
      exp.phi[isowv == w, , r] *
      ev[isowv == w, , r] *
      mhat[isowv == w, , r] *
      (ev[isowv == w, , r] *
        mhat[isowv == w, , r] +
        exp.phi[isowv == w, , r])^(-2)
    base <- sweep(base, 2, lca$beta3, '*')

    v1 <- sapply(1:R, function(s) {
      colSums(base, na.rm = TRUE) * (r == s)
    })
    v2 <- colSums(
      sweep(
        base,
        1,
        lca$kappa2[
          isoyear(as.Date(rownames(base))) -
            min(isoyear(as.Date(rownames(dv)))) +
            1,
          r
        ],
        '*'
      ),
      na.rm = TRUE
    )
    vv3 <- sapply(1:R, function(s) {
      c(rowSums(sweep(base, 2, lca$beta2, '*'), na.rm = TRUE) * (r == s))
    })
    dimnames(vv3) <- list(
      isoyear(as.Date(rownames(vv3))) -
        min(isoyear(as.Date(rownames(dv)))) +
        1,
      dimnames(dv)[[3]]
    )

    v3 <- matrix(
      0,
      nrow = length(unique(isoyv)),
      ncol = R,
      dimnames = list(unique(isoyv), dimnames(dv)[[3]])
    )
    v3[rownames(vv3), colnames(vv3)] <- vv3

    v4 <- colSums(base * lca$kappa3[w, r], na.rm = TRUE) +
      colSums(deriv, na.rm = TRUE)
    v5 <- sapply(1:R, function(s) {
      unname(c(
        rep(0, (w - 1)),
        sum(sweep(base, 2, lca$beta3, '*'), na.rm = TRUE) * (r == s),
        rep(0, (nrow(lca$kappa3) - w))
      ))
    })

    # Return
    c(v1, v2, v3, v4, v5)
  })
}
