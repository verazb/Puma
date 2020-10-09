radiusatet <- function (y, d, pscore, xadd = NULL, Rad = 3, quantile = 0.9, 
    biascorr = 1, ytype = 1, scoreweight = 1) 
{
    xx <- cbind(pscore, xadd)
    n <- length(pscore)
    xx <- as.matrix(xx)
    xxerr = 0
    if (dim(xx)[2] > 1) {
        index <- c(1:dim(xx)[1])
        x1 <- xx[d == 1, ]
        x1 <- as.matrix(x1)
        index1 <- index[d == 1]
        x0 <- xx[d == 0, ]
        x0 <- as.matrix(x0)
        index0 <- index[d == 0]
        n1 <- length(index1)
        n0 <- length(index0)
        k1 <- dim(x1)[2]
        maxdist <- 0
        indexcontrol = c(rep(0, n1))
        x_match = matrix(0, n1, k1)
        check <- try(solve(cov(x0)), silent = TRUE)
        if (class(check) == "try-error") 
            err = 1
        else err = 0
        if (err == 0) {
            invcovmat = solve(cov(x0))
            invcovmat[1] = invcovmat[1] * scoreweight
            mindist <- c()
            for (i in 1:n1) {
                dif <- c()
                for (j in 1:k1) dif <- cbind(dif, x1[i, j] - 
                  x0[, j])
                tempo <- (dif %*% invcovmat) * dif
                mdist <- rowSums(tempo)
                temp <- seq(1:n0)
                zv0 <- temp[mdist == min(mdist)]
                zv0 <- zv0[1]
                indexcontrol[i] <- index0[zv0]
                x_match[i, ] <- x0[zv0, ]
                if (mdist[zv0] > maxdist) 
                  maxdist <- mdist[zv0]
                mindist = c(mindist, mdist[zv0])
            }
            if (quantile > 0 & quantile < 1) 
                maxdist <- quantile(mindist, quantile)
            if (quantile < 0 | quantile > 1) 
                print("Argument quantile must be larger than 0 and smaller than 1. Maximum distance is taken instead.")
            maxdist = maxdist * Rad
            wt <- c()
            matchindex0 <- c()
            for (i in 1:n1) {
                dif <- c()
                for (j in 1:k1) dif <- cbind(dif, x1[i, j] - 
                  x0[, j])
                tempo <- (dif %*% invcovmat) * dif
                mdist2 <- rowSums(tempo)
                dist <- cbind(mdist2, index0)
                e <- dist[, 1] <= maxdist
                if (sum(e) < 1.5) {
                  index00 <- indexcontrol[i]
                  w00 <- 1
                }
                if (sum(e) > 1.5) {
                  dist <- dist[e == 1, ]
                  dist[, 1] <- dist[, 1] * (dist[, 1] > 0) + 
                    (dist[, 1] <= 0) * 1e-10
                  if (maxdist == 0) 
                    maxdist <- 1e-10
                  dist[, 1] <- dist[, 1]/maxdist
                  w00 <- 1/dist[, 1]
                  w00 <- w00/sum(w00)
                  index00 = dist[, 2]
                }
                wt <- c(wt, w00)
                matchindex0 <- c(matchindex0, index00)
            }
            weights <- c(rep(0, n))
            for (i in 1:n) {
                e <- (matchindex0 == index[i]) * 1
                if (sum(e) > 0) 
                  weights[i] = t(wt) %*% e
            }
            wgt00 <- weights[d == 0]
            weights <- weights/sum(weights) + d * (1/n1)
            weightsneu <- weights
            y1 <- d * y * weights
            y0 <- (1 - d) * y * weights
            effect.unadj <- sum(y1) - sum(y0)
            effect <- effect.unadj
            bias <- 0
            if (biascorr == 1) {
                y11 = y[d == 1]
                y00 = y[d == 0]
                x0 <- cbind(x0, x0[, 1]^2)
                x1 <- cbind(x1, x1[, 1]^2)
                x00 <- cbind(1, x0)
                x11 <- cbind(1, x1)
                y_p0 = c()
                y_p1 = c()
                bincheck = 0
                if (ytype != 1) {
                  bincheck = sum(y00 == 1) + sum(y00 == 0) - 
                    length(y00)
                  if (bincheck != 0) {
                    print("Outcome is not binary. Resuming with linear bias correction as if ytype==1.")
                  }
                  if (bincheck == 0) {
                    b <- glm(y00 ~ x0, family = binomial(logit), 
                      weights = wgt00)$coef
                    y_p0 <- 1/(1 + exp(-x00 %*% b))
                    y_p1 <- 1/(1 + exp(-x11 %*% b))
                  }
                }
                if (ytype == 1 | bincheck != 0) {
                  b <- lm(y00 ~ x0, weights = wgt00)$coef
                  y_p0 <- x00 %*% b
                  y_p1 <- x11 %*% b
                  wgt000 = wgt00/sum(wgt00)
                  B = sqrt(wgt000) * x00
                  B = (t(B)) %*% B
                  check <- try(solve(B), silent = TRUE)
                  if (class(check) == "try-error") 
                    errB = 1
                  else errB = 0
                  if (errB == 0) {
                    B = solve(B)
                    mean.x11 <- c()
                    for (j in 1:(dim(x11)[2])) mean.x11 <- c(mean.x11, 
                      mean(x11[, j]))
                    wneu0 = t((t(wgt000) %*% x00) %*% B %*% (t(wgt000 * 
                      x00)))
                    wneu1 = t(mean.x11 %*% B %*% (t(wgt000 * 
                      x00)))
                    indi <- seq(1:n)[d == 0]
                    weightsneu <- c(rep(1/n1, n))
                    weightsneu[indi] <- wgt000 - wneu0 + wneu1
                  }
                }
                if ((length(y_p1) > 1) & (length(y_p0) > 1)) {
                  y_p1 <- mean(y_p1)
                  y_p0 <- weighted.mean(x = y_p0, w = wgt00)
                  bias <- y_p1 - y_p0
                  effect <- effect.unadj - bias
                }
                else {
                  effect <- effect.unadj
                  print("Something is wrong with the bias correction. Resuming without bias correction.")
                  weightsneu <- weights
                }
            }
        }
        else {
            print("Covariance matrix of variables in Mahalanobis distance is not invertable. Resuming with propensity score matching instead of Mahalanobis matching.")
            xx <- pscore
            xx <- as.matrix(xx)
        }
    }
    if (dim(xx)[2] == 1) {
        index <- c(1:length(xx))
        x1 <- xx[d == 1]
        x1 <- as.matrix(x1)
        index1 <- index[d == 1]
        x0 <- xx[d == 0]
        x0 <- as.matrix(x0)
        index0 <- index[d == 0]
        n1 <- length(index1)
        n0 <- length(index0)
        k1 <- 1
        maxdist <- 0
        indexcontrol = c(rep(0, n1))
        x_match = matrix(0, n1, k1)
        mindist = c()
        for (i in 1:n1) {
            dif <- c()
            for (j in 1:k1) dif <- cbind(dif, x1[i, j] - x0[, 
                j])
            mdist <- dif^2
            temp <- seq(1:n0)
            zv0 <- temp[mdist == min(mdist)]
            zv0 <- zv0[1]
            indexcontrol[i] <- index0[zv0]
            x_match[i] <- x0[zv0]
            if (mdist[zv0] > maxdist) 
                maxdist <- mdist[zv0]
            mindist = c(mindist, mdist[zv0])
        }
        if (quantile > 0 & quantile < 1) 
            maxdist <- quantile(mindist, quantile)
        if (quantile < 0 | quantile > 1) 
            print("Argument quantile must be larger than 0 and smaller than 1. Maximum distance is taken instead.")
        maxdist <- maxdist * Rad
        wt <- c()
        matchindex0 <- c()
        for (i in 1:n1) {
            dif <- c()
            for (j in 1:k1) dif <- cbind(dif, x1[i, j] - x0[, 
                j])
            mdist2 <- dif^2
            dist <- cbind(mdist2, index0)
            e <- dist[, 1] <= maxdist
            if (sum(e) < 1.5) {
                index00 <- indexcontrol[i]
                w00 <- 1
            }
            if (sum(e) > 1.5) {
                dist <- dist[e == 1, ]
                dist[, 1] <- dist[, 1] * (dist[, 1] > 0) + (dist[, 
                  1] <= 0) * 1e-10
                if (maxdist == 0) 
                  maxdist <- 1e-10
                dist[, 1] <- dist[, 1]/maxdist
                w00 <- 1/dist[, 1]
                w00 <- w00/sum(w00)
                index00 = dist[, 2]
            }
            wt <- c(wt, w00)
            matchindex0 <- c(matchindex0, index00)
        }
        weights <- c(rep(0, n))
        for (i in 1:n) {
            e <- (matchindex0 == index[i]) * 1
            if (sum(e) > 0) 
                weights[i] = t(wt) %*% e
        }
        wgt00 <- weights[d == 0]
        weights <- weights/sum(weights) + d * (1/n1)
        weightsneu <- weights
        y1 <- d * y * weights
        y0 <- (1 - d) * y * weights
        effect.unadj <- sum(y1) - sum(y0)
        effect <- effect.unadj
        bias <- 0
        if (biascorr == 1) {
            y11 = y[d == 1]
            y00 = y[d == 0]
            x0 <- cbind(x0, x0[, 1]^2)
            x1 <- cbind(x1, x1[, 1]^2)
            x00 <- cbind(1, x0)
            x11 <- cbind(1, x1)
            y_p0 = c()
            y_p1 = c()
            bincheck = 0
            if (ytype != 1) {
                bincheck = sum(y00 == 1) + sum(y00 == 0) - length(y00)
                if (bincheck != 0) {
                  print("Outcome is not binary. Resuming with linear bias correction as if ytype==1.")
                }
                if (bincheck == 0) {
                  b <- glm(y00 ~ x0, family = binomial(logit), 
                    weights = wgt00)$coef
                  y_p0 <- 1/(1 + exp(-x00 %*% b))
                  y_p1 <- 1/(1 + exp(-x11 %*% b))
                }
            }
            if (ytype == 1 | bincheck != 0) {
                b <- lm(y00 ~ x0, weights = wgt00)$coef
                y_p0 <- x00 %*% b
                y_p1 <- x11 %*% b
                wgt000 = wgt00/sum(wgt00)
                B = sqrt(wgt000) * x00
                B = (t(B)) %*% B
                check <- try(solve(B), silent = TRUE)
                if (class(check) == "try-error") 
                  errB = 1
                else errB = 0
                if (errB == 0) {
                  B = solve(B)
                  mean.x11 <- c()
                  for (j in 1:(dim(x11)[2])) mean.x11 <- c(mean.x11, 
                    mean(x11[, j]))
                  wneu0 = t((t(wgt000) %*% x00) %*% B %*% (t(wgt000 * 
                    x00)))
                  wneu1 = t(mean.x11 %*% B %*% (t(wgt000 * x00)))
                  indi <- seq(1:n)[d == 0]
                  weightsneu <- c(rep(1/n1, n))
                  weightsneu[indi] <- wgt000 - wneu0 + wneu1
                }
            }
            if ((length(y_p1) > 1) & (length(y_p0) > 1)) {
                y_p1 <- mean(y_p1)
                y_p0 <- weighted.mean(x = y_p0, w = wgt00)
                bias <- y_p1 - y_p0
                effect <- effect.unadj - bias
            }
            else {
                effect <- effect.unadj
                print("Something is wrong with the bias correction. Resuming without bias correction.")
                weightsneu <- weights
            }
        }
    }
    list(distances = mindist, radiusdist = maxdist, atet = effect, 
        atet.unadj = effect.unadj, weights = weightsneu, maxdist = maxdist, 
        y1 = sum(y1), y0 = sum(y0), y0cor = sum(y0) + bias)
}