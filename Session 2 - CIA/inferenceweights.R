make_var <- function (y, w, bws) {
    N = length(y)
    if ((sum(y == 0) + sum(y == 1)) == N) 
        y_binary = 1
    if ((sum(y == 0) + sum(y == 1)) != N) 
        y_binary = 0
    if (1 * (is.null(bws)) == 1) {
        intquart = quantile(w, 0.75) - quantile(w, 0.25)
        sigma = min(sd(w), (intquart/1.349))
        bws = c(2.34 * sigma * N^(-1/5))
        bws = max(bws, 1e-05)
    }
    if (y_binary == 1) {
        erw_x = fitted(npreg(y ~ w, bws = bws, ckertype = "epanechnikov"))
        var_x = erw_x * (1 - erw_x)
    }
    if (y_binary == 0) {
        erw_x = fitted(npreg(y ~ w, bws = bws, ckertype = "epanechnikov"))
        var_x = fitted(npreg((y - erw_x)^2 ~ w, bws = bws, ckertype = "epanechnikov"))
    }
    var = (mean((w^2) * var_x)) * N
    list(bandwidth = bws, variance = var)
}


inference.weights <- function (y, d, w, bw1 = NULL, bw0 = NULL, estimand = "ATET") {
    w1 = w[d == 1]
    w0 = w[d == 0]
    y1 = y[d == 1]
    y0 = y[d == 0]
    if (estimand == "ATET" | estimand == "ATENT") {
        temp0 = make_var(y0, w0, bw0)
        var0 = temp0$variance
        var1 = var(y1)/length(y1)
        bw0 <- temp0$bandwidth
    }
    if (estimand == "ATE") {
        temp1 <- make_var(y1, w1, bw1)
        var1 = temp1$variance
        temp0 <- make_var(y0, w0, bw0)
        var0 = temp0$variance
        bw1 <- temp1$bandwidth
        bw0 <- temp0$bandwidth
    }
    list(var0 = var0, var1 = var1, bw1 = bw1, bw0 = bw0)
}


bootstrap_ATET <- function (x, xadd = NULL, d, y, Rad = 3, quantile = 0.9, biascorr = biascorr,
    ytype = ytype, logit = 0, index = 0, boot = 299, effect, 
    se) 
{
    obs <- length(y)
    xx <- cbind(rep(1, obs), xadd)
    xx <- as.matrix(xx)
    x <- as.matrix(x)
    rboot1 <- c()
    rboot2 <- c()
    while (length(rboot1) < boot) {
        sboot <- sample(1:obs, obs, TRUE)
        if (dim(x)[2] == 1) 
            xb <- x[sboot]
        if (dim(x)[2] > 1) 
            xb <- x[sboot, ]
        db <- d[sboot]
        yb <- y[sboot]
        if ((dim(xx)[2]) == 1) {
            xaddb <- NULL
        }
        if ((dim(xx)[2]) == 2) {
            xaddb <- xadd[sboot]
        }
        if ((dim(xx)[2]) > 2) {
            xaddb <- xadd[sboot, ]
        }
        if (index != 1) {
            if (logit != 1) 
                pscoreb <- glm(db ~ xb, family = binomial(probit))$fitted
            if (logit == 1) 
                pscoreb <- glm(db ~ xb, family = binomial(logit))$fitted
        }
        if (index == 1) 
            pscoreb <- lm(db ~ xb)$fitted
        temp1 <- radiusatet(y = yb, d = db, pscore = pscoreb, 
            xadd = xaddb, Rad = Rad, quantile = quantile, biascorr = biascorr, 
            ytype = ytype)
        rboot1 = c(rboot1, temp1$atet)
        rboot2 = c(rboot2, temp1$atet.unadj)
    }
    se1 <- sd(rboot1)
    mean1 <- mean(rboot1)
    se2 <- sd(rboot2)
    mean2 <- mean(rboot2)
    tboot1 <- (rboot1 - effect)/se
    pval.t <- mean(abs(tboot1) > rep(abs(effect/se), length(tboot1)))
    list(se1 = se1, se2 = se2, pval.t = pval.t, rboot1 = rboot1, 
        rboot2 = rboot2)
}


bootstrap_ATE  <- function (x, xadd = NULL, d, y, Rad = 3, quantile = 0.9, biascorr = biascorr,
    ytype = ytype, logit = 0, index = 0, boot = 299, effect, 
    se) 
{
    obs <- length(y)
    xx <- cbind(rep(1, obs), xadd)
    xx <- as.matrix(xx)
    x <- as.matrix(x)
    rboot1 <- c()
    rboot2 <- c()
    while (length(rboot1) < boot) {
        sboot <- sample(1:obs, obs, TRUE)
        if (dim(x)[2] == 1) 
            xb <- x[sboot]
        if (dim(x)[2] > 1) 
            xb <- x[sboot, ]
        db <- d[sboot]
        yb <- y[sboot]
        if ((dim(xx)[2]) == 1) {
            xaddb <- NULL
        }
        if ((dim(xx)[2]) == 2) {
            xaddb <- xadd[sboot]
        }
        if ((dim(xx)[2]) > 2) {
            xaddb <- xadd[sboot, ]
        }
        if (index != 1) {
            if (logit != 1) 
                pscoreb <- glm(db ~ xb, family = binomial(probit))$fitted
            if (logit == 1) 
                pscoreb <- glm(db ~ xb, family = binomial(logit))$fitted
        }
        if (index == 1) 
            pscoreb <- lm(db ~ xb)$fitted
        temp1 <- radiusatet(y = yb, d = db, pscore = pscoreb, 
            xadd = xaddb, Rad = Rad, quantile = quantile, biascorr = biascorr, 
            ytype = ytype)
        notreat = 1 - db
        yyb = -yb
        ppscoreb = 1 - pscoreb
        temp2 <- radiusatet(y = yyb, d = notreat, pscore = ppscoreb, 
            xadd = xaddb, Rad = Rad, quantile = quantile, biascorr = biascorr, 
            ytype = ytype)
        eff1 <- mean(db) * temp1$atet + mean(notreat) * (temp2$atet)
        eff2 <- mean(db) * temp1$atet.unadj + mean(notreat) * 
            (temp2$atet.unadj)
        weights <- temp1$weights + temp2$weights
        weights <- weights * db/sum(weights * db) + weights * 
            (1 - db)/sum(weights * (1 - db))
        rboot1 = c(rboot1, eff1)
        rboot2 = c(rboot2, eff2)
    }
    se1 <- sd(rboot1)
    mean1 <- mean(rboot1)
    se2 <- sd(rboot2)
    mean2 <- mean(rboot2)
    tboot1 <- (rboot1 - effect)/se
    pval.t <- mean(abs(tboot1) > rep(abs(effect/se), length(tboot1)))
    list(se1 = se1, se2 = se2, pval.t = pval.t, rboot1 = rboot1, 
        rboot2 = rboot2)
}