radiusmatch <- function (y, d, x, xadd = NULL, radius = 3, quantile = 0.9, biascorr = 1, 
    ynonbin = 1, estimand = "ATET", commonsup = 1, maxwgt = 1, 
    logit = 0, index = 0, boot = 299, scoreweight = 5, se.kernel = 1) 
{
    # checknp <- is.installed("np")
    # if (checknp == FALSE) 
    #     install.packages("np")
    # checknploaded <- packageLoaded("np")
    # if (checknploaded == FALSE) 
    #     library(np)
    ytype <- ynonbin
    xx <- x
    xxadd <- xadd
    lxxadd <- length(xxadd)
    ly <- length(y)
    dim <- (lxxadd > ly) * 1
    dimx <- (length(xx) > ly) * 1
    se.boot = NULL
    pval.t.boot = NULL
    if (maxwgt < 0 | maxwgt > 1) 
        print("Argument maxwgt must be larger than zero and smaller than one to be implemented. Resuming without using maxwgt.")
    if (estimand == "ATET" | estimand == "ATENT") {
        if (estimand == "ATET") {
            yy <- y
            dd <- d
        }
        if (estimand == "ATENT") {
            dd = 1 - d
            yy = y
        }
        if (index != 1) {
            if (logit != 1) 
                pscore <- glm(dd ~ xx, family = binomial(probit))$fitted
            if (logit == 1) 
                pscore <- glm(dd ~ xx, family = binomial(logit))$fitted
        }
        if (index == 1) 
            pscore = cbind(rep(1, length(dd)), xx) %*% glm(dd ~ 
                xx, family = binomial(probit))$coef
        temp = max(pscore) + 1
        if (commonsup == 1) {
            temp <- sort(pscore[dd == 0], decreasing = TRUE)
            temp = temp[1]
        }
        if (commonsup == 2) {
            temp <- sort(pscore[dd == 0], decreasing = TRUE)
            temp = temp[2]
        }
        yy = yy[pscore <= temp]
        dd = dd[pscore <= temp]
        isnull = (is.null(xadd)) * 1
        if (isnull == 0) {
            if (dim == 0) 
                xxadd = xxadd[pscore <= temp]
            if (dim == 1) 
                xxadd = xxadd[pscore <= temp, ]
        }
        if (dimx == 0) 
            xx = xx[pscore <= temp]
        if (dimx == 1) 
            xx = xx[pscore <= temp, ]
        pscore = pscore[pscore <= temp]
        if (maxwgt > 0 & maxwgt < 1) {
            if (index != 1) 
                score = pscore
            if (index == 1) {
                if (logit != 1) 
                  score <- glm(dd ~ xx, family = binomial(probit))$fitted
                if (logit == 1) 
                  score <- glm(dd ~ xx, family = binomial(logit))$fitted
            }
            if (estimand == "ATET") {
                temp = sum((1 - dd) * score/(1 - score))
                relwgt = ((1 - dd) * score/(1 - score))/temp
                maxscore = 1
                if (sum(relwgt >= maxwgt) > 0) 
                  maxscore = min(score[relwgt >= maxwgt])
                yy = yy[score <= maxscore]
                dd = dd[score <= maxscore]
                if (isnull == 0) {
                  if (dim == 0) 
                    xxadd = xxadd[score <= maxscore]
                  if (dim == 1) 
                    xxadd = xxadd[score <= maxscore, ]
                }
                if (dimx == 0) 
                  xx = xx[score <= maxscore]
                if (dimx == 1) 
                  xx = xx[score <= maxscore, ]
                pscore = pscore[score <= maxscore]
            }
            if (estimand == "ATENT") {
                temp = sum((1 - dd) * (1 - score)/score)
                relwgt = ((1 - dd) * (1 - score)/score)/temp
                minscore = 0
                if (sum(relwgt >= maxwgt) > 0) 
                  minscore = max(score[relwgt >= maxwgt])
                yy = yy[score >= minscore]
                dd = dd[score >= minscore]
                if (isnull == 0) {
                  if (dim == 0) 
                    xxadd = xxadd[score >= minscore]
                  if (dim == 1) 
                    xxadd = xxadd[score >= minscore, ]
                }
                if (dimx == 0) 
                  xx = xx[score >= minscore]
                if (dimx == 1) 
                  xx = xx[score >= minscore, ]
                pscore = pscore[score >= minscore]
            }
        }
        effect <- radiusatet(y = yy, d = dd, pscore = pscore, 
            xadd = xxadd, Rad = radius, quantile = quantile, 
            biascorr = biascorr, ytype = ytype, scoreweight = scoreweight)
        if (estimand == "ATENT") {
            effect1 <- -effect$atet
            effect2 <- -effect$atet.unadj
            y1 <- effect$y0cor
            y0 <- effect$y1
        }
        else {
            effect1 <- effect$atet
            effect2 <- effect$atet.unadj
            y1 <- effect$y1
            y0 <- effect$y0cor
        }
        weights <- effect$weights
        radiusdist = effect$radiusdist
        distances = effect$distances
        if (se.kernel == 1) 
            inf <- inference.weights(y = yy, d = dd, w = weights, 
                estimand = estimand)
        else inf <- inference.weights.2(y = yy, d = dd, w = weights, 
            estimand = estimand)
        se <- sqrt(inf$var0 + inf$var1)
        if (se.kernel == 1) 
            bw1 <- inf$bw1
        bw0 <- inf$bw0
        if (boot > 0) {
            bs <- bootstrap_ATET(x = xx, xadd = xxadd, d = dd, 
                y = yy, Rad = radius, quantile = quantile, biascorr = biascorr, 
                ytype = ytype, logit = logit, index = index, 
                boot = boot, effect = effect1, se = se)
        }
    }
    if (estimand == "ATE") {
        yy <- y
        dd <- d
        if (index != 1) {
            if (logit != 1) 
                pscore <- glm(dd ~ xx, family = binomial(probit))$fitted
            if (logit == 1) 
                pscore <- glm(dd ~ xx, family = binomial(logit))$fitted
        }
        if (index == 1) 
            pscore = cbind(rep(1, length(dd)), xx) %*% glm(dd ~ 
                xx, family = binomial(probit))$coef
        temp <- sort(pscore[dd == 0], decreasing = TRUE)
        temp = max(pscore) + 1
        if (commonsup == 1) {
            temp <- sort(pscore[dd == 0], decreasing = TRUE)
            temp = temp[1]
        }
        if (commonsup == 2) {
            temp <- sort(pscore[dd == 0], decreasing = TRUE)
            temp = temp[2]
        }
        temp2 = min(pscore) - 1
        if (commonsup == 1) {
            temp2 <- sort(pscore[dd == 1])
            temp2 = temp2[1]
        }
        if (commonsup == 2) {
            temp2 <- sort(pscore[dd == 1])
            temp2 = temp2[2]
        }
        yy = yy[pscore <= temp]
        dd = dd[pscore <= temp]
        isnull = (is.null(xadd)) * 1
        if (isnull == 0) {
            if (dim == 0) 
                xxadd = xxadd[pscore <= temp]
            if (dim == 1) 
                xxadd = xxadd[pscore <= temp, ]
        }
        if (dimx == 0) 
            xx = xx[pscore <= temp]
        if (dimx == 1) 
            xx = xx[pscore <= temp, ]
        pscore = pscore[pscore <= temp]
        yy = yy[pscore >= temp2]
        dd = dd[pscore >= temp2]
        isnull = (is.null(xadd)) * 1
        if (isnull == 0) {
            if (dim == 0) 
                xxadd = xxadd[pscore >= temp2]
            if (dim == 1) 
                xxadd = xxadd[pscore >= temp2, ]
        }
        if (dimx == 0) 
            xx = xx[pscore >= temp2]
        if (dimx == 1) 
            xx = xx[pscore >= temp2, ]
        pscore = pscore[pscore >= temp2]
        if (maxwgt > 0 & maxwgt < 1) {
            if (index != 1) 
                score = pscore
            if (index == 1) {
                if (logit != 1) 
                  score <- glm(dd ~ xx, family = binomial(probit))$fitted
                if (logit == 1) 
                  score <- glm(dd ~ xx, family = binomial(logit))$fitted
            }
            temp = sum((1 - dd) * score/(1 - score))
            relwgt = ((1 - dd) * score/(1 - score))/temp
            maxscore = 1
            if (sum(relwgt >= maxwgt) > 0) 
                maxscore = min(score[relwgt >= maxwgt])
            yy = yy[score <= maxscore]
            dd = dd[score <= maxscore]
            if (isnull == 0) {
                if (dim == 0) 
                  xxadd = xxadd[score <= maxscore]
                if (dim == 1) 
                  xxadd = xxadd[score <= maxscore, ]
            }
            if (dimx == 0) 
                xx = xx[score <= maxscore]
            if (dimx == 1) 
                xx = xx[score <= maxscore, ]
            pscore = pscore[score <= maxscore]
            temp = sum((dd) * (1 - score)/score)
            relwgt = ((dd) * (1 - score)/score)/temp
            minscore = 0
            if (sum(relwgt >= maxwgt) > 0) 
                minscore = max(score[relwgt >= maxwgt])
            yy = yy[score >= minscore]
            dd = dd[score >= minscore]
            if (isnull == 0) {
                if (dim == 0) 
                  xxadd = xxadd[score >= minscore]
                if (dim == 1) 
                  xxadd = xxadd[score >= minscore, ]
            }
            if (dimx == 0) 
                xx = xx[score >= minscore]
            if (dimx == 1) 
                xx = xx[score >= minscore, ]
            pscore = pscore[score >= minscore]
        }
        temp1 <- radiusatet(y = yy, d = dd, pscore = pscore, 
            xadd = xxadd, Rad = radius, quantile = quantile, 
            biascorr = biascorr, ytype = ytype, scoreweight = scoreweight)
        radiusdist = temp1$radiusdist
        distances = temp1$distances
        dneg = 1 - dd
        if (index != 1) 
            pscore = 1 - pscore
        if (index == 1) 
            pscore = lm(dneg ~ xx)$fitted
        temp2 <- radiusatet(y = yy, d = dneg, pscore = 1 - pscore, 
            xadd = xxadd, Rad = radius, quantile = quantile, 
            biascorr = biascorr, ytype = ytype, scoreweight = scoreweight)
        effect1 <- mean(dd) * temp1$atet + mean(1 - dd) * (-temp2$atet)
        effect2 <- mean(dd) * temp1$atet.unadj + mean(1 - dd) * 
            (-temp2$atet.unadj)
        y1 <- mean(dd) * temp1$y1 + mean(1 - dd) * (temp2$y0cor)
        y0 <- mean(dd) * temp1$y0cor + mean(1 - dd) * (temp2$y1)
        weights <- temp1$weights + temp2$weights
        weights <- weights * dd/sum(weights * dd) + weights * 
            (1 - dd)/sum(weights * (1 - dd))
        if (se.kernel == 1) 
            inf <- inference.weights(y = yy, d = dd, w = weights, 
                estimand = estimand)
        else inf <- inference.weights.2(y = yy, d = dd, w = weights, 
            estimand = estimand)
        se <- sqrt(inf$var0 + inf$var1)
        if (se.kernel == 1) 
            bw1 <- inf$bw1
        bw0 <- inf$bw0
        if (boot > 0) {
            bs <- bootstrap_ATE(x = xx, xadd = xxadd, d = dd, 
                y = yy, Rad = radius, quantile = quantile, biascorr = biascorr, 
                ytype = ytype, logit = logit, index = index, 
                boot = boot, effect = effect1, se = se)
        }
    }
    if (boot > 0) {
        se.boot = bs$se1
        pval.t.boot = bs$pval.t
    }
    dropped = length(y) - length(yy)
    list(pscore = pscore, effect = effect1, effect.unadjust = effect2, 
        se = se, se.boot = se.boot, pval.t.boot = pval.t.boot, 
        weights = weights, dropped = dropped, radius = radiusdist, 
        y1 = y1, y0 = y0)
}
