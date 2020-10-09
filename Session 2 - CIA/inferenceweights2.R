make_var2 <- function (y, w) {
    nn <- length(w)
    m = max(10, round(sqrt(length(w[w > 0])) * 2, digits = 0))
    var.w <- matrix(0, nn, 1)
    for (i in 1:nn) {
        if (w[i] > 0) {
            dif <- (w[i] - w[w > 0])^2
            dif.sort <- sort(dif, decreasing = FALSE)
            dif <- (w[i] - w)^2
            y.matched <- y[dif <= dif.sort[m + 1] & w > 0]
            var.w[i] <- var(y.matched)
        }
    }
    list(var = sum(w^2 * c(var.w)))
}

inference.weights.2 <- function (y, d, w, estimand = "ATET") 
{
    w1 = w[d == 1]
    w0 = w[d == 0]
    y1 = y[d == 1]
    y0 = y[d == 0]
    if (estimand == "ATET" | estimand == "ATENT") {
        var0 <- make_var2(y0, w0)$var
        var1 <- var(y1)/length(y1)
    }
    if (estimand == "ATE") {
        var1 <- make_var2(y1, w1)$var
        var0 <- make_var2(y0, w0)$var
    }
    list(var0 = var0, var1 = var1)
}