#' A direct copy of Leo Guelman's uplift::qini function
#' available in the R uplift package at commit
#' 95965272e71c312623c95c439fb0b84f95c185b7:
#' https://github.com/cran/uplift/blob/95965272e71c312623c95c439fb0b84f95c185b7/R/qini.R#L5
#'
#' @param p1 vector of numeric uplift predictions.  Some uplift models produce two predictions: if-treated and if-control.  In this case, if-treated predictions can be provided as p1, and if-control predictions can be provided as p0.
#' @param W vector of {0,1} treatment indicators
#' @param Y vector of {0,1} outcomes
#' @param p0 vector of numeric control predictions (default 0)
#' @param plotit boolean plot the Qini chart
#' @param direction 1: calculate the differential response as p1-p0, 2: p0-p1
#' @param groups 5, 10, or 20: the number of quantiles in which to divide the population
#'
#' @import stats
#'
#' @examples
#'
#' set.seed(0)
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p <- 3
#' beta <- -0.5
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- rl(pmax(beta+X[,1], 0) * W + X[,2])
#' p1 <- 1/(1+exp(-(beta+X[,1])))
#' plot_uplift_guelman(p1, W, Y, groups=10, plotit=TRUE)
#'
#' \donttest{
#' library(grf)
#' set.seed(123)
#'
#' alpha <- 0.1
#' n <- 1000
#' W <- rbinom(n, 1, 0.5)
#' Y <- W
#' p1 <- Y + alpha*rnorm(n)
#' plot_uplift_guelman(p1, W, Y, groups=10)
#'
#'
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p = 10
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.2)
#' Y <- rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
#' tau.forest <- causal_forest(X, Y, W)
#' tau.hat <- predict(tau.forest, X)
#' plot_uplift_guelman(tau.hat$predictions, W, Y)
#'
#' }
#' @export

plot_uplift_guelman <- function(p1,
                                W,
                                Y,
                                p0=rep(0,length(p1)),
                                plotit=TRUE,
                                direction = 1,
                                groups=10
                                ){

  performance <- function(p1, p0, Y, W,
                          direction = 1, groups = 10) {

    ### check valid arguments
    if (!direction %in% c(1, 2))
      stop("uplift: direction must be either 1 or 2")
    if (!groups %in% c(5, 10, 20))
      stop("uplift: groups must be either 5, 10 or 20")

    ### check for NAs.
    if (any(is.na(p1))) stop("uplift: NA not permitted in p1")
    if (any(is.na(p0))) stop("uplift: NA not permitted in p0")
    if (any(is.na(Y))) stop("uplift: NA not permitted in Y")
    if (any(is.na(W))) stop("uplift: NA not permitted in W")

    ### check classes
    if(!is.numeric(Y))
      stop("uplift: Y must be a numeric vector")
    if(!is.numeric(W))
      stop("uplift: W must be a numeric vector")

    ### check valid values for Y and W
    if (!all(Y %in% c(0, 1)))
      stop("uplift: Y must be either 0 or 1")
    if (!all(W %in% c(0, 1)))
      stop("uplift: W must be either 0 or 1")

    ### check length of arguments
    if (length (p1) != length(p0) |
        length (Y) != length(W)                |
        length(p1) != length(Y))
      stop("uplift: arguments p1, p0, Y and W must all have the same length")

    ### define dif.pred based on direction
    if (direction == 2) {
      dif.pred = p0 - p1} else {
        dif.pred = p1 - p0
      }


    mm <- cbind(dif.pred = dif.pred, Y = Y, W = W, dif.pred_r = rank(-dif.pred))
    bk <- unique(quantile(mm[, 4],
                          probs = seq(0, 1, 1 / groups)))
    print(bk)
    if ((length(bk)-1) != groups)
      warning("uplift: due to ties in uplift predictions, the number of groups is less than ", groups)

    mm <- cbind(mm, decile = cut(mm[, 4], breaks = bk, labels = NULL,
                                 include.lowest = TRUE))

    n.y1_ct0 <- tapply(mm[mm[, 3] == 0, ][, 2], mm[mm[, 3] == 0, ][, 5], sum)
    n.y1_ct1 <- tapply(mm[mm[, 3] == 1, ][, 2], mm[mm[, 3] == 1, ][, 5], sum)
    r.y1_ct0 <- tapply(mm[mm[, 3] == 0, ][, 2], mm[mm[, 3] == 0, ][, 5], mean)
    r.y1_ct1 <- tapply(mm[mm[, 3] == 1, ][, 2], mm[mm[, 3] == 1, ][, 5], mean)
    n.ct0 <- tapply(mm[mm[, 3] == 0, ][, 2], mm[mm[, 3] == 0, ][, 5], length)
    n.ct1 <- tapply(mm[mm[, 3] == 1, ][, 2], mm[mm[, 3] == 1, ][, 5], length)

    df <- merge(cbind(n.y1_ct0, r.y1_ct0, n.ct0), cbind(n.y1_ct1, r.y1_ct1, n.ct1), by= "row.names", all = TRUE)

    df$Row.names <- as.numeric(df$Row.names)
    df[, c(2, 4, 5, 7)][is.na(df[, c(2, 4, 5, 7)])] <- 0 # missing implies 0 counts

    if (direction == 2) {
      df$uplift = df$r.y1_ct0 - df$r.y1_ct1} else {
        df$uplift = df$r.y1_ct1 - df$r.y1_ct0
      }
    df <- df[order(df$Row.names), ]

    res <- cbind(group   = df$Row.names,
                 n.ct1    = df$n.ct1,
                 n.ct0    = df$n.ct0,
                 n.y1_ct1 = df$n.y1_ct1,
                 n.y1_ct0 = df$n.y1_ct0,
                 r.y1_ct1 = df$r.y1_ct1,
                 r.y1_ct0 = df$r.y1_ct0,
                 uplift   = df$uplift)

    res <- round(res, 6)
    class(res) <- "performance"
    return(res)
  }

  qini <- function(x, direction = 1, plotit = TRUE, ...) {

    if (!inherits(x, "performance"))
      stop("uplift: x is not of class performance")

    ### check valid arguments
    if (!direction %in% c(1, 2))
      stop("uplift: direction must be either 1 or 2")

    perf <- x
    groups <- nrow(perf)

    if (direction == 1) {

      ### Model Incremental gains
      inc.gains <- cumsum(perf[, 4] - perf[, 5] * sum(perf[, 2]) / sum(perf[, 3])) / sum(perf[, 2])

      ### Overall incremental gains
      overall.inc.gains <- sum(perf[, 4]) / sum(perf[, 2]) - sum(perf[, 5]) / sum(perf[, 3])


    } else {

      ### Model Incremental gains
      inc.gains <- cumsum(-1 * (perf[, 4] - perf[, 5] * sum(perf[, 2]) / sum(perf[, 3]))) / sum(perf[, 2])

      ### Overall incremental gains
      overall.inc.gains <- sum(perf[, 5]) / sum(perf[, 3]) - sum(perf[, 4]) / sum(perf[, 2])

    }

    ### Random incremental gains
    random.inc.gains <- cumsum(rep(overall.inc.gains / groups, groups))

    ### Compute area under the model incremental gains (uplift) curve
    x <- seq(1 / groups, 1, 1 / groups)
    y <- inc.gains

    auuc <- 0
    for (i in 2:length(x)) {
      auuc <- auuc + 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
    }

    ### Compute area under the random incremental gains curve
    y.rand <- random.inc.gains
    auuc.rand <- 0
    for (i in 2:length(x)) {
      auuc.rand <- auuc.rand + 0.5 * (x[i] - x[i-1]) * (y.rand[i] + y.rand[i-1])
    }

    ### Compute the difference between the areas (Qini coefficient)
    Qini <- auuc - auuc.rand
    miny <- 100 * min(c(random.inc.gains, inc.gains))
    maxy <- 100 * max(c(random.inc.gains, inc.gains))

    if (plotit) {
      plot(inc.gains * 100 ~ seq(100 / groups, 100, 100 / groups), type ="b",
           col = "blue", lty = 2, xlab = "Proportion of population targeted (%)",
           ylab = "Cumulative incremental gains (pc pt)", ylim = c(miny, maxy), ...)
      lines(random.inc.gains * 100 ~ seq(100 / groups, 100, 100 / groups), type = "l", col = "red", lty = 1)
      legend("topright", c("Model", "Random"),
             col=c("blue", "red"), lty=c(2,1))
    }

    res <- list(Qini = Qini,
                inc.gains = inc.gains,
                random.inc.gains = random.inc.gains)

    return(res)

  }

  perf_dr <- performance(p1, p0, Y, W, direction=direction, groups)
  Q <- qini(perf_dr, plotit=plotit)
  return(Q)

}
