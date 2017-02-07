##
## Function for multiple criteria risk parity optimization
## with respect to higher moments
##
#' Multiple criteria risk parity optimization
#'
#' This function conducts a multiple criteria optimization with respect to the assets'
#' risk contributions of the portfolio higher moments (variance, skewness and kurtosis).
#' The three-element weight vector \code{lambda} determines with respect to which higher moments the
#' optimization is conducted. \code{NA} entry(ies) in this vector indicate the skipping of
#' a higher moment. That is \code{lambda = c(1, NA, NA)} would conduct a pure ERC optimization
#' with respect to the assets' contribution to the portfolio variance.
#'
#' @param start \code{numeric}, a \eqn{(N x 1)} vector of starting values.
#' @param returns \code{matrix}, a \eqn{(T x N)} matrix of asset returns.
#' @param lambda \code{numeric}, a \eqn{(3 x 1)} weight vector for the partial objectives.
#' @param ... ellipsis argument, passed to \code{nlminb()}.
#'
#' @return Object of S4-class \code{\link[FRAPO]{PortSol-class}}.
#'
#'
#' @references Baitinger, E. and Dragosch, A. and Topalova, A. (2017),
#' Extending the Risk Parity Approach to Higher Moments: Is there Any Value Added?,
#' \emph{The Journal of Portfolio Managament}, \bold{43}(2), 24--36.
#'
#' Boudt, K. and Peterson, B. and Croux, C. (2008/09), Estimation and decomposition of downside risk
#' for portfolios with non-normal returns, \emph{The Journal of Risk}, \bold{11}(2), Winter 2008/09, 79--103.
#'
#' Jondeau, E. and Rockinger, M. (2006), Optimal portfolio allocation under higher moments,
#' \emph{European Financial Management}, \bold{12}(1), 29--55.
#'
#' @examples
#' data(MultiAsset)
#' MA <- as.timeSeries(MultiAsset[, 1:4])
#' r <- na.omit(diff(log(MA)) * 100)
#' N <- ncol(r)
#' erc <- mcrp(start = runif(N), returns = r, lambda = c(1, NA, NA), lower = rep(0, N))
#' w <- Weights(erc)
#' PortRiskContrib(r, w)
#'
#' @export
mcrp <- function(start, returns, lambda = c(1, 1, 1), ...){
    l <- lambda
    r <- returns
    me2 <- M2(r)
    if (!is.na(l[2])) {
        me3 <- M3(r)
        }
    if (!is.na(l[3])) {
        me4 <- M4(r)
        }
    f <- function(x) {
        prisk <- c(crossprod(x, me2) %*% x)
        ans <- 0
        if (!is.na(l[1])) {
            ctb <- 2 * x * me2 %*% x
            pctb <- ctb / prisk
            ans <- ans + l[1] * stats::var(pctb)
        }
        if (!is.na(l[2])) {
            pm3 <- c(crossprod(x, me3 %*% kronecker(x, x)))
            pskew <- pm3 / (prisk ^ ( 3 / 2 ))
            dm3 <- 3 * (me3 %*% kronecker(x, x))
            term1 <- prisk ^ ( 3 / 2 ) * dm3
            term2 <- 2 * pm3 * sqrt(prisk) * me2 %*% x
            pderiv <- (term1 - term2) / prisk ^ 3
            ctb <- x * pderiv
            pctb <- ctb / pskew
            ans <- ans + l[2] * stats::var(pctb)
        }
        if (!is.na(l[3])) {
            pm4 <- c(crossprod(x, me4) %*% kronecker(kronecker(x, x), x))
            pkurt <- pm4 / prisk ^ 2
            dm4 <- 4 * (me4 %*% kronecker(kronecker(x, x), x))
            term1 <- prisk * dm4
            term2 <- 2 * pm4 * me2 %*% x
            pderiv <- (term1 - term2) / (2 * prisk ^ 3)
            ctb <- x * pderiv
            pctb <- ctb / pkurt
            ans <- ans + l[3] * stats::var(pctb)
        }
        ans
    }
    opt <- stats::nlminb(start = start, objective = f, ...)
    w <- opt$par / sum(abs(opt$par))
    idx <- !is.na(lambda)
    momchar <- c("variance", "skewness", "kurtosis")
    title <- paste("Multiple criteria objective(s):",
                   paste(momchar[idx], collapse = ", "))
    ans <- new("PortSol", weights = w,
               opt = opt,
               call = match.call(),
               type = title)
    ans
}
