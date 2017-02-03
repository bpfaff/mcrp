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
#' @return \code{list}
#'
#'
#' @references Baitinger, E. and Dragosch, A. and Topalova, A. (2017),
#' Extending the Risk Parity Approach to Higher Moments: Is there Any Value Added?,
#' \emph{The Journal of Portfolio Managament}, \bold{43}(2), 24--36.
#'
#' Boudt, K. and Peterson, B. and Croux, C. (2008/09), Estimation and decomposition of downside risk
#' for portfolios with non-normal returns, \emph{The Journal of Risk}, \bold{11}(2), Winter 2008/09, 79--103.
#'
#' Jondeau, E. and Rockinger, M. (2006), Optimal portfolio allocation under higher moments, \
#' emph{European Financial Management}, \bold{12}(1), 29--55.
#'
#' @examples
#' data(EuStockMarkets)
#' r <- diff(log(EuStockMarkets)) * 100
#' N <- ncol(r)
#' w <- rep(1 / N, N) ## start values
#' erc <- mcrp(start = w, returns = r, lambda = c(1, NA, NA), lower = rep(0, N))
#' w <- Weights(erc)
#' PortRiskContrib(r, w)
#'
#' @export
mcrp <- function(start, returns, lambda = c(1, 1, 1), ...){
    f <- function(x, l = lambda, r = returns){
        ans <- 0
        if (!is.na(l[1])) {
            ans <- ans + l[1] * stats::var(PortRiskContrib(r = r, w = x))
        }
        if (!is.na(l[2])) {
            ans <- ans + l[2] * stats::var(PortSkewContrib(r = r, w = x))
        }
        if (!is.na(l[3])) {
            ans <- ans + l[3] * stats::var(PortKurtContrib(r = r, w = x))
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
