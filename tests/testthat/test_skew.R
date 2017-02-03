context("skew")
library(mcrp)
data(EuStockMarkets)
r <- diff(log(EuStockMarkets)) * 100
N <- ncol(r)
w <- rep(1 / N, N) ## equal weight allocation

test_that("Co-skewness is matrix (N x N^2)", {
    a <- M3(r)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, N ^ 2))
})

test_that("Portfolio skewness is scalar", {
    a <- pm3(r, w)
    expect_identical(length(a), 1L)
})

test_that("Partial derivatives of portfolio skewness is matrix", {
    a <- dm3(r, w)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, 1))
    b <- PortSkewDeriv(r, w)
    expect_true(is.matrix(b))
    expect_equal(dim(b), c(N, 1))
})

test_that("Skewness contributions sum to one or
are equal to portfolio skewness", {
    a <- PortSkewContrib(r, w, percentage = TRUE)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, 1))
    expect_equal(sum(a), 1.0)
    b <- PortSkewContrib(r, w, percentage = FALSE)
    pskew1 <- PortSkew(r, w)
    expect_true(is.matrix(b))
    expect_equal(dim(b), c(N, 1))
    expect_equal(sum(b), pskew1)
    d <- cm3(r, w, percentage = TRUE)
    expect_true(is.matrix(d))
    expect_equal(dim(d), c(N, 1))
    expect_equal(sum(d), 1.0)
    e <- cm3(r, w, percentage = FALSE)
    pskew2 <- pm3(r, w)
    expect_true(is.matrix(e))
    expect_equal(dim(e), c(N, 1))
    expect_equal(sum(e), pskew2)
})
