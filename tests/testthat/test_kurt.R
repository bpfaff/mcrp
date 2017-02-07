context("kurt")
library(mcrp)
data(MultiAsset)
MA <- as.timeSeries(MultiAsset[, 1:4])
r <- na.omit(diff(log(MA)) * 100)
N <- ncol(r)
w <- rep(1 / N, N) ## equal weight allocation

test_that("Kurtosis is matrix (N x N^3)", {
    a <- M4(r)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, N ^ 3))
})

test_that("Portfolio kurtosis is scalar and positive", {
    a <- pm4(r, w)
    expect_identical(length(a), 1L)
    expect_true(a > 0)
})

test_that("Partial derivatives of portfolio kurtosis is matrix", {
    a <- dm4(r, w)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, 1))
    b <- PortKurtDeriv(r, w)
    expect_true(is.matrix(b))
    expect_equal(dim(b), c(N, 1))
})

test_that("Kurtosis contributions sum to one or
are equal to portfolio kurtosis", {
    a <- PortKurtContrib(r, w, percentage = TRUE)
    expect_true(is.matrix(a))
    expect_equal(dim(a), c(N, 1))
    expect_equal(sum(a), 1.0)
    b <- PortKurtContrib(r, w, percentage = FALSE)
    pkurt1 <- PortKurt(r, w)
    expect_true(is.matrix(b))
    expect_equal(dim(b), c(N, 1))
    expect_equal(sum(b), pkurt1)
    d <- cm4(r, w, percentage = TRUE)
    expect_true(is.matrix(d))
    expect_equal(dim(d), c(N, 1))
    expect_equal(sum(d), 1.0)
    e <- cm4(r, w, percentage = FALSE)
    pkurt2 <- pm4(r, w)
    expect_true(is.matrix(e))
    expect_equal(dim(e), c(N, 1))
    expect_equal(sum(e), pkurt2)
})
