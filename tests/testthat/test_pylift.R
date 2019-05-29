context("placeholder tests")
library(whisker)
library(testthat)

test_that("plUpliftEval factor treatments error", {
  n <- 100;
  W <- rbinom(n, 1, 0.5)
  W <- as.factor(W)
  Y <- rbinom(n, 1, 0.5)
  p <- rnorm(n)
  expect_error(plUpliftEval(W,Y,p))
})

test_that("plUpliftEval factor outcomes error", {
  n <- 100;
  W <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  Y <- as.factor(Y)
  p <- rnorm(n)
  expect_error(plUpliftEval(W,Y,p))
})


test_that("plUpliftEval factor predictions error", {
  n <- 100;
  W <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  p <- as.factor(rnorm(n))
  expect_error(plUpliftEval(W,Y,p))
})

test_that("plUpliftEval different length treatments vs. predictions error", {
  n <- 100;
  W <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)
  p <- as.factor(rnorm(2*n))
  expect_error(plUpliftEval(W,Y,p))
})

test_that("plUpliftEval different length treatments vs. outcomes error", {
  n <- 100;
  W <- rbinom(n, 1, 0.5)
  Y <- rbinom(2*n, 1, 0.5)
  p <- as.factor(rnorm(n))
  expect_error(plUpliftEval(W,Y,p))
})

