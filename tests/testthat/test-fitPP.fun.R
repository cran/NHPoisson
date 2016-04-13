context("test-fitPP.fun.R")

data(BarTxTn)
covB <- cbind(cos(2*pi*BarTxTn$dia/365),
              sin(2*pi*BarTxTn$dia/365),
              BarTxTn$TTx,BarTxTn$Txm31,
              BarTxTn$Txm31**2)

invisible(
  capture.output(
    BarEv <- POTevents.fun(
      T = BarTxTn$Tx, thres = 318,
      date = cbind(BarTxTn$ano, BarTxTn$mes, BarTxTn$dia))
  ))
invisible(
  capture.output(
    mod1B <- fitPP.fun(
      covariates = covB,
      posE = BarEv$Px, inddat = BarEv$inddat,
      tit ="BAR Tx; cos, sin, TTx, Txm31, Txm31**2",
      start = list(b0 = -100, b1 = 1, b2 = -1, b3 = 0, b4 = 0, b5 = 0))
  ))

test_that("fitPP.fun returns output that conforms to the CRAN library", {
  expect_output(str(mod1B), "22 slots")
  expect_output(str(mod1B), "-89.28928 2.53403 1.42455 -0.00554 0.55688")
})

