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

test_that("BarEv returns output that conforms to the CRAN library", {
  expect_output(str(BarEv), "214 220 222 226 687 902 1019 1163 1303 1581")
  expect_output(str(BarEv), "thres : num 318")
})