set.seed(10)
test_that("power_compare function: binary case", {
  expect_equal(power_compare(family="binary",Niter=1, sgdata=sgdata, bmrpars=log(BMR),
                             Nsample=50, betaf0=0,
                             beta_gc=c(1,Fe), par=c(0.8,0.2)),testdataBinary)

})


test_that("power_compare function: continuous case", {
  expect_equal(power_compare(family="cont",Niter=1, sgdata=sgdata, bmrpars=log(BMR),
                             Nsample=50, betaf0=0,
                             beta_gc=c(1,Fe), par=c(0,5,0.1,1)),testdataCont)

})






