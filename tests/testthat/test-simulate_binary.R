set.seed(10)
test_that("power_compare function: binary case", {
  expect_equal(power_compare(family="binary",Niter=1, sgdata=sgdata, bmrpars=log(BMR),
                             Nsample=50, betaf0=0,
                             beta_gc=c(1,Fe), par=c(0.8,0.2),
                             hotspot=c(0.9992105264,0.0007894736,0.9995310769,
                                       0.0004689231, 0.5935003662,0.4064996338,
                                       0.1473376680, 1.0000000000,
                                       0, 1.0000000000)),testdataBinary)

})


test_that("power_compare function: continuous case", {
  expect_equal(power_compare(family="cont",Niter=1, sgdata=sgdata, bmrpars=log(BMR),
                             Nsample=50, betaf0=0,
                             beta_gc=c(1,Fe), par=c(0,5,0.1,1),
                             hotspot=c(0.9992105264,0.0007894736,0.9995310769,
                                       0.0004689231, 0.5935003662,0.4064996338,
                                       0.1473376680, 1.0000000000,
                                       0, 1.0000000000)),testdataCont)

})






