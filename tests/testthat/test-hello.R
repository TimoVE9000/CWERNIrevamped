test_that("test that testing works, by performing multiplication", {
  expect_equal(2 * 2, 4)
})



test_that("test generate_spat_abund", {
  result = length(generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)))
    expect_true(result >= 50 && result <= 1000)
})

test_that("test nulsim", {
  result = nulsim(10, 0.6, 0.1, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 500)
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= 9 && max (result[1,]) <= 10)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
})

test_that("test CWRsim", {
  tmax=sample(7:10, 1)
  result = CWRsim(tmax,0.6 , 0.1, 0.05, 0.1, 0, 0.0005, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 500)
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= tmax-1 && max (result[1,]) <= tmax)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
})


test_that("test CWRsim", {
  tmax=sample(7:10, 1)
  result = simbvar( tmax, 0.6, 0.1, 16000, 0.05,0, 5, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 500)
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= tmax-1 && max (result[1,]) <= tmax)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
})




test_that("test nulsimSPEC", {
  result = nulsimSPEC(10, 0.6, 0.1, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 20000000, c(1,3,5))
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= 10-0.1 && max (result[1,]) <= 10)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
  expect_true (ncol (result) == 5)
  expect_true (sum(round(result[1,], 1) == c(0,1,3,5,10))==5)

})


test_that("test CWRsimSPEC", {
  result = CWRsimSPEC(10,0.6 , 0.1, 0.05, 0.1, 0, 0.0005, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),  200,c(1,3,5))
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= 10-0.1 && max (result[1,]) <= 10)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
  expect_true (ncol (result) == 5)
  expect_true (sum(round(result[1,], 1) == c(0,1,3,5,10))==5)

})


test_that("test simbvarSPEC", {
  result = simbvarSPEC( 10, 0.6, 0.1, 16000, 0.05,0, 5, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),200,c(1,3,5))
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= 10-0.1 && max (result[1,]) <= 10)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
  expect_true (ncol (result) == 5)
  expect_true (sum(round(result[1,], 1) == c(0,1,3,5,10))==5)

})



test_that("test multinulsim", {
  multinulsim(2, 0.6, 0.1, 16000,generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 2,"nulsim.tiff",  200)
  expect_true(file.exists("nulsim.tiff"))
})

test_that("test multisimbvar", {
  multisimbvar( 2, 0.6, 0.1, 16000, 0.05,0, 5, generate_spat_abund (theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 2, "bvar.tiff", 200)
  expect_true(file.exists("bvar.tiff"))
})

test_that("test multiCWRsim", {
  multiCWRsim(2,0.6 , 0.1, 0.05, 0.1, 0, 0.0005, 16000,generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 2, "CWR.tiff", 200)
  expect_true(file.exists("CWR.tiff"))
})

test_that("test multinulsimSPEC", {
multinulsimSPEC(6, 0.6, 0.1, 1600, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 2, 200, c(1,3,5))
expect_true(file.exists("1 simdata_nulsim.csv"))
expect_true(file.exists("2 simdata_nulsim.csv"))
expect_true(file.exists("1 time_rac_nulsim.tiff"))
expect_true(file.exists("3 time_rac_nulsim.tiff"))
expect_true(file.exists("5 time_rac_nulsim.tiff"))
expect_true(file.exists("5 time_rac_nulsim.tiff"))
expect_true(file.exists("Trajectory_totalcommsize_nulsim.tiff"))
})



