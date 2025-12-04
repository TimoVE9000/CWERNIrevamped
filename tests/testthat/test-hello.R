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
  expect_true(max (result[1,]) >= tmax-0.1 && max (result[1,]) <= tmax)
  expect_true(min (result[1,])==0)
  expect_true (nrow (result) > 1)
  expect_true (ncol (result) == 5)
  expect_true (sum(round(result[1,], 1) == c(0,1,3,5,10))==5)

})





