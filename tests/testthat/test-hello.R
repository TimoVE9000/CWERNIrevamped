test_that("test that testing works, by performing multiplication", {
  expect_equal(2 * 2, 4)
})



test_that("test generate_spat_abund", {
  result = length(generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)))
    expect_true(result >= 50 && result <= 1000)
})



test_that("test nulsim", {
  result = nulsim(10, 0.6, 0.1, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)), 1000)
  expect_true (is.matrix(result))
  expect_true(max (result[1,]) >= 9 && max (result[1,]) <= 10)
  expect_true(min (result[1,])==0)
})
