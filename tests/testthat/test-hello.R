test_that("test that testing works, by performing multiplication", {
  expect_equal(2 * 2, 4)
})



test_that("test generate_spat_abund", {
  result = length(generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)))
    expect_true(result >= 50 && result <= 1000)
})
