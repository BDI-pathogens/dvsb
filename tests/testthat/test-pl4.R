test_that("PL4_inv(PL4(x)) == x", {
  f <- c(0.95, 0.01, 3.33, 2.44)
  xlog_values <- -1:1
  for (xlog in xlog_values) {
    expect_equal(xlog,
                 xlog %>%
                   PL4(f[1], f[2], f[3], f[4]) %>%
                   PL4_inv(f[1], f[2], f[3], f[4]))
  }
})

test_that("PL4(PL4_inv(y)) == y", {
  f <- c(0.95, 0.01, 3.33, 2.44)
  y_values <- f[4] + (f[3] - f[4]) * 1:3 / 4
  for (y in y_values) {
    expect_equal(y,
                 y %>%
                   PL4_inv(f[1], f[2], f[3], f[4]) %>%
                   PL4(f[1], f[2], f[3], f[4]))
  }
})
