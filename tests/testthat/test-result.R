x <- c(1, 5, 3, 2, 4, 6, 7, 5)
y <- c(3, 5, 7, 3, 8, 4, 6, 7)

# Results should be the same although dccpp::dcor is faster
mine <- dcor(x, y)

expect_true(round(mine, 7) == 0.6288829)