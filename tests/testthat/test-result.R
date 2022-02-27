x <- c(1, 45, 5, 6, 6, 23, 12, 4, 56, 546, 45)
y <- c(12, 4, 6, 6, 12, 34, 634, 6, 56, 43, 45)

# Results should be the same although dccpp::dcor is faster
mine <- dcor(x, y)
expect_true(round(mine, 8) == 0.05947956)