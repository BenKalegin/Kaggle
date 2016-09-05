library(forecast)
x <- c(1, 1, 2, 3, 5, 8, 13, 21)
fit <- nnetar(x, size = 2, repeats = 300, )
fcast <- forecast(fit)
fcast