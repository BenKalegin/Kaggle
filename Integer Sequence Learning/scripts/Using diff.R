#library(inline)
#library(Rcpp)
#fast_equal = cxxfunction(signature(x = 'numeric', y = 'numeric'), '
#  NumericVector var(x);
#  double precision = as<double>(y);
#
#  for (int i = 0, size = var.size(); i < size; ++i) {
#    if (var[i] - var[0] > precision || var[0] - var[i] > precision)
#      return Rcpp::wrap(false);
#  }

  #return Rcpp::wrap(true);', plugin = 'Rcpp')

equals <- function(x) { diff(range(x)) < .Machine$double.eps ^ 0.5 }

ratio <- function(x) {
    result <- vector(mode = "numeric", length = length(x) - 1)
    for (i in 1:length(x) - 1)
        result[i] <- x[i + 1] / x[i]
    result
}


ratio2 <- function(x) {
    result <- vector(mode = "numeric", length = length(x) - 1)
    for (i in 2:(length(x) - 1))
        result[i] <- x[i + 1] / x[i-1]
    result
}

s1 <- c(0, 1, 3, 17, 75, 361, 1683, 7937, 37275, 175321, 824163, 3875057, 18218475, 85655881, 402715443, 1893393377, 8901903675, 41852858041, 196773803523, 925144274897, 4349623252875) # 20450023957801
ratio(s1)

4.701562 3 8

test$ratio2 <- lapply(test$Sequence, function(x) {
    r1 = x[length(x)] / x[length(x) - 1]
    r2 = x[length(x)- 1] / x[length(x)-2]
    r1 / r2 - 1
})

sum(abs(unlist(test$ratio2)) < 0.01, na.rm = TRUE) # 36770
