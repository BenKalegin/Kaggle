# x[j]   = a * x[j-1] + b * x[j-2]
# x[j+1] = a * x[j]   + b * x[j-1]

solveRecurrent <- function(x, depth, offset) {
    A <- matrix(0, nrow = depth, ncol = depth)
    b <- matrix(0, nrow = depth, ncol = 1)

    for (r in 1:depth) { 
        for (c in 1:depth) {
            A[r, c] <- x[offset + (r - 1) + (c-1)]
        }
        b[r] <- x[offset + depth + (r-1)] 
    }
    result <- try(round(solve(A, b)));
    if ("matrix" != class(result)) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    }
    t(result)
}


#x11 <- c(1, 1, 2, 3, 5)
#x12 <- c(1, 3, 7, 17)
#x5 <- c(18, 24, 30, 36, 42, 54)
#solveRecurrent(x5, depth = 2, offset = 3)

isRecurrent <- function(x, depth) {
    s <- solveRecurrent(x, depth = depth, offset = 1)
    result = TRUE
    if (!is.na(s)) {
      for(i in 2:(length(x) - depth)) {
        value <- 0
        for (j in 1:depth) {
          value <- value + x[i+j] * s[j]    
        }
        if (value != x[i+depth]) {
          result = FALSE
          break;
        }
      }
    }else
      result = FALSE
    
}


#isRecurrent(x11, depth = 2)
#isRecurrent(x12, depth = 2)

train <- read.csv("../input/train.csv", stringsAsFactors = FALSE, nrow = 1000)
sequences <- lapply(strsplit(train$Sequence, split = ","), FUN = as.numeric)

rec2 <- lapply(sequences, isRecurrent, depth = 2)




