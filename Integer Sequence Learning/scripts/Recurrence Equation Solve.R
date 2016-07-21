# x[j]   = a * x[j-1] + b * x[j-2]
# x[j+1] = a * x[j]   + b * x[j-1]

solveRecurrent <- function(x, depth, offset) {
    A <- matrix(0, nrow = depth, ncol = depth)
    b <- matrix(0, nrow = depth, ncol = 1)

    for (r in 1:depth) { 
        for (c in 1:depth) {
            A[r, c] <- x[offset + (r - 1) + (depth - c)]
        }
        b[r] <- x[offset + depth + (r-1)] 
    }
    result <- try(round(solve(A, b)));
    if ("matrix" != class(result)) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    }
    result
}


#x11 <- c(1, 1, 2, 3, 5)
#x12 <- c(1, 3, 7, 17)
#x5 <- c(18, 24, 30, 36, 42, 54)
#solveRecurrent(x5, depth = 2, offset = 3)

isRecurrent <- function(x, depth) {
    s1 <- solveRecurrent(x, depth = depth, offset = 1)
    s2 <- solveRecurrent(x, depth = depth, offset = 2)
    identical(s1, s2)
}


#isRecurrent(x11, depth = 2)
#isRecurrent(x12, depth = 2)

train <- read.csv("../input/train.csv", stringsAsFactors = FALSE, nrow = 1000)
sequences <- lapply(strsplit(train$Sequence, split = ","), FUN = as.numeric)

rec2 <- lapply(sequences, isRecurrent, depth = 2)




