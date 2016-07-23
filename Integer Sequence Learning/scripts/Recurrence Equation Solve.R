# x[j]   = a * x[j-1] + b * x[j-2]
# x[j+1] = a * x[j]   + b * x[j-1]

solveRecurrent <- function(x, depth) {
    A <- matrix(0, nrow = depth, ncol = depth)
    b <- matrix(0, nrow = depth, ncol = 1)

    for (r in 1:depth) { 
        for (c in 1:depth) {
            A[r, c] <- x[1 + (r - 1) + (c-1)]
        }
        b[r] <- x[1 + depth + (r-1)] 
    }
    result <- try(round(solve(A, b)));
    if ("matrix" != class(result)) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    }
    t(result)[1,]
}


#x11 <- c(1, 1, 2, 3, 5)
#x12 <- c(1, 3, 7, 17)
#x5 <- c(18, 24, 30, 36, 42, 54)
#solveRecurrent(x5, depth = 2, offset = 3)

isRecurrent <- function(x, depth, s) {
    result = TRUE
    if (!anyNA(s)) {
      for(i in (depth+1):(length(x) - depth)) {
        value <- 0
        for (j in 1:depth) {
          value <- value + x[i+j-1] * s[j]    
        }
        if (is.na(x[i+depth]) || abs(value - x[i+depth]) > 0.01 ) {
          result = FALSE
          break;
        }
      }
    }else
      result = FALSE
    result
}
#isRecurrent(x11, depth = 2)
#isRecurrent(x12, depth = 2)

predictNext <- function(x, depth, s, isRecurrent) {
    if (isRecurrent) {
        i <- length(x) - depth + 1;
        value <- 0
        for (j in 1:depth) {
            value <- value + x[i + j - 1] * s[j]
        }
    } else
        value <- NA
    value
}


#train <- read.csv("../input/train.csv", stringsAsFactors = FALSE, nrow = 10000)
#sequences <- lapply(strsplit(train$Sequence, split = ","), FUN = as.numeric)

#rec2 <- lapply(sequences, isRecurrent, depth = 2)
#rec3 <- lapply(sequences, isRecurrent, depth = 3)
#system.time(rec4 <- lapply(sequences, isRecurrent, depth = 4)) # 1.7 on 10000, then 3.07 when value == x[] changed to abs(value - x) < 0.01, but records found 6120 instead of 227

test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.numeric)
# length(testseq) == 113845
test$Solve2 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 2)
test$Check2 <- mapply(isRecurrent, test$Sequence, 2, test$Solve2)
# sum(test$Check2 == TRUE) 2510
# sum(test$Check2 == FALSE) 111335

#head(test[test$Check2 == TRUE,])

test$Last <- mapply(predictNext, SIMPLIFY = TRUE, test$Sequence, 2, test$Solve2, test$Check2)
test$Last[is.na(test$Last)] <- 0
test$Last[abs(test$Last) < 1E8] <- as.integer(test$Last[abs(test$Last) < 1E8])
write.csv(test[, c("Id", "Last")], "../output/recurrent2.csv", row.names = FALSE) # submission 1
