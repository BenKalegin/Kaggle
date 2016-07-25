# x[j]   = a * x[j-1] + b * x[j-2]
# x[j+1] = a * x[j]   + b * x[j-1]

solveRecurrent <- function(x, depth) {
    if (length(x) < depth + 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        A <- matrix(0, nrow = depth, ncol = depth)
        b <- matrix(0, nrow = depth, ncol = 1)

        for (r in 1:depth) {
            for (c in 1:depth) {
                A[r, c] <- x[1 + (r - 1) + (c - 1)]
            }
            b[r] <- x[1 + depth + (r - 1)]
        }
        result <- try(round(solve(A, b)));
        if ("matrix" != class(result)) {
            result <- matrix(NA, nrow = depth, ncol = 1);
        }
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

test$Solve3 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 3)
test$Check3 <- mapply(isRecurrent, test$Sequence, 3, test$Solve3)
# sum(test$Check3 == TRUE) 2586
# sum(test$Check3 == TRUE & test$Check2 == FALSE) 2551
# test[test$Check3 == TRUE & test$Check2 == TRUE,]

test$Solve4 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 4)
test$Check4 <- mapply(isRecurrent, test$Sequence, 4, test$Solve4)
# sum(test$Check4 == TRUE) 2505
# sum(test$Check4 == TRUE & test$Check2 == FALSE & test$Check3 == FALSE) 2435

test$Solve5 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 5)
test$Check5 <- mapply(isRecurrent, test$Sequence, 5, test$Solve5)
# sum(test$Check5 == TRUE) 1637
# sum(test$Check5 == TRUE & test$Check4 == FALSE & test$Check2 == FALSE & test$Check3 == FALSE) 1578

test$Solve6 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 6)
test$Check6 <- mapply(isRecurrent, test$Sequence, 6, test$Solve6)
# sum(test$Check6 == TRUE) 1351
# sum(test$Check6 == TRUE & test$Check5 == FALSE & test$Check4 == FALSE & test$Check2 == FALSE & test$Check3 == FALSE) 1293

test$Solve7 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 7)
test$Check7 <- mapply(isRecurrent, test$Sequence, 7, test$Solve7)
# sum(test$Check7 == TRUE) 884
# sum(test$Check7 == TRUE & test$Check6 == FALSE & test$Check5 == FALSE & test$Check4 == FALSE & test$Check2 == FALSE & test$Check3 == FALSE) 840

n <- function(prefix, number) {
    paste(prefix, i, sep = "")
}

sink("../output/iterations.txt")
for (i in 62:100) {
    print(paste("---- iteration ", i, " -----"))
    test[[n("Solve", i)]] <- lapply(test$Sequence, FUN = solveRecurrent, depth = i)

    test[[n("Check", i)]] <- mapply(isRecurrent, test$Sequence, i, test[[n("Solve", i)]])
    found <- sum(test[[n("Check", i)]] == TRUE)
    print(paste("---- found ", found, " -----"))

    if (found == 0)
        break;
}
sink()
unlink("../output/iterations.txt")

#save(test, file = "../output/temp62.rda")
#load(file = "../output/temp62.rda")
#for (i in 8:62) {
#    names(test)[names(test) == paste("Solve", i)] <- n("Solve", i)
#    names(test)[names(test) == paste("Check", i)] <- n("Check", i)
#}

for (i in 2:62) {
    test[[n("Last", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$Sequence, i, test[[n("Solve", i)]], test[[n("Check", i)]])
}



coalesce1a <- function(...) {
    ans <- ..1
    for (elt in list(...)[-1]) {
        i <- which(is.na(ans))
        ans[i] <- elt[i]
    }
    ans
}

test$Last <-NA
for (cn in sapply(c(2:62), function(x) paste("Last", x, sep = ""))) {
    test$Last <- coalesce1a(test$Last, test[[cn]])
}

test$Last[is.na(test$Last)] <- 0 # 16936 for 2..62
test$Last[abs(test$Last) < 1E8] <- as.integer(test$Last[abs(test$Last) < 1E8])
write.csv(test[, c("Id", "Last")], "../output/recurrent62.csv", row.names = FALSE) 
# submission 1 : 0.04014 submission 2 (rec23) 0.06096, submission 3 (up to 7) 0.11148
# submission 4 2-62 0.15547

