library(gmp)

solveRecurrent <- function(x, depth) {
    A <- matrix(NA, nrow = depth + 1, ncol = depth + 1, byrow = TRUE)
    b <- matrix(NA, nrow = depth + 1, ncol = 1)

    for (r in 1:(depth + 1)) {
        for (c in 1:depth) {
            A[r, c] <- x[1 + (r - 1) + (c - 1)]
        }
        A[r, depth + 1] <- 1
        b[r] <- x[1 + depth + (r - 1)]
    }
    result <- try(round(solve(A, b), 1));
    if ("matrix" != class(result)) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    }
    t(result)[1,]
}

isRecurrent <- function(x, depth, s) {
    result = TRUE
    if (!anyNA(s)) {
        for (i in (depth + 1):(length(x) - depth)) {
            value <- as.bigz(0)
            for (j in 1:depth) {
                value <- add.bigz(value, div.bigz(x[i + j - 1] * as.bigz(round(s[j] * 10)), as.bigz(10)))
            }
            value <- value + as.bigz(round(s[j + 1]))
            if (is.na(value))
                stop()
            if ( value != x[i + depth]) {
                result = FALSE
                break;
            }
        }
    } else
        result = FALSE
    result
}

predictNext <- function(x, depth, s, isRecurrent) {
    if (isRecurrent) {
        i <- length(x) - depth + 1;
        value <- as.bigz(0)
        for (j in 1:depth) {
            value <- add.bigz(value, div.bigz(x[i + j - 1] * as.bigz(round(s[j] * 10)), as.bigz(10)))
        }
        value <- value + as.bigz(round(s[j + 1]))
        if (is.na(value))
            stop()
    } else
        value <- NA
    as.character(value)
}

#system.time(rec4 <- lapply(sequences, isRecurrent, depth = 4)) # 1.7 on 10000, then 3.07 when value == x[] changed to abs(value - x) < 0.01, but records found 6120 instead of 227

test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
test$BigSequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.bigz)
test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.numeric)
# length(testseq) == 113845
#test$Solve2 <- lapply(test$Sequence, FUN = solveRecurrent, depth = 2)
#test$Check2 <- mapply(isRecurrent, test$Sequence, 2, test$Solve2)
# sum(test$Check2 == TRUE) 2510
# sum(test$Check3 == TRUE) 2586
# sum(test$Check4 == TRUE) 2505
# sum(test$Check5 == TRUE) 1637

n <- function(prefix, number) {
    paste(prefix, i, sep = "")
}

limitDepth <- 2

sink("../output/iterations.txt")
maxDepth <- limitDepth
for (i in 1:limitDepth) {
    print(paste("---- iteration ", i, " ----- ", Sys.time()))
    test[[n("Solve", i)]] <- lapply(test$Sequence, FUN = solveRecurrent, depth = i)

    test[[n("Check", i)]] <- mapply(isRecurrent, test$BigSequence, i, test[[n("Solve", i)]])
    found <- sum(test[[n("Check", i)]] == TRUE)
    print(paste("---- found ", found, " ----- ", Sys.time()))

    if (found == 0) {
        maxDepth <- i
        break;
    }
}
sink()
unlink("../output/iterations.txt")


for (i in 1:maxDepth) {
    test[[n("Last", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$BigSequence, i, test[[n("Solve", i)]], test[[n("Check", i)]])
}


coalesce1a <- function(...) {
    ans <- ..1
    for (elt in list(...)[-1]) {
        i <- which(is.na(ans))
        ans[i] <- elt[i]
    }
    ans
}

test$Last <- NA
for (cn in sapply(c(1:maxDepth), function(x) paste("Last", x, sep = ""))) {
    test$Last <- coalesce1a(test$Last, test[[cn]])
}


test$Last[is.na(test$Last)] <- 0 # 16936 for 2..62
#test$Last[abs(test$Last) < 1E9] <- as.integer(test$Last[abs(test$Last) < 1E9])
write.csv(test[, c("Id", "Last")], "../output/recurrent-exp2.csv", row.names = FALSE)

# experiment 1
#---- iteration  1  -----  2016-08-03 23:45:35
#---- found  2402  -----  2016-08-03 23:45:52
#---- iteration  2  -----  2016-08-03 23:45:52
#---- found  6739  -----  2016-08-03 23:46:11

# experiment 2: changed check and predict using bigz
#[1] "---- iteration  1  -----  2016-08-05 00:37:22"
#[1] "---- found  823  -----  2016-08-05 00:38:14"
#[1] "---- iteration  2  -----  2016-08-05 00:38:14"
#[1] "---- found  3000  -----  2016-08-05 00:39:27"
