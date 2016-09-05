# tried both first and last approx
library(gmp)
#library(lapack, warn.conflicts = FALSE, quietly = TRUE)

solveRecurrentBias <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    if (length(x) - 2 * depth < 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        A <- matrix(NA, nrow = depth + 1, ncol = depth + 1, byrow = TRUE)
        b <- matrix(NA, nrow = depth + 1, ncol = 1)
        offset <- 1
        if (takeLast == TRUE) {
            offset <- length(x) - 2 * depth
            while ((offset > 1) & (abs(x[offset]) > 1e4))
                offset <- offset - 1
            }
        for (r in 1:(depth + 1)) {
            A[r, 1:depth] <- x[(offset + (r - 1)):(offset + (r - 1) + depth - 1)]
            #for (c in 1:depth) A[r, c] <- x[offset + (r - 1) + (c - 1)]
            A[r, depth + 1] <- 1
            b[r] <- x[offset + depth + (r - 1)]
        }
        solved <- try(solve(A, b))
        if ("matrix" != class(solved))
            result <- matrix(NA, nrow = depth, ncol = 1)
        else
            result <- round(solved, 1);
        }
    t(result)[1,]
}

solveRecurrentNobias <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    if (length(x) - 2 * depth < 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        A <- matrix(NA, nrow = depth, ncol = depth, byrow = TRUE)
        b <- matrix(NA, nrow = depth, ncol = 1)
        offset <- 1
        if (takeLast == TRUE) {
            offset <- length(x) - 2 * depth
            while ((offset > 1) & (abs(x[offset]) > 1e4))
                offset <- offset - 1
            }
        for (r in 1:depth) {
            A[r, 1:depth] <- x[(offset + (r - 1)):(offset + (r - 1) + depth - 1)]
            b[r] <- x[offset + depth + (r - 1)]
        }
        solved <- try(solve(A, b))
        if ("matrix" != class(solved))
            result <- matrix(NA, nrow = depth, ncol = 1)
        else
            result <- round(solved, 1);
        }
    t(result)[1,]
}

solveRecurrentOld <- function(x, depth) {
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

lmRecurrent <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    if (length(x) - depth < 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        if (takeLast) {
            df <- data.frame(y = tail(x, - depth))
        }
        formulaString <- "y~"
        for (i in 1:depth) {
            df[[paste0("x", i)]] <- x[i:(length(x) - depth + i - 1)]
            formulaString <- paste0(formulaString, "+x", i)
        }
        formulaString <- sub("~\\+", "~", formulaString)

        fit <- lm(formula(formulaString), df)
        maxResidual <- max(abs(fit$residuals))

        df <- list()
        for (i in 1:depth) {
            df[[paste0("x", i)]] <- x[length(x) - depth + i]
        }
        df <- as.data.frame(df)
        prediction <- predict(fit, df)

        prediction <- round(prediction)
    }
    prediction
}

isRecurrent <- function(x, depth, s, biased = FALSE) {
    result = FALSE
    if (!anyNA(s)) {
        ts <- t(head(s, depth))
        bias = ifelse(biased, tail(s, 1), 0)
        match <- 0
        notmatch <- 0
        for (i in depth:(length(x) - depth)) {
            value <- ts %*% x[i:(i + depth - 1)] + bias
            if (is.na(value))
                stop()
            if (abs(value - x[i + depth]) < abs(value) / 1E8)
                match <- match + 1
            else
                notmatch <- notmatch + 1
            }
        result <- match > notmatch * 2
    }
    result
}

predictNext <- function(x, depth, s, isRecurrent, biased = FALSE) {
    #print(id)
    s <- unlist(s)
    if (isRecurrent) {
        i <- length(x) - depth + 1;
        value <- as.bigz(0)
        for (j in 1:depth) {
            value <- add.bigz(value, div.bigz(x[i + j - 1] * as.bigz(round(s[j] * 10)), as.bigz(10)))
        }
        if (biased == TRUE)
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

limitDepth <- 70

sink("../output/iterations.txt")
maxDepth <- limitDepth
for (i in 1:limitDepth) {
    print(paste("---- iteration ", i, " ----- ", Sys.time()))
    test[[n("Solve0_", i)]] <- lapply(test$Sequence, FUN = solveRecurrentNobias, depth = i, takeLast = FALSE)
    test[[n("Solve1_", i)]] <- lapply(test$Sequence, FUN = solveRecurrentBias, depth = i, takeLast = FALSE)
    test[[n("SolveL0_", i)]] <- lapply(test$Sequence, FUN = solveRecurrentNobias, depth = i, takeLast = TRUE)
    test[[n("SolveL1_", i)]] <- lapply(test$Sequence, FUN = solveRecurrentBias, depth = i, takeLast = TRUE)

    test[[n("Check0_", i)]] <- mapply(isRecurrent, test$Sequence, i, test[[n("Solve0_", i)]], FALSE)
    test[[n("Check1_", i)]] <- mapply(isRecurrent, test$Sequence, i, test[[n("Solve1_", i)]], TRUE)
    test[[n("CheckL0_", i)]] <- mapply(isRecurrent, test$Sequence, i, test[[n("SolveL0_", i)]], FALSE)
    test[[n("CheckL1_", i)]] <- mapply(isRecurrent, test$Sequence, i, test[[n("SolveL1_", i)]], TRUE)
    test[[n("Check_", i)]] <- test[[n("Check0_", i)]] | test[[n("Check1_", i)]] | test[[n("CheckL0_", i)]] | test[[n("CheckL1_", i)]]

    found <- sum(test[[n("Check_", i)]] == TRUE)
    print(paste("---- found ", found, " ----- ", Sys.time()))
}
sink()
#unlink("../output/iterations.txt")


coalesce1a <- function(...) {
    ans <- ..1
    for (elt in list(...)[-1]) {
        i <- which(is.na(ans))
        ans[i] <- elt[i]
    }
    ans
}


for (i in 1:maxDepth) {
    print(paste("Predicting ", i))
    test[[n("Last0_", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$BigSequence, i, test[[n("Solve0_", i)]], test[[n("Check0_", i)]], FALSE)
    test[[n("Last1_", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$BigSequence, i, test[[n("Solve1_", i)]], test[[n("Check1_", i)]], TRUE)
    test[[n("LastL0_", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$BigSequence, i, test[[n("SolveL0_", i)]], test[[n("CheckL0_", i)]], FALSE)
    test[[n("LastL1_", i)]] <- mapply(predictNext, SIMPLIFY = TRUE, test$BigSequence, i, test[[n("SolveL1_", i)]], test[[n("CheckL1_", i)]], TRUE)
    test[[n("Last_", i)]] <- mapply(coalesce1a, test[[n("Last0_", i)]], test[[n("Last1_", i)]], test[[n("LastL1_", i)]], test[[n("LastL0_", i)]])
}

test$Last <- NA
for (cn in sapply(c(1:maxDepth), function(x) paste("Last_", x, sep = ""))) {
    test$Last <- coalesce1a(test$Last, test[[cn]])
}


test$Last[is.na(test$Last)] <- 0
write.csv(test[, c("Id", "Last")], "../output/recurrent-exp6.csv", row.names = FALSE)



