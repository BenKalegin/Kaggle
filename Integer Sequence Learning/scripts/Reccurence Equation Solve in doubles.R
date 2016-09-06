library(gmp)

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
        result <- match > notmatch
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

makeRecurrentMatrix <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    A <- matrix(NA, nrow = depth, ncol = depth, byrow = TRUE)
    b <- matrix(NA, nrow = depth, ncol = 1)

    if (length(x) - 2 * depth > 1) {
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

    }
    list(A, b)
}


solveSvd <- function(x, depth, takeLast = TRUE) {
    Ab <- makeRecurrentMatrix(x, i, takeLast)
    result <- matrix(NA, nrow = depth, ncol = 1)
    if (!anyNA(Ab)) {
        A <- Ab[[1]]
        b <- Ab[[2]]
        a.svd <- svd(A)
        if (!(0 %in% a.svd$d)) {
            ds <- diag(1 / a.svd$d[1:depth])
            u <- a.svd$u
            v <- a.svd$v
            us <- as.matrix(u[, 1:depth])
            vs <- as.matrix(v[, 1:depth])
            bs <- unlist(b)
            xs <- vs %*% ds %*% (t(us) %*% bs)
            round(xs, 1)
            result <- unlist(a.ginv <- vs %*% ds %*% (t(us) %*% bs))
        }
    }
    result
}


#system.time(rec4 <- lapply(sequences, isRecurrent, depth = 4)) # 1.7 on 10000, then 3.07 when value == x[] changed to abs(value - x) < 0.01, but records found 6120 instead of 227

test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
test$Sequence[210] <- gsub("38280596832649216", "38280596832649217", test$Sequence[210])
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

maxDepth <- limitDepth
for (i in 1:limitDepth) {
    print(paste("---- iteration ", i, " ----- ", Sys.time()))
    for (j in 1:nrow(test)) {
        solve0 <- solveSvd[test$Sequence[j]], depth = i, takeLast = FALSE)
        if (isRecurrent(test$Sequence[j], i, test[j, n("Solve0_", i)], FALSE))
        solve1 <- solveSvd[test$Sequence[j]], depth = i, takeLast = TRUE)
        
        test[j, n("Check1_", i)] <- isRecurrent(test$Sequence[j], i, test[j, n("Solve1_", i)], TRUE)
        test[[n("Check_", i)]] <- test[[n("Check0_", i)]] | test[[n("Check1_", i)]]
    }
    found <- sum(test[[n("Check_", i)]] == TRUE)
    print(paste("---- found ", found, " ----- ", Sys.time()))
}


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
    test[[n("Last_", i)]] <- mapply(coalesce1a, test[[n("Last0_", i)]], test[[n("Last1_", i)]])
}
test$Last <- NA
for (cn in sapply(c(1:maxDepth), function(x) paste("Last_", x, sep = ""))) {
    test$Last <- coalesce1a(test$Last, test[[cn]])
}



test$Repeated <- sapply(test$Sequence, function(x) {
    x <- unlist(x)
    value <- (length(x) > 1) & (x[length(x)] == x[length(x) - 1])
    value
})

sum(test$Repeated) #5524
r <- which(test$Repeated == TRUE & is.na(test$Last))
test$Last[r] <- sapply(test$Sequence[r], tail, 1)

test$Last[is.na(test$Last)] <- 0
write.csv(test[, c("Id", "Last")], "../output/recurrent-exp7.csv", row.names = FALSE)

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

# experiment 3: same as 2 but solve get numbers from the end
#[1] "---- iteration  1  -----  2016-08-06 18:46:15"
#[1] "---- found  1267  -----  2016-08-06 18:47:29"
#[1] "---- iteration  2  -----  2016-08-06 18:47:29"
#[1] "---- found  2484  -----  2016-08-06 18:49:11"

# experiment 4: same as 3 but move index down until less 1e6
# [1] "---- iteration  1  -----  2016-08-07 01:57:02"
# [1] "---- found  1564  -----  2016-08-07 01:57:51"
# [1] "---- iteration  2  -----  2016-08-07 01:57:51"
# [1] "---- found  2863  -----  2016-08-07 01:59:00"

# experiment 4-1: same as 3 but move index down until less 1e5
#[1] "---- iteration  1  -----  2016-08-07 03:07:27"
#[1] "---- found  1564  -----  2016-08-07 03:08:14"
#[1] "---- iteration  2  -----  2016-08-07 03:08:14"
#[1] "---- found  2900  -----  2016-08-07 03:09:25"

# experiment 4-2: same as 3 but move index down until less 1e4
#[1] "---- iteration  1  -----  2016-08-07 13:41:28"
#[1] "---- found  1563  -----  2016-08-07 13:42:16"
#[1] "---- iteration  2  -----  2016-08-07 13:42:16"
#[1] "---- found  2901  -----  2016-08-07 13:43:29"

#submission 4-2: score 0.15105

#submission 5-1: score 0.15371
#submission 5-2(15): score 0.15453 ???? ?? 5-1 ??? 

#exp6: isRecurrent 1e-3 tolerance



