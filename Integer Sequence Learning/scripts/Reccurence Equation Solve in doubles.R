library(gmp)

solveRecurrent <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    if (length(x) - 2 * depth < 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        A <- matrix(NA, nrow = depth + 1, ncol = depth + 1, byrow = TRUE)
        b <- matrix(NA, nrow = depth + 1, ncol = 1)
        offset <- 1
        if (takeLast == TRUE) {
            offset <- length(x) - 2 * depth
            while ((offset > 1) & (abs(x[offset]) > 1e4)) offset <- offset - 1
        }
        for (r in 1:(depth + 1)) {
            for (c in 1:depth) {
                if (offset + (r - 1) + (c - 1) > length(x)) {
                    xxx <- 1;
                    stop();
                }
                A[r, c] <- x[offset + (r - 1) + (c - 1)]
            }
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
        for(i in 1:depth)
        {
            df[[paste0("x",i)]] <- x[length(x)-depth+i]
        }
        df <- as.data.frame(df)
        prediction <- predict(fit, df)

        prediction <- round(prediction)
    }
    prediction
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
    test[[n("Solve", i)]] <- sapply(test$Sequence, FUN = solveRecurrent, depth = i, takeLast = TRUE)

    test[[n("Check", i)]] <- mapply(isRecurrent, test$BigSequence, i, test[[n("Solve", i)]])
    found <- sum(test[[n("Check", i)]] == TRUE)
    print(paste("---- found ", found, " ----- ", Sys.time()))

    if (found == 0) {
        maxDepth <- i
        break;
    }
}
sink()
#unlink("../output/iterations.txt")


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


test$Last[is.na(test$Last)] <- 0 
write.csv(test[, c("Id", "Last")], "../output/recurrent-exp4.csv", row.names = FALSE)

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
exp2 <- read.csv("../output/recurrent-exp2.csv")
exp4 <- read.csv("../output/recurrent-exp4.csv")
exp42 <- exp2[exp2$Last != exp4$Last & exp4$Last == 0,]
write.csv(exp42, "../output/exp4-exp2.csv", row.names = TRUE)
"C:\Users\Ben\Downloads\Solve recurrent exp4.csv"

test$table <- sapply(test$Sequence, FUN = table)
test$repeats <- mapply(function(table, seq) length(table) < sqrt(length(seq)), test$table, test$Sequence)
test$periodic <- mapply(function(table, seq) {
    result <- length(table) < sqrt(length(seq))
    if (result) {
        result <- all(abs(diff(unlist(table))) < 2)
    }
    result
}, test$table, test$Sequence)
#sum(test$periodic == TRUE) 351