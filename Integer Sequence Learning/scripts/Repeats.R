#exact repeat at the end
test$Repeated <- sapply(test$Sequence, function(x) {
    x <- unlist(x)
    value <- (length(x) > 1) & (x[length(x)] == x[length(x) - 1])
    value
})
sum(test$Repeated == TRUE, na.rm = TRUE) #5524
r <- which(test$Repeated == TRUE & submission17$Last == 0)
submission <- read.csv("../submissions/21 (18358)/submission.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character"))

submission$Last[r] <- sapply(test$Sequence[r], tail, 1)
r <- which(submission$Last == 0)
submission$Last[r] <- 1
write.csv(submission, "../submissions/21 (18358)/submissionWithRepeats.csv", row.names = FALSE)

#repeatable groups
test$table <- sapply(test$Sequence, FUN = table)
test$repeats <- mapply(function(table, seq) length(table) < sqrt(length(seq)), test$table, test$Sequence)
test$periodic <- mapply(function(table, seq) {
    result <- length(table) < sqrt(length(seq))
    if (result) {
        result <- all(abs(diff(unlist(table))) < 2)
    }
    result
}, test$table, test$Sequence)

test$diff1 <- sapply(1:nrow(test), function(i) {
    length(which(diff(x, 1) == 0)) > length(x) / 2
})
sum(test$diff1 == TRUE, na.rm = TRUE) #0

test$diff2 <- sapply(1:nrow(test), function(i) {
    x <- unlist(test$Sequence[i])
    length(which(diff(x, 2) == 0)) > length(x) / 2
})
sum(test$diff2 == TRUE, na.rm = TRUE) #3207
head(test[test$diff2 == TRUE,])

test$diff3 <- sapply(1:nrow(test), function(i) {
    x <- unlist(test$Sequence[i])
    length(which(diff(x, 3) == 0)) > length(x) / 2
})
sum(test$diff3 == TRUE & test$diff2 == FALSE, na.rm = TRUE) #521
head(test[test$diff3 == TRUE,])
diff(unlist(test$Sequence[178]))