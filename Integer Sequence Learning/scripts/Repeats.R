#exact repeat at the end
test$Repeated <- sapply(test$Sequence, function(x) {
    x <- unlist(x)
    value <- (length(x) > 1) & (x[length(x)] == x[length(x) - 1])
    value
})

sum(test$Repeated) #5524
r <- which(test$Repeated == TRUE & submission17$Last == 0)
submission17$Last[r] <- sapply(test$Sequence[r], tail, 1)
r <- which(submission17$Last == 0)
submission17$Last[r] <- 1
write.csv(submission17, "../submissions/17 (17772)/submissionWithRepeats.csv", row.names = FALSE)

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