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
