submission17 <- read.csv("../submissions/17 (17772)/submission.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character"))

test$Repeated <- sapply(test$Sequence, function(x) {
    x <- unlist(x)
    value <- (length(x) > 1) & (x[length(x)] == x[length(x) - 1])
    value
})

sum(test$Repeated) #5524
r <- which(test$Repeated == TRUE)
submission17$Last[r] <- sapply(test$Sequence[r], tail, 1)
sum(submission17$Last != 0)
write.csv(submission17, "../submissions/17 (17772)/submissionWithRepeats.csv", row.names = FALSE)

same as 17 plus repeats