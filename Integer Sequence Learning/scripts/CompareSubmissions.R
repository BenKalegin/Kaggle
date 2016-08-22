exp62 <- read.csv("../output/recurrent62.csv")
exp5 <- read.csv("../output/recurrent-exp5.csv")
exp562 <- exp62[exp62$Last != exp5$Last & exp5$Last == 0,]
write.csv(exp562, "../output/exp5-62.csv", row.names = TRUE)

exp5 <- read.csv("../output/recurrent-exp5.csv")
exlm <- read.csv("../output/linearPrevious10WithModeFallback.csv")
exp5lm <- exlm[exlm$Last != exp5$Last & exp5$Last == 0,]
write.csv(exp5lm, "../output/exp5-lm.csv", row.names = TRUE)

submission17 <- read.csv("../submissions/17 (17772)/submission.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character"))
verified <- data.table(read.csv("../download/lookup.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character")))
setkey(verified, Id)
i <- 1
s <- sapply(1: nrow(submission17), function(i) {
    if (i %% 1000 == 0) {
        print(paste("Done", i, "sequences"))
    }
    verified[J(submission17$Id[i])]$Last
})
submission17$Verified <- s
sum(submission17$Verified == submission17$Last) #19918