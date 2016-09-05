exp62 <- read.csv("../output/recurrent62.csv")
exp5 <- read.csv("../output/recurrent-exp5.csv")
exp562 <- exp62[exp62$Last != exp5$Last & exp5$Last == 0,]
write.csv(exp562, "../output/exp5-62.csv", row.names = TRUE)

exp5 <- read.csv("../submissions/17 (17772)/submission.csv")
exlm <- read.csv("../output/linearPrevious10WithModeFallback.csv")
exp5lm <- exlm[exlm$Last != exp5$Last & exp5$Last == 0,]


write.csv(exp5lm, "../output/exp17-lm.csv", row.names = TRUE)

submission <- read.csv("../submissions/22 (18735)/submission.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character"))
library(data.table)
verified <- data.table(read.csv("../download/lookup.csv", stringsAsFactors = FALSE, colClasses = c("integer", "character")))
sum(submission$Last == verified$Last, na.rm = TRUE) #21638
sum(exlm$Last == verified$Last, na.rm = TRUE) #18811

sum(exlm$Last == verified$Last & submission$Last != verified$Last, na.rm = TRUE) #5728
missed <- which((exlm$Last == verified$Last) & (submission$Last != verified$Last))
missedInfo <- test[missed,]
missedInfo$lm <- exlm$Last[missed]
missedInfo$verified <- verified$Last[missed]
missedInfo$solved <- submission$Last[missed]

x <- missedInfo[2,]$Sequence
solveRecurrentOld(x, 5)
missedInfo[5:10,]


