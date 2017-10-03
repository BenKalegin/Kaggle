library(gmp)
test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 1000000)
#test$Sequence[210] <- gsub("38280596832649216", "38280596832649217", test$Sequence[210])
strSequence <- strsplit(test$Sequence, split = ",")
test$BigSequence <- sapply(strSequence, FUN = as.bigz)
test$Sequence <- sapply(strSequence, FUN = as.numeric)
test$Solved <- NA
test$Formula <- NA

for (order in 1:10) {
    print(paste("---- Order ", order, " ----- ", Sys.time()))

    formulaString <- "y ~ x"
    for (i in 2:order) {
        formulaString <- paste0(formulaString, "+I(x^", i, ")")
    }

    for (j in 1:(nrow(test))) {
        if (j %% 1000 == 0)
            print(j)

        if (is.na(test$Solved[j])) {
            y <- data.frame(y = unlist(test$Sequence[j]), x = 1:length(test$Sequence[j]))
            fit <- lm(formula(formulaString), y)
            if (sum(abs(fit$residuals)) < nrow(test) / 2) {
                test$Solved[j] = "Polyfit"
                test$Formula[j] = paste(as.character(fit$coefficients), collapse = ", ")
                #prediction <- predict(fit)
            }
        }
    }
}
