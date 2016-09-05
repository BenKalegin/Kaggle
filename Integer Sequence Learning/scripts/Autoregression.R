test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.numeric)
test$Predict <- sapply(1:nrow(test), function(i) {
    if ((i %% 10) == 0)
        print(i)
    rm <- suppressWarnings(try(ar.ols(unlist(test$Sequence[i]))));
    #print(class(rm))
    result <- NA
    if ("try-error" != class(rm))
        result <- predict(rm)$pred[1]
    result
})

sum(test$Predict == verified$Last & submission$Last != verified$Last, na.rm = TRUE) #5728
arfound <- which(test$Predict == verified$Last & submission$Last != verified$Last)
test[arfound,]
x <- unlist(test$Sequence[1361])
str(ar.ols(x))