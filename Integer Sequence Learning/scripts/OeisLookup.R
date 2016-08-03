install.packages("rjson")
test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.character)

#list <- list.files("../download", pattern = "*.json")
library(rjson)

test$isOrder2 <- FALSE;
for (i in 1:1000) {
    print(i)
    json_file <- paste(sep = "", "../download/", test$Id[i], ".json")
    # 357.json
    # 484.json
    # 1741.json
    
    if (file.exists(json_file)) {
        json <- fromJSON(file = json_file)
        if (json$count > 0 & !is.null(json$results$link)) {
            for (j in 1:length(json$results$link)) {
                print("links")
                isRec <- grep("Rec#order_02", json$results$link[j])
                if (isRec) {
                    print("FOUND")
                    test$isOrder2[i] <- TRUE
                }
            }
        }
    }
}
