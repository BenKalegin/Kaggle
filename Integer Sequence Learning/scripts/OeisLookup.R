prepareOesi <- function() {
    #oeis.org/classic/stripped.gz
    oesi <- readLines("../download/stripped.csv")
    oesi <- oesi[ - (1:4)]
    oesi <- sapply(strsplit(oesi, split = ","), FUN = function(x) unlist(x))
    trim.trailing <- function(x) sub("\\s+$", "", x)
    oesiId <- sapply(oesi, FUN = function(x) trim.trailing(head(x, 1)))
    oesiSequence <- sapply(oesi, FUN = function(x) paste(unlist(x[-1]), collapse = ","))
    oesiId <- trim.trailing(oesiId)
    df <- data.frame(Id = oesiId, Sequence = oesiSequence)
    write.csv(df, "../download/oesi.csv", row.names = FALSE)
}

library(data.table)
matchTestSet <- function() {
    lookup <- data.table(read.csv("../download/oesi.csv", stringsAsFactors = FALSE))
    lookup <- lookup[Sequence != ""]
    setkey(lookup, Sequence)
    test <- data.table(read.csv("../input/test.csv", stringsAsFactors = FALSE))
    test$Last <- sapply(test$Sequence, FUN = searchLookup, lookup = lookup)
    test[, Sequence := NULL]
    write.csv(test, "../download/lookup.csv", row.names = FALSE)
}

startsWith = function(x, prefix) {
    if (!is.character(x) || !is.character(prefix))
        stop("non-character object(s)")
    suppressWarnings(substr(x, 1L, nchar(prefix)) == prefix)
}

lookupLast <- function(id, lookup) {

    lookup$
}

searchLookup <- function(value, lookup) {
    low = 1
    high = nrow(lookup)
    mid <- 1
    while (low <= high) {
          # invariants:value > A[i] for all i < low, value < A[i] for all i > high
        mid = ceiling((low + high) / 2)
        s <- lookup$Sequence[mid]
        if (s > value)
            high = mid - 1
        else if (s < value)
            low = mid + 1
        else {
            break;
        }
    }
    if ((mid < nrow(lookup) - 1) & (startsWith(lookup$Sequence[mid + 1], value)))
        mid <- mid + 1
    diff <- substring(lookup$Sequence[mid], nchar(value)+2)
    diffSeq <- unlist(strsplit(diff, split = ","))
    diffSeq[1]
}


findOrder2 <- function() {
    install.packages("rjson")
    test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
    test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = unlidt)
    empty <- lapply(test$Sequence, length)

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
}
