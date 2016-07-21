# load the data 
train <- read.csv("../input/train.csv", stringsAsFactors = FALSE)
sequences <- strsplit(train$Sequence, split = ",")

# What is the most frequent final element of each sequence in the training set?

finalElements <- sapply(sequences, tail, n = 1)
frequencies <- table(finalElements)
index <- which.max(frequencies)
mostFrequentFinalElement <- names(frequencies)[index]
mostFrequentFinalElement