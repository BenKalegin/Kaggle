# lm fit
formulaString <- "y~"
for (i in 1:depth) {
    df[[paste0("x", i)]] <- x[i:(length(x) - depth + i - 1)]
    formulaString <- paste0(formulaString, "+x", i)
}
formulaString <- sub("~\\+", "~", formulaString)
fit <- lm(formula(formulaString), df)