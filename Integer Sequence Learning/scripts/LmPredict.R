lmRecurrent <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    if (length(x) - depth < 1) {
        result <- matrix(NA, nrow = depth, ncol = 1);
    } else {
        if (takeLast) {
            df <- data.frame(y = tail(x, - depth))
        }
        formulaString <- "y~"
        for (i in 1:depth) {
            df[[paste0("x", i)]] <- x[i:(length(x) - depth + i - 1)]
            formulaString <- paste0(formulaString, "+x", i)
        }
        formulaString <- sub("~\\+", "~", formulaString)

        fit <- lm(formula(formulaString), df)
        maxResidual <- max(abs(fit$residuals))

        df <- list()
        for (i in 1:depth) {
            df[[paste0("x", i)]] <- x[length(x) - depth + i]
        }
        df <- as.data.frame(df)
        prediction <- predict(fit, df)

        prediction <- round(prediction)
    }
    prediction
}

