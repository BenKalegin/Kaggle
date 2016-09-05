# svd https://www.ecse.rpi.edu/~qji/CV/svd_review.pdf

makeRecurrentMatrix <- function(x, depth, takeLast = TRUE) {
    x <- unlist(x)
    A <- matrix(NA, nrow = depth, ncol = depth, byrow = TRUE)
    b <- matrix(NA, nrow = depth, ncol = 1)

    if (length(x) - 2 * depth > 1) {
        offset <- 1
        if (takeLast == TRUE) {
            offset <- length(x) - 2 * depth
            while ((offset > 1) & (abs(x[offset]) > 1e4))
                offset <- offset - 1
            }
        for (r in 1:depth) {
            A[r, 1:depth] <- x[(offset + (r - 1)):(offset + (r - 1) + depth - 1)]
            b[r] <- x[offset + depth + (r - 1)]
        }

    }
    list(A, b)
}


#test <- read.csv("../input/test.csv", stringsAsFactors = FALSE, nrow = 10000000)
#test$Sequence <- sapply(strsplit(test$Sequence, split = ","), FUN = as.numeric)

depth = 4
Ab <- makeRecurrentMatrix(test$Sequence[2], depth, F)
A <- Ab[[1]]
b <- Ab[[2]]
a.svd <- svd(A)
ds <- diag(1 / a.svd$d[1:depth])
u <- a.svd$u
v <- a.svd$v
us <- as.matrix(u[, 1:depth])
vs <- as.matrix(v[, 1:depth])
bs <- unlist(b)
xs <- vs %*% ds %*% (t(us) %*% bs)
round(xs, 1)
A %*% xs - bs
(a.ginv <- vs %*% ds %*% (t(us) %*% bs))
a.ginv %*% A
#(0, 10, 0, -15, 0, 7, 0, -1)

