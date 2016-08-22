# svd https://www.ecse.rpi.edu/~qji/CV/svd_review.pdf
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

k = 15
x[k] * x[k - 2] / (x[k - 1] ^ 2)