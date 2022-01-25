new_jive = function (y, X, Z) 
{
    n <- length(y)
    k <- dim(X)[2]
    Xj <- matrix(0, n, k)
    iZZ <- solve(t(Z) %*% Z)
    Ga.hat <- (iZZ) %*% (t(Z) %*% X)
    h <- vector("numeric", n)
    for (i in 1:n) {
        h[i] <- t(Z[i, ]) %*% iZZ %*% Z[i, ]
        Xj[i, ] <- (t(Z[i, ]) %*% Ga.hat - h[i] * X[i, ])/(1 - 
            h[i])
    }
    XXj <- t(Xj) %*% X
    iXXj <- solve(XXj)
    B <- (iXXj) %*% (t(Xj) %*% y)
    return(B = list(beta = B, x_hat = Xj))
}