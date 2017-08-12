makesym <- function(mat, N){
    for(i in 1:N){
        sym <- (mat[i, , ] + t(mat[i, ,])) / 2
        mat[i, ,] <- sym
        diag(mat[i,,]) <- 1
    }
    return(mat)
}