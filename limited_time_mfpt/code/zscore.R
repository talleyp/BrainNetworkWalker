zscore <- function(data){
    rows = dim(data)[2]
    cols = dim(data)[3]
    steps = dim(data)[1]
    for(k in 1:steps){
        for(i in 1:cols){
            m <- mean(data[k,i,], na.rm = T)
            s <- sd(data[k,i,], na.rm = T)
            
            for(j in 1:rows){
                data[k,i,j] = (data[k,i,j] - m) / s
            }
        }
    }
    return(data)
}