makedelta <- function(data){
    rows = dim(data)[2]
    cols = dim(data)[3]
    steps = dim(data)[1]
    for(k in 1:steps){
        for(i in 1:cols){
            rt <- data[k,i,i]
            data[k,,i] <- data[k,,i] - rt
        }
    }
    return(data)
}