test_functions <- function(numSteps, FC, fArr_sym_den, mfpt_sym_den){
    FC <- FC[upper.tri(FC)]

    div <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- mfpt_sym_den[steps, , ] / fArr_sym_den[steps, , ] 
        mat <- mat[upper.tri(mat)]
        div[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    
    
    subt <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- abs(fArr_sym_den[steps, , ] - mfpt_sym_den[steps, , ])
        mat <- mat[upper.tri(mat)]
        subt[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    
    expo <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- fArr_sym_den[steps, , ] ^ mfpt_sym_den[steps, , ] 
        mat <- mat[upper.tri(mat)]
        expo[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    

    
    stepVec <- seq(1, numSteps)
    
    
    cor_matrix = as.data.frame(cbind(stepVec, abs(div),
                                     abs(subt), abs(expo)))
    colnames(cor_matrix) <- c("stepVec","div", "subt", "expo")

    
    
    library(ggplot2)
    ggplot(data = cor_matrix, aes(x= stepVec, y = value, color = variable)) + 
        geom_point(aes(y = div, col = "divide")) +  
        geom_line(aes(y = div, col = "divide")) + 
        geom_point(aes(y = subt, col = "subtract")) + 
        geom_line(aes(y = subt, col = "subtract")) + 
        geom_point(aes(y = expo, col = "exponential")) + 
        geom_line(aes(y = expo, col = "exponential")) + 
        labs(x="Step", y="Correlation")
    
    ggsave(paste(getwd(),"/limited_time_mfpt","/results","/corr_funct.png", sep = '/'))
    return(cor_matrix)
}
