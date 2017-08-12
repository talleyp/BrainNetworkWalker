cor_matrices <- function(numSteps, FC, fArr_sym_den, mfpt_sym_den, n_SCDEN, z_zym, delta){
    FC <- FC[upper.tri(FC)]
    
    cor_FC_FA <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- fArr_sym_den[steps, , ]
        mat <- mat[upper.tri(mat)]
        cor_FC_FA[steps] <- abs(cor(FC, mat))
    }
    
    cor_FC_mfpt <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- mfpt_sym_den[steps, , ]
        mat <- mat[upper.tri(mat)]
        cor_FC_mfpt[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }

    cor_n_scden <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- n_SCDEN[steps, , ]
        mat <- mat[upper.tri(mat)]
        cor_n_scden[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    
    cor_z <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- z_zym[steps, , ]
        mat <- mat[upper.tri(mat)]
        cor_z[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    
    cor_delta <- rep(NULL, numSteps)
    for(steps in 1:numSteps){
        mat <- delta[steps, , ]
        mat <- mat[upper.tri(mat)]
        cor_delta[steps] <- abs(cor(FC, mat, use = "pairwise.complete.obs"))
    }
    
    per_na <- rep(NULL,numSteps)
    for(steps in 1:numSteps){
        total <- length(c(mfpt_sym_den[steps, , ]))
        per_na[steps] <- sum(is.na(c(mfpt_sym_den[steps, , ]))) / total
    }
    
    stepVec <- seq(1, numSteps)

    
    cor_matrix = as.data.frame(cbind(stepVec, abs(cor_FC_FA), abs(cor_FC_mfpt), abs(cor_n_scden), 
                                     abs(cor_z), abs(cor_delta), per_na))
    colnames(cor_matrix) <- c("stepVec","cor_FC_FA","cor_FC_mfpt", "cor_n_scden","cor_z","cor_delta","per_na")
    
#     maxMFPT <- max(abs(cor_matrix$cor_FC_mfpt), na.rm=TRUE)
#     locMFPT <- which(abs(cor_matrix$cor_FC_mfpt) == max(abs(cor_matrix$cor_FC_mfpt), na.rm=TRUE), arr.ind = TRUE)
#     
#     maxFA <- max(abs(cor_matrix$cor_FC_FA), na.rm=TRUE)
#     locFA <- which(abs(cor_matrix$cor_FC_FA) == max(abs(cor_matrix$cor_FC_FA), na.rm=TRUE), arr.ind = TRUE)
#     
#     maxN <- max(abs(cor_matrix$cor_n_scden), na.rm=TRUE)
#     locN <- which(abs(cor_matrix$cor_n_scden) == max(abs(cor_matrix$cor_n_scden), na.rm=TRUE), arr.ind = TRUE)    

    
    library(ggplot2)
    ggplot(data = cor_matrix, aes(x= stepVec, y = value, color = variable)) + 
        geom_point(aes(y = cor_FC_FA, col = "fracA")) + 
        geom_line(aes(y = cor_FC_FA, col = "fracA")) + 
        geom_point(aes(y = cor_FC_mfpt, col = "mfpt")) +  
        geom_line(aes(y = cor_FC_mfpt, col = "mfpt")) + 
        geom_point(aes(y = cor_n_scden, col = "nArr")) + 
        geom_line(aes(y = cor_n_scden, col = "nArr")) + 
        geom_point(aes(y = cor_z, col = "Z")) + 
        geom_line(aes(y = cor_z, col = "Z")) + 
        geom_point(aes(y = cor_delta, col = "Delta")) + 
        geom_line(aes(y = cor_delta, col = "Delta")) + 
        geom_point(aes(y = per_na, col = "% NA")) + 
        geom_line(aes(y = per_na, col = "% NA")) + 
#         annotate("text", x = locFA, y = maxFA, label = paste("Max FA = ", maxFA)) +
#         annotate("text", x = locMFPT, y = maxMFPT, label = paste("Max FA = ", maxMFPT)) +
#         annotate("text", x = locMFPT, y = maxN, label = paste("Max FA = ", maxN)) +
        labs(x="Step", y="Correlation")
    
    ggsave(paste(getwd(),"/limited_time_mfpt","/results","/corr_to_fc.png", sep = '/'))
    return(cor_matrix)
}
