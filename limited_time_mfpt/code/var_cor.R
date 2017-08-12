
cor_mfpt_fa <- rep(NULL, numSteps)
for(steps in 1:numSteps){
    cor_mfpt_fa[steps] <- cor(c(fracArr_SCden[steps, , ]), c(mfpt_SCden[steps, ,]),use = "pairwise.complete.obs")
}

cor_mfpt_narr <- rep(NULL, numSteps)
for(steps in 1:numSteps){
    cor_mfpt_narr[steps] <- cor(c(nArr_SCden[steps, , ]), c(mfpt_SCden[steps, ,]),use = "pairwise.complete.obs")
}

cor_narr_fa <- rep(NULL, numSteps)
for(steps in 1:numSteps){
    cor_narr_fa[steps] <- cor(c(fracArr_SCden[steps, , ]), c(nArr_SCden[steps, ,]),use = "pairwise.complete.obs")
}

stepVec <- seq(1, numSteps)

cor_matrix = as.data.frame(cbind(stepVec, cor_mfpt_fa, cor_mfpt_narr, cor_narr_fa))
colnames(cor_matrix) <- c("stepVec","cor_mfpt_fa","cor_mfpt_narr", "cor_narr_fa")


library(ggplot2)
ggplot(data = cor_matrix, aes(x= stepVec, y = value, color = variable)) + 
    geom_point(aes(y = cor_mfpt_fa, col = "fracA-mfpt")) + 
    geom_line(aes(y = cor_mfpt_fa, col = "fracA-mfpt")) + 
    geom_point(aes(y = cor_mfpt_narr, col = "mfpt-narr")) +  
    geom_line(aes(y = cor_mfpt_narr, col = "mfpt-narr")) + 
    geom_point(aes(y = cor_narr_fa, col = "nArr-fracA")) + 
    geom_line(aes(y = cor_narr_fa, col = "nArr-fracA")) + 
    labs(x="Step", y="Correlation")

ggsave(paste(getwd(),"/limited_time_mfpt","/results","/corr_vars.png", sep = '/'))

