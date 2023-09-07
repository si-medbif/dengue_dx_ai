library(tidyverse)
library(readr)

#Load data
dat <- as.data.frame(read_csv("RawData/data_continuous.csv"))

#Specify and code factor variables
fcol <- c("d01_sex", "d04_abdominal"          
         ,"d04_diarrhea"   ,         "d04_itching"            
         , "d04_limbus"   ,           "d04_liver" 
         ,"d04_maculo"    ,          "d04_rash"   
         ,"d04_uri" )

for(col in fcol){
  dat[[col]] <- factor(dat[[col]])
}

#Fit a logistic regression model
form <- formula(paste0("dhf_dx ~ ",paste(names(dat)[3:38],collapse = "+")))
mall <- glm(form , data = dat, family = "binomial")

#Extract odds ratios
multivar <- data.frame(exp(cbind(Odds_Ratio = coef(mall), confint(mall))))

#Extract a p-value for each independent variable
sm <- summary(mall)
pval <- data.frame(sm$coefficients)
pval$variable <- row.names(pval)
multivar$variable <- rownames(multivar)
multivar <-merge(multivar,pval, by = "variable", all = T)

#Save results
write.csv(multivar,"Results/Multivar_Odds_ratio_continuous.csv", row.names = F)
