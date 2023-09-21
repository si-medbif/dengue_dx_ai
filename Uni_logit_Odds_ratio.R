#R script for a univariate logistic regression analysis for deriving odds ratios for independent variables to be used in DHF prognosis models.

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


write.csv(multivar,"Results/Multivar_Odds_ratio_continuous.csv", row.names = F)


library(tidyverse)
library(readr)

#Load data
dat <- as.data.frame(read_csv("RawData/data_continuous_forOddRatio.csv"))

result <- data.frame()

# Correct variables with skewed distribution
dat$d05_lft_ast <- log10(dat$d05_lft_ast)
dat$d05_lft_alt <- log10(dat$d05_lft_alt)
dat$d05_sgot_platelet_ratio <- log10(dat$d05_sgot_platelet_ratio)
dat$d05_platelet <- log10(dat$d05_platelet * 10^3)

cols <- names(dat)[2:41]

for(col in cols){
form1 <- formula(paste(col,"~ dhf_dx"))
form2 <- formula(paste("dhf_dx ~",col))

#Find median and IQR for each variable
ares <- as.matrix(aggregate(form1, data =dat, FUN = quantile, probs = c(0.5,0.25,0.75), na.rm = T))
ares <- round(ares,3)
ares_notDHF <- paste0(ares[1,][2], " (", ares[1,3],", ",ares[1,4],")")
ares_DHF <- paste0(ares[2,][2], " (", ares[2,3],", ",ares[2,4],")")

#Fit a logistic regression model         
m <- glm(form2,data = dat[which(dat[[col]] != -Inf & !is.na(dat[[col]]) & !is.nan(dat[[col]])),], family =  "binomial")
sm <- summary(m)

#Extract odds ratio and p-value
p.val <- sm$coefficients[2,4]
ors <- format(round(as.matrix(exp(cbind(coef(m),confint(m)))),3), nsmall = 3)
ors <- paste0(ors[2,1], " (",ors[2,2],", ",ors[2,3],")")

t_res <- data.frame(Variale = col, NonDHF = ares_notDHF, DHF = ares_DHF , OddsRatio = ors, Pval = p.val)
result <- rbind(result,t_res)
}

#Format p-value
result$FormatP.val <- format(round(result$Pval,3), nsmall = 3)
result$FormatP.val[result$Pval < 0.001] <- "<0.001" 

#Save results
write.csv(result,file = "Results/Univar_Odds_ratio_continuous.xlsx")
