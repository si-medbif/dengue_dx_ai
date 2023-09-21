library(readr)
library(lubridate)
library(tidyverse)
library(stringr)
library(ggpubr)
library(rstatix)

#Load data
dat <- read_csv("RawData/REDCap_raw_data.csv")
dat1 <- read_csv("RawData/REDCap_label_data.csv")

theme_set(theme_bw())

# Format Dx and Day0
DX_Day0 <- dat1[!is.na(dat1$`Final Dx`),c("SNO:","Final Dx","Final D0")]
DX_Day0$simple_Dx <- "DHF"
DX_Day0$simple_Dx[which(DX_Day0$`Final Dx` == "Not dengue")] <- "OFI"
DX_Day0$simple_Dx[which(DX_Day0$`Final Dx` == "DF")] <- "DF"
DX_Day0$final_day0 <- ymd(DX_Day0$`Final D0`)
DX_Day0$`Final Dx` <- DX_Day0$`Final D0` <- NULL

dat <- merge(dat,DX_Day0, by.x = "sno",by.y="SNO:")

# Format Site
dat$site <- factor((substr(dat$sno,3,3)),levels = 1:2, labels = c("KK","SK"))

# Format to_day0
dat$date <- ymd(dat$date_d04)
dat$to_day0 <- dat$date - dat$final_day0

# Format study day
dat$study_day <- as.numeric(str_match(dat$redcap_event_name, ".?study_day_([0-9])_arm_1")[,2])

# Format fever day
fever <- dat[!is.na(dat$fever_date), c("sno","fever_date")]
fever$f_date <- ymd(fever$fever_date)
dat <- merge(dat,fever, by = "sno", all.x = T)
dat$fever_day <- dat$date - dat$f_date

# Format WBC
dat$log.wbc <- log10(1+dat$wbc_d05*10^3)

# Format LYMP
dat$log.lymp <- log10(1+(dat$lymp * dat$wbc_d05*10^3)/100)

# Show only DF versus DHF
dat <- dat[which(dat$simple_Dx != "OFI"),]
dat$simple_Dx <- factor(dat$simple_Dx, levels = c("DHF","DF"))

# Select fever day 3 to 6
fdat2 <- dat[which(dat$simple_Dx != "OFI" & dat$fever_day %in% 3:6 & dat$log.pmn > 0),]
fdat2$fever_day_F <- paste0("Fever day ",fdat2$fever_day)

# Compare data at each fever day
stat.test <- fdat2 %>%
  group_by(site,fever_day_F) %>%
  wilcox_test(log.wbc ~ simple_Dx) %>%
  add_xy_position(x = "simple_Dx")
stat.test$p.format <- sprintf("%.3f",stat.test$p)
stat.test$p.format[stat.test$p < 0.001] <- "<0.001"

# Plot graph
g <- ggplot(fdat2, aes(y = log.wbc, x = simple_Dx, col = simple_Dx)) +
  geom_boxplot() + 
  geom_jitter(width = 0.3, shape = 1, size = 0.5)+
  facet_grid(site ~ fever_day_F)+
  stat_pvalue_manual(
    stat.test
    , bracket.nudge.y = 0.25
    , hide.ns = F,
    label = "{p.format}"
  ) + 
  scale_y_continuous(name = expression("Log10[White blood cells count(10"^{3}~"cells/ml)]"),expand = expansion(mult = c(0.05, 0.1))) +
  xlab("Diagnosis")+
  scale_color_discrete(name = "Diagnosis")

#Save graph
ggsave(g,file = "Plots/KK_v_SK_WBC_by_FeverDay.eps",width = 10, height = 8, units = "in")
