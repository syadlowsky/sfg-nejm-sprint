library(survival)
library(causalTree)
library(dplyr)

setwd("~/phd/sfg/nejm-sprint")
baseline <- read.csv("data/baseline.csv")
safety <- read.csv("data/safety.csv")
outcomes <- read.csv("data/outcomes.csv")
merged.data <- merge.data.frame(baseline, safety, by="MASKID")
merged.data <- merge.data.frame(merged.data, outcomes, by="MASKID")

outcome <- Surv(merged.data$SAE_DAYS, merged.data$SAE_EVNT)
interaction.model <- coxph(outcome ~ INTENSIVE*RACE4 + INTENSIVE*FEMALE + INTENSIVE*BMI + INTENSIVE*ASPIRIN + INTENSIVE*SMOKE_3CAT + INTENSIVE*STATIN,
                           data=merged.data)

tree <- causalTree(SAE_EVNT ~ RACE4 + FEMALE + BMI + ASPIRIN + SMOKE_3CAT + STATIN, data = merged.data, treatment = merged.data$INTENSIVE,
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                   xval = 5, cp = 0, minsize = 50, propensity = 0.5)
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
svg("honest_tree_pruned.svg", height=6, width=8)
rpart.plot(opfit)
dev.off()
# See how this relationship looks
merged.data$bmi.quantile <- ntile(merged.data$BMI, 10)
rate <- tapply(merged.data$SAE_EVNT & merged.data$INTENSIVE, merged.data$bmi.quantile, mean)
rate <- rate - tapply(merged.data$SAE_EVNT & !merged.data$INTENSIVE, merged.data$bmi.quantile, mean)

se <- tapply(merged.data$SAE_EVNT & merged.data$INTENSIVE, merged.data$bmi.quantile, sd)
se <- (se + tapply(merged.data$SAE_EVNT & !merged.data$INTENSIVE, merged.data$bmi.quantile, sd))/sqrt(tapply(merged.data$SAE_EVNT & !merged.data$INTENSIVE, merged.data$bmi.quantile, NROW))

bmi <- tapply(merged.data$BMI, merged.data$bmi.quantile, mean)
setEPS()
postscript("sae_bmi_interaction.eps", height=6, width=8)
require(plotrix)
plotCI(bmi, rate, ui=rate+1.96*se, li=rate-1.96*se, xlab="BMI", ylab="Change in SAE incidence with treatment", main="BMI vs effect of treatment on primary outcome")
dev.off()
