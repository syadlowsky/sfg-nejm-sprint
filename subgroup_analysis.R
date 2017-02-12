library(causalTree)
library(dplyr)

basedir <- "~/phd/sfg/nejm-sprint"
setwd(basedir)

figdir <- file.path(basedir, "figs")

baseline <- read.csv("data/baseline.csv")
safety <- read.csv("data/safety.csv")
outcomes <- read.csv("data/outcomes.csv")
merged.data <- merge.data.frame(baseline, safety, by="MASKID")
merged.data <- merge.data.frame(merged.data, outcomes, by="MASKID")

outcome <- merged.data$EVENT_PRIMARYORDEATH & merged.data$T_PRIMARYORDEATH <= 365.25*5
mean(outcome)
tree <- causalTree(outcome ~ RACE4 + FEMALE + BMI + ASPIRIN + SMOKE_3CAT + STATIN, data = merged.data, treatment = merged.data$INTENSIVE,
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                   xval = 5, cp = 0, minsize = 50, propensity = 0.5)
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
opfit <- prune(tree, opcp)
svg(file.path(figdir, "honest_tree_pruned.svg"), height=6, width=8)
rpart.plot(opfit)

# See how this relationship looks
merged.data$bmi.quantile <- ntile(merged.data$BMI, 10)
rate <- tapply(outcome & merged.data$INTENSIVE, merged.data$bmi.quantile, mean)
rate <- rate - tapply(outcome & !merged.data$INTENSIVE, merged.data$bmi.quantile, mean)
bmi <- tapply(merged.data$BMI, merged.data$bmi.quantile, mean)
svg(file.path(figdir, "ate_bmi_interaction.svg"), height=6, width=8)
plot(bmi, rate, xlab="BMI", ylab="Average effect of intensive treatment", main="BMI vs effect of treatment on primary outcome")
