library(survival)

setwd("~/phd/sfg/nejm-sprint")

baseline <- read.csv("data/baseline.csv")
safety <- read.csv("data/safety.csv")
outcomes <- read.csv("data/outcomes.csv")

merged.data <- merge.data.frame(baseline, safety, by="MASKID")
merged.data <- merge.data.frame(merged.data, outcomes, by="MASKID")

survival.outcome <- Surv(time=merged.data$SAE_DAYS, event=merged.data$SAE_EVNT)
cox.fit <- coxph(survival.outcome~merged.data$INTENSIVE * factor(merged.data$NEWSITEID))
summary(cox.fit)
