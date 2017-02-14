library(survival)
library(dplyr)

basedir <- "~/phd/sfg/nejm-sprint"
setwd(basedir)

figdir <- file.path(basedir, "figs")
setwd("~/phd/sfg/nejm-sprint")

baseline <- read.csv("data/baseline.csv")
safety <- read.csv("data/safety.csv")
outcomes <- read.csv("data/outcomes.csv")
sites <- read.csv("sites/networkXid.csv")

merged.data <- merge.data.frame(baseline, safety, by="MASKID")
merged.data <- merge.data.frame(merged.data, outcomes, by="MASKID")
merged.data <- merge.data.frame(merged.data, sites, by="NEWSITEID", all.x = T)

data.with.site <- merged.data
net.levels <- levels(data.with.site$NET)
data.with.site$NET <- ifelse(is.na(data.with.site$NET), "MIS", net.levels[data.with.site$NET])
data.with.site$NET <- as.factor(data.with.site$NET)

print(plyr::ddply(data.with.site, ~NET, count))

survival.outcome <- with(data.with.site, Surv(time=SAE_DAYS, event=SAE_EVNT))
cox.fit <- with(data.with.site, coxph(survival.outcome ~ INTENSIVE * NET))
print(summary(cox.fit))

sae.events.all.patients <- merged.data[, c("REL_SAE_EVNT", "HYP_SAE_EVNT", "SYN_SAE_EVNT", "BRA_SAE_EVNT",
                                               "ELE_SAE_EVNT", "INJ_SAE_EVNT", "AKI_SAE_EVNT", "HYP_ERS_EVNT", "SYN_ERS_EVNT",
                                               "BRA_ERS_EVNT", "ELE_ERS_EVNT", "INJ_ERS_EVNT", "AKI_ERS_EVNT", "LON_MCE_EVNT",
                                               "HIN_MCE_EVNT", "LOK_MCE_EVNT", "HIK_MCE_EVNT", "ALO_OHO_EVNT", "WDZ_OHO_EVNT")]

southeast <- !is.na(merged.data$NET) & merged.data$NET=="SOUTHEAST"
southeast.patients <- sae.events.all.patients[southeast, ]
(sum(sae.events.all.patients)/sum(southeast.patients)) * colSums(southeast.patients) / colSums(sae.events.all.patients)