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

merged.data <- merge(baseline, safety, by="MASKID")
merged.data <- merge(merged.data, outcomes, by="MASKID")
merged.data <- merge(merged.data, sites, by="NEWSITEID", all.x = T)

data.with.site <- merged.data
net.levels <- levels(data.with.site$NET)
data.with.site$NET <- ifelse(is.na(data.with.site$NET), "MIS", net.levels[data.with.site$NET])
data.with.site$NET <- as.factor(data.with.site$NET)

print(plyr::ddply(data.with.site, ~NET, count))

survival.outcome <- with(data.with.site, Surv(time=SAE_DAYS, event=SAE_EVNT))
cox.fit <- with(data.with.site, coxph(survival.outcome ~ INTENSIVE * NET))
print(summary(cox.fit))

treated.arm <- data.with.site[data.with.site$INTENSIVE==1,]
control.arm <- data.with.site[data.with.site$INTENSIVE==0,]
rate <- with(treated.arm, tapply(SAE_EVNT, NET, mean))
rate <- rate - with(control.arm, tapply(SAE_EVNT & !INTENSIVE, NET, mean))

se <- with(treated.arm, tapply(SAE_EVNT, NET, sd) / sqrt(tapply(SAE_EVNT, NET, NROW)))
se <- sqrt(se^2 + with(control.arm, tapply(SAE_EVNT, NET, sd) / sqrt(tapply(SAE_EVNT, NET, NROW)))^2)

network <- rownames(se)
setEPS()
postscript("sae_network_interaction.eps", height=6, width=8)
require(plotrix)
par(mar=c(5,5,3,3))
plotCI(1:5, rate, ui=rate+1.28*se, li=rate-1.28*se, ylim=c(-0.1, 0.1), xlab="Network", ylab="Change in SAE incidence with treatment", main="Network vs effect of treatment on SAE",
       axes=F, mgp=c(3,0,0))
axis(2)
axis(1, at=1:5, labels=as.character(network), las=1)
dev.off()
