rm(list=ls())
###Question 1###
library(survival)

#read data
setwd("Desktop/Grad2Winter/NJEM/sprint_pop/data")
outcome = read.csv("outcomes.csv",header=TRUE)
baseline = read.csv("baseline.csv",header=TRUE)

#check if the maskID is the same
for (i in 1:dim(outcome)[1]){
  if (outcome[i,1]!=baseline[i,1]){
    print (i)
  }
}

#extract the intensive treatment/primary outcome/T_days columns for likelihood ratio
primary = data.frame(cbind(baseline$MASKID,baseline$INTENSIVE,outcome$EVENT_PRIMARY,outcome$T_PRIMARY))
colnames(primary) = c("ID","Intensive","Event_primary","T_primary")

#Fit COX PH model
model1 = coxph(Surv(T_primary,Event_primary) ~ Intensive, data=primary)




###Question 2####
bp = read.csv("bp.csv")
bp = bp[,1:3]
bp2 = cbind(bp, rep(0,nrow(bp)))
colnames(bp2) = c("ID","VISITCODE","sbp","Intensive")

for (i in 1:length(unique(bp$MASKID))){
  bp2[bp2$ID==baseline$MASKID[i],]$Intensive = rep(baseline$INTENSIVE[i],nrow(bp2[bp2$ID==baseline$MASKID[i],]))
}

bp2 = na.omit(bp2)
month = rep(0,nrow(bp2))
#transform factors to numeric for comparison
for (i in 1:nrow(bp2)){
  if (bp2$VISITCODE[i]=="RZ"){
    month[i] = 0
  }
  else{
    level = levels(bp2$VISITCODE)[bp2[i,2]]
    month[i] = as.numeric(strsplit(level,"M"))
  }
}

bp3 = data.frame(cbind(bp2$ID,month,bp2$sbp,bp2$Intensive))
colnames(bp3) = c("ID","month","sbp","Intensive")

#select the last visit for each patient
bplast = c()
for (i in unique(bp3$ID)){
  extract = bp3[bp3$ID==i,]
  lastVisit = extract[extract$month==max(extract$month),]
  bplast = rbind(bplast,lastVisit)
}

#get those who go to follow-up sbp checks
bplast2 = bplast[bplast$month!=0,]

#calculate the sample size and mean sbp for the two arms
size1 = nrow(bplast2[bplast2$Intensive==1,])
size0 = nrow(bplast2[bplast2$Intensive==0,])
mean1 = mean(bplast2[bplast2$Intensive==1,]$sbp)
mean0 = mean(bplast2[bplast2$Intensive==0,]$sbp)

