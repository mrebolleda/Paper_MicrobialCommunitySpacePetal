---
paper: ""
title: "Growth Curves"
author: "Maria Rebolleda-Gomez (mariarebolleda@gmail.com)"
date: "June 2nd 2020"
output: html_document
---

Load packages and open data

```r

# Library
library(tidyverse)
library(data.table)
library(mgcv)
library(grofit)
library(emmeans)


#Open growth curves data as data.table
dt=fread("GrowthCurves_Analysis/GC_noblank_LowHiSucrose.csv")

#Open metadata
metadata=fread("GrowthCurves/input/GC_metadata.csv")
```

Re-format data.table to long format and add strain names, low and hi sucrose columns. 
```r
dt.long=melt(dt, measure.vars = 3:ncol(dt))
dt.long[, c("Strain", "Concentration") := tstrsplit(variable, "-", fixed=TRUE)]

#Function to calculate time intervals from time in hr:min:sec format
time_min=function(time){
  n=length(time)
  M=strsplit(time, ":")%>%
  unlist%>%as.numeric%>%matrix(ncol=n)

  #Convert hours and seconds to minutes
  M[1,]=M[1,]*60
  M[3,]=M[3,]/60
  
  #Add all minutes and remove initial time from all 
  time.min=colSums(M)
  time.min-time.min[1]
}


dt.long[, Minutes:=time_min(Time),variable]

```

Plot data for exploration
```r
ggplot(dt.long, aes(x=Minutes, y=value, group=variable, color=Concentration))+
  geom_line()+
  theme_bw()

#There are variable initial absorbance values, to normalize and be able to compare better, remove value at time 0

dt.long[,norm_value:=value-value[Minutes==0],variable]
```


Format data for analysis
```r
#Merge data and metadata
setnames(dt.long, "Strain", "Strain.ID")
data_all=merge(dt.long[,variable:norm_value],metadata,all.x = T)

#Identify curves with no significant growth (i.e. OD600(48hrs)<=0.01) 
lowgrowth=data_all$variable[data_all$Minutes==max(data_all$Minutes)&data_all$norm_value<=0.01]
data_all_growth=data_all[!data_all$variable%in%lowgrowth,]

ggplot(data_all_growth, aes(x=Minutes, y=norm_value, group=variable, color=Concentration))+
  geom_point()+
  theme_bw()

#Re-format data for grofit
samples=length(unique(data_all_growth$variable))
time=t(matrix(data_all_growth$Minutes,ncol=samples))

dt_grofit=dcast(data_all_growth, Strain.ID+variable +Concentration ~ Minutes, value.var = "norm_value")

```

Fit growth curves with `gcfit` and plot curves with fitted models 
```r
growth_results=gcFit(time,dt_grofit)
gc_results=growth_results$gcTable
write.csv(gc_results,"~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/GrowthCurves/Output/GC_fit_summary.csv")

#Loop to extract values from model in order to plot:
#1. Create empty list 
temp=vector(mode = "list",samples)

#2. Save number of time points to use in the loop
times=dim(time)[2]

#Loop through each sample to get the fitted data (model and spline)
for (i in 1:samples){
  if (length(growth_results$gcFittedModels[[i]]$fit.data)==0){
    temp[[i]]=data.table(fit_model=NA,
                     fit_spline=growth_results$gcFittedSplines[[i]]$fit.data)
    warning("No model fit, replaced with NA")
  } else {
     temp[[i]]=data.table(fit_model=growth_results$gcFittedModels[[i]]$fit.data,
                     fit_spline=growth_results$gcFittedSplines[[i]]$fit.data)
  }
  temp[[i]][,Strain.ID:=rep(growth_results$gcFittedModels[[i]]$gcID[1],times)]
  temp[[i]][,Concentration:=rep(growth_results$gcFittedModels[[i]]$gcID[3],times)]
  temp[[i]][,Reliability:=rep(growth_results$gcFittedModels[[i]]$reliable,times)]
  temp[[i]][,Model:=rep(growth_results$gcFittedModels[[i]]$model,times)]
}

#Rbind through list to have one data frame
dt.fit=rbindlist(temp)

#Add columns with time and data
dt.fit[,Time:=data_all_growth$Minutes] 
dt.fit[,Data:=data_all_growth$norm_value] 
dt.fit[,Family:=data_all_growth$Strain.family] 


#Plot data:
fams_colors=c("#a6d96a","#d73027","#35978f","#fdae61","#bf812d")
ggplot(dt.fit, aes(x=Time/60, y=Data, group=interaction(Strain.ID,Concentration), color=Family))+
  facet_grid(~Concentration)+
  xlab("Time (hours)")+
  ylab("OD (600nm)")+
  geom_point(size=0.5, alpha=0.5)+
  geom_line(aes(y=fit_spline), linetype=5)+
  geom_line(aes(y=fit_model),size=1)+
  scale_color_manual(values=fams_colors)+
  theme_bw()

```

Determine growth rate to use- if no parametric model was able to fit the data, then use splines. 
```r
gc_results=as.data.table(gc_results)
gc_results[,mu.cons:=mu.model]
gc_results[is.na(mu.model),mu.cons:=mu.spline]
```

Determine the impact of petal location on growth rates 
```r
#Add back no-growth data
data_all_sm=data_all[Minutes==0,c(-3,-5,-6)]
setnames(gc_results, "concentration", "Concentration")
setnames(gc_results, "TestId", "Strain.ID")
growth_meta_all=gc_results%>%select(c("Strain.ID":"integral.spline","mu.cons"))%>%
  merge(data_all_sm,by=c("Strain.ID","Concentration"),all.y = T)
growth_meta_all[is.na(growth_meta_all$mu.cons),"mu.cons"]=0

#Remove E.coli control
growth_meta_all=growth_meta_all[Strain.ID!="EC",]

ggplot(growth_meta_all, aes(x=Transect, y=mu.cons, colour=Flower))+
  facet_grid(~Concentration)+
  xlab("Transect")+
  ylab("Growth rate")+
  geom_point(size=1)+
  stat_smooth(method="lm",se=T)+
  #stat_summary(geom="point",size=2.5, fun.y=mean)+
  scale_color_manual(values=c("#cb181d","#006ba4"))+
  theme_bw()

m1=lm(mu.cons~Transect+Flower+Concentration+Transect:Flower+Transect:Concentration+Flower:Concentration, data=growth_meta_all)
summary(m1)


marg_means=emtrends(m3,c("Flower","Concentration"),var="Transect")
ggplot(as.data.frame(marg_means), aes(x=paste(Flower,Concentration), y=Transect.trend, color=Flower, shape=Concentration))+
  geom_point()+
  geom_errorbar(aes(ymin=Transect.trend-SE,ymax=Transect.trend+SE), width=0.2)+
  geom_hline(yintercept = 0, linetype=2)+
  scale_color_manual(values=c("#cb181d","#006ba4"))+
  scale_shape_manual(values=c(19,21))+
  ylab("Slope")+
  geom_point()+
  theme_bw()


```





