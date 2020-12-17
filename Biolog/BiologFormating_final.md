---
title: "BiologAnalysis_20191126"
author: "Maria Rebolleda-Gomez"
date: "26/11/2019"
output: html_document
---
Load packages
```{r}
library(tidyr)
library(dplyr)
```

Open files
```{r}
#1. Write your pathway to files
setwd("C:/Users/rah10/Dropbox/Flower-microbe Color Curiculum Dev/Analyses/Biolog/Input/")

files="C:/Users/rah10/Dropbox/Flower-microbe Color Curiculum Dev/Analyses/Biolog/Input/"

#2. Load in carbon list
carbon <- read.csv(file="Carbon List.csv")


#3. Load in strain list
strain_list <- read_xlsx("C:/Users/rah10/Dropbox/Flower-microbe Color Curiculum Dev/Analyses/Biolog/Input/Strain List.xlsx")

  #remove unnecessary cols
  strain_list$Per.ID=NULL
  strain_list$Notes=NULL
```

Function to read each file and transform into long table and calculate standarized value
```{r}
# Create dataframe
data_all=c("A","13",0,"DELETE","DELETE",0)

#for each strain in strain_list, read files that start "Biolog_" strain.ID ".csv" and make a table
for (i in strain_list$Strain.ID){
  biologtable=read.csv(paste(files,"Biolog_",i,".csv",sep=""))
biolog_long <- gather(biologtable, column, abs550, X1:X12, factor_key=TRUE)
  biolog_long$Strain=rep(i,nrow(biolog_long))
  minuscontrol=biolog_long$abs550-biolog_long$abs550[biolog_long$X=="A"&biolog_long$column=="X1"]
  is0=minuscontrol==0
  biolog_long$abs550_std=0
  biolog_long$abs550_std[!is0&!is.na(is0)]=minuscontrol[!is0&!is.na(is0)]/minuscontrol[biolog_long$X=="A"&biolog_long$column=="X10"]
  data_all=rbind(data_all,biolog_long)
}


data_all$qualitative="borderline"
data_all$qualitative[data_all$abs550_std<=0]="negative"
data_all$qualitative[data_all$abs550_std>=1]="positive"


#Change abs550 to numeric
data_all$abs550=as.numeric(data_all$abs550)


```

Change format to match carbon sources and merge with carbon file
```{r}
colnames(data_all)[1]="row"

#remove empty 1st row
data_all <- data_all[-c(1),]

#change "Strain" to "Strain.ID"
colnames(data_all)[4]="Strain.ID"

data_all$column=substr(biolog_long$column,2,3)

head(data_all)
```


merge carbon source and strain data
```{r}
data_all_carbon=merge(data_all,carbon)

data_final=merge(strain_list, data_all_carbon)
```




