---
title: "Obtain relative abundances and ASV composition"
author: "Maria Rebolleda-Gomez"
date: "`r Sys.Date()`"
output:
  pdf_document: default
  html_document:
    df_print: paged
---

1. Load packages. 
```r
# Knitr options
knitr::opts_chunk$set(
  cache = TRUE,
  echo = FALSE,
  fig.align = "center",
  fig.height = 4,
  fig.width = 6)

# Library
library(tidyverse)
library(data.table)
#library(DECIPHER)
```

2. Open data, format ASV table and calculate relative abundances. 
```r
otu_tab=fread("OTU_nochim_all_dada2.txt")
dim(otu_tab)
taxa_tab=fread("taxa_dada2_Sylva")
dim(taxa_tab)
mymap <-fread("metadata.csv")
dim(mymap)


#Transpose OTU table
inverse=t(otu_tab)%>%as.data.table
colnames(inverse)=as.character(inverse[1,])
inverse=inverse[-1,]

#Convert to numeric`
inverse=inverse[,lapply(inverse,as.numeric)]

#Get relative abundance of taxa 
relative_abundance=function(x,samplesrows){
  if (samplesrows==F){
    x/rep(colSums(x),each=nrow(x))
  } else {x/rowSums}
}

relab=relative_abundance(inverse,samplesrows = F)
sums=colSums(relab) #Sanity check
#One sample had no ASVs after rarefaction. To remove: 
relab=relab[,!is.na(sums),with=F]
relab[,mean_relab:=rowMeans(relab)]


#Format table (change names, make long format and split the sample names in relevant columns)
relab[,ASV:=colnames(otu_tab)[-1]]
dim(relab)
names(taxa_tab)[1]="ASV"
otu_taxa=merge(relab,taxa_tab)
otu_taxa_long=melt(otu_taxa, measure.vars =2:(ncol(otu_taxa)-7),variable.name = "Sample", value.name = "Relative_abundance")
samplesID=as.character(otu_taxa_long$Sample)%>%strsplit(split="_")%>%
  lapply(matrix, nrow=1) 
samplesID=do.call(rbind,samplesID) 
#NOTE: For some reason do.call did not work in pipe and rbindlist needs the elements of the list to be data.frames
colnames(samplesID)=c("spp","rep","date","sample_type")
otu_taxa_long=cbind(otu_taxa_long,samplesID)

#Write file
#fwrite(otu_taxa_long, "~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/ESVs.csv")   
```

3. Summarize data by family
```r
family_data=otu_taxa_long[,mean_relab:=sum(Relative_abundance),by=.(Family,Sample)]

#Format colums to order the Families according to their Class and add an Other category for families with >0.05
family_data[,TempFamily:=Family]
family_data[mean_relab<=0.05,TempFamily:="Other"]

```


3. Make plot with relative abundances. 
```r
#Re-format number of replicates and "Family" levels according to "Class"
family_data$TempFamily <- factor(family_data$TempFamily, levels = unique(family_data$TempFamily[order(family_data$Class)]))%>%
  relevel("Other")
otu_taxa_long[spp=="HETU"&rep=="6",rep:=3]

ggplot(family_data,aes(x=rep,y=Relative_abundance,fill=TempFamily, Class))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=c("gray20",'#99000d','#cb181d',
                        '#ef3b2c','#fb6a4a','#fc9272',
                        '#fcbba1','#fee5d9','#8c510a',
                        '#bf812d','#dfc27d',"#66c2a4",
                        "#addd8e",'#bdd7e7','#6baed6',
                        '#3182bd','#08519c'))+
  facet_wrap(spp~sample_type)+
  theme_bw()
       


```

4. Check if there are genera that belong only to the bottom of the petal. 
*There was no genus present in the bottom but not present in the top, however this excludes weird genera that could not be assigned, so I also checked families.*
```r
#Check by genus
bottom=otu_taxa_long$Genus[otu_taxa_long$sample_type=="B"]
rest=otu_taxa_long$Genus[otu_taxa_long$sample_type!="B"]
otu_taxa_long$Genus[!(bottom%in%rest)]

#Check by family 
bottom_fm=paste(otu_taxa_long$Family[otu_taxa_long$sample_type=="B"], 
                otu_taxa_long$Genus[otu_taxa_long$sample_type=="B"])
rest_fm=paste(otu_taxa_long$Family[otu_taxa_long$sample_type!="B"], 
                otu_taxa_long$Genus[otu_taxa_long$sample_type!="B"])
otu_taxa_long$Genus[!(bottom_fm%in%rest_fm)]

```

5. Compare with isolates
```r
#Open isolates taxonomic data
isolates=fread("../Input/Isolates_taxonomy.csv")
ngenera <- isolates$Genus %>% 
  unique %>% length
nfamilies <- isolates$`Bacteria Family` %>% 
  unique %>% length

ESVgenus=unique(otu_taxa_long$Genus)
iso_genus=unique(isolates$Genus)

simple_venndiagram=function(x,y){
  ux <- unique(x)
  uy <- unique(y)
  
  venn <- list(only_x=NULL, 
               only_y=NULL,
               shared=NULL)
  
  venn$only_x <- ux[!ux%in%uy]
  venn$only_y <- uy[!uy%in%ux]
  
  if (length(ux)>length(uy)){
    venn$shared=ux[ux%in%uy]
  } else {
    venn$shared=uy[uy%in%ux]
  }
  venn
}

genus_venn <- simple_venndiagram(otu_taxa_long$Genus, isolates$Genus)
lapply(genus_venn, length)

#Of those genera- Kluyviera, Pseudoescherichia, Pantoea, Lelliotia are not distinguishable! (see tree; I also use BLAST with the amplicons for all of these genera and they match with 100% sequences of all three genera). 

family_venn <- simple_venndiagram(otu_taxa_long$Family, isolates$`Bacteria Family`)
lapply(family_venn, length)

```

Calculate relative abundance of families
```r
iso_relab=isolates[,.(.N), by=.(`Bacteria Family`)]
iso_relab[,relab:=N/sum(N)]

esv_family=otu_taxa_long[,.(mean(mean_relab)),by=.(Family)]
esv_family[,TempFamily:=Family]
esv_family[!Family%in%iso_relab$`Bacteria Family`,TempFamily:="Other"]

iso_relab[,TempFamily:=`Bacteria Family`]
iso_relab[!`Bacteria Family`%in%esv_family$Family,TempFamily:="Other"]

plot_data=data.frame(Family=c(iso_relab$TempFamily,esv_family$TempFamily),
                     Relab=c(iso_relab$relab, esv_family$V1),
                     Type=c(rep("cultivable",nrow(iso_relab)),
                            rep("amplicon",nrow(esv_family))))
plot_data$Family=factor(plot_data$Family, c("Other","Nocardiaceae",
                                                "Microbacteriaceae","Rhizobiaceae",
                                                "Bacillaceae","Pseudomonadaceae",
                                                "Enterobacteriaceae","Xanthomonadaceae"))

ggplot(plot_data,aes(x=Type,y=Relab,fill=Family))+
  geom_bar(stat="identity",color="black")+
  ylab("Relative abundance")+
  scale_fill_manual(values=c("gray20",'#99000d','#cb181d',
                        '#bf812d','#41ae76','#bdd7e7',
                        '#6baed6','#08519c'))+
  theme_bw()


data_chi=as.data.table(plot_data)
data_chi=data_chi[,.(sum(Relab)),by=.(Type,Family)]
data_chi=dcast(data_chi,Family~Type, value.var="V1")
data_chi=data_chi[,2:3]%>%as.matrix
data_chi[,1]=data_chi[,1]*nrow(otu_taxa_long)
data_chi[,2]=data_chi[,2]*90

chisq.test(t(data_chi))
```
