---
title: "Biolog NMDS Plot"
author: "Rebecca Hayes"
date: "12/2/2019"
output: html_document
editor_options: 
  chunk_output_type: console
---

### Load packages and data

```r
library(tidyr)
library(vegan)
library(dplyr)
library(ggpubr)
library(viridis)
library(cowplot)
library(ggrepel)

```

Data is available in [Dryad](https://doi.org/10.5061/dryad.2v6wwpzjt). 

```r

data.long <- read.csv("CC2018_Biolog Data_12-2-2019 RH.csv")
uv <- read.csv("CC2018_LD50 Calculator pooled_12-2-2019 RH.csv")
gro <- read.csv("Growth curve data.csv")

#remove unnecessary cols
uv <- uv %>%
  dplyr::select(-X,-Strain.family, -Strain.species, -Flower,-Transect, -se)

```

### Format data
Merge row and column into one, remove separate columns
```r

#remove text columns and filters to include only the carbon assays
carb.long <- data.long %>%
  dplyr::select(-row,-column, -abs550, -qualitative,-X) %>%
  filter(Type == "Carbon")

#remove text columns and filters to include only the chemical assays
chem.long <- data.long %>%
  dplyr::select(-row,-column, -abs550, -qualitative,-X) %>%
  filter(Type == "Chemical")


```

Transform data from long to wide to make table with Strain.ID as rows and Assay as columns
```r

#change carb.long to wide format
carb.wide <- carb.long %>%
    group_by(Strain.ID) %>%
    spread(key = Assay, value = abs550_std)

#change carb.long to wide format
chem.wide <- chem.long %>%
    group_by(Strain.ID) %>%
    spread(key = Assay, value = abs550_std)

#preserve text columns in own object
meta.carb <- carb.wide[,1:6]
meta.chem <- chem.wide[,1:6]


```

### Make nmds plots 
Plot with carbon sources

```r



#set negative values to zero
carb.wide[carb.wide < 0] <- 0

#create matrix from cols 7:77 (only assay values)
carb.comm <- carb.wide[,7:77]

#run NMDS
carb.nmds <- metaMDS(carb.comm,
          k = 2,
          maxit = 999,
          trymax = 250,
          wascores = TRUE, na.rm=TRUE)


####numbered steps must be done in order

#1. plot just points
plot(carb.nmds$points)


#make data frame from the NMDS points
points <- as.data.frame(carb.nmds$points)


carb.nmdv <- envfit(carb.nmds$points, carb.comm) #create ordination
carb.arrows <- carb.nmdv$vectors$arrows #pull out arrows
carb.arrows$carbon <- rownames(carb.arrows) #name all vectors for the assay type
carb.arrows <- as.data.frame(carb.nmdv$vectors$arrows*sqrt(carb.nmdv$vectors$r)) #use only really significant arrows
carb.arrows$p.values <- carb.nmdv$vectors$pvals #add in column for p vals
carb.arrows$carbon <- rownames(carb.arrows)
points$family <- as.character(meta.carb$Strain.family) #add in family metadata
points$Strain.ID <- as.character(meta.carb$Strain.ID) #add in strain id meta
points$host <- as.character(meta.carb$Flower) #add in host meta
points$trans <- as.character(meta.carb$Transect) #add in transect meta


#make scatter plot add arrows where p<0.01 by strain family
plot <- ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$family, size=12)) +
  geom_segment(data=carb.arrows[carb.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=carb.arrows[carb.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=carbon, size=12), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Strain Family") +
   theme(legend.position = "left", legend.text=element_text(size=30),legend.title=element_text(size=36))

View(legend)
##get legend
legend <- get_legend(plot)
as_ggplot(legend)
legend.text=element_text(size=30),legend.title=element_text(size=36)

#make scatter plot add arrows where p<0.01 by transect
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$trans, size=1.5)) +
  geom_segment(data=carb.arrows[carb.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=carb.arrows[carb.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=carbon), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Transect")


#make scatter plot add arrows where p<0.01 by host
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$host, size=1.5)) +
  geom_segment(data=carb.arrows[carb.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=carb.arrows[carb.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=carbon), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Host")


```

Chemical NMDS
```r



#set negative values to zero
chem.wide[chem.wide < 0] <- 0

#create matrix from cols 7:29 (only assay values)
chem.comm <- chem.wide[,7:29]

#run NMDS
chem.nmds <- metaMDS(chem.comm,
          k = 2,
          maxit = 999,
          trymax = 250,
          wascores = TRUE, na.rm=TRUE)
####numbered steps must be done in order

#1. plot just points
#plot(chem.nmds$points)

#2. add sites (only numbers atm)
#orditorp(chem.nmds$points,display="sites",cex=1.25,air=0.01)

#2.make spiderplot
#with(meta.data, ordispider(chem.nmds$points, Transect, kind = "se", conf = 0.95))

#make data frame from the NMDS points
points <- as.data.frame(chem.nmds$points)

 

chem.nmdv <- envfit(chem.nmds$points, chem.comm)
chem.arrows <- chem.nmdv$vectors$arrows
chem.arrows <- as.data.frame(chem.nmdv$vectors$arrows*sqrt(chem.nmdv$vectors$r))
chem.arrows$p.values <- chem.nmdv$vectors$pvals
chem.arrows$chem <- rownames(chem.arrows)

#add in revelant metadata
points$family <- as.character(meta.chem$Strain.family)
points$Strain.ID <- as.character(meta.chem$Strain.ID)
points$host <- as.character(meta.chem$Flower)
points$trans <- as.character(meta.chem$Transect)

#make scatter plot add arrows where p<0.01 by strain family
plot <- ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$family, size=12)) +
  geom_segment(data=chem.arrows[chem.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=chem.arrows[chem.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=chem, size=12), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Strain Family") +
   theme(legend.position = "left", legend.text=element_text(size=30),legend.title=element_text(size=36))


plot
##get legend
legend <- cowplot::get_legend(plot)
as_ggplot(legend)


#make scatter plot add arrows where p<0.01 by transect
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$trans, size=1.5)) +
  geom_segment(data=chem.arrows[chem.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=chem.arrows[chem.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=chem), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Transect")

#make scatter plot add arrows where p<0.01 by host
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$host, size=1.5)) +
  geom_segment(data=chem.arrows[chem.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=chem.arrows[chem.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=chem), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Host")


```


Adding in UV data
```r

points$uv <- as.numeric(uv$LD50)
###only do once! removes data from row 17 (D11) b/c no UV assay
points <- points[-c(17),] 
points <- merge(points, uv, by="Strain.ID")


###only do once! removes data from row 17 (D11) b/c no UV assay



points$uv <- as.numeric(uv.metachem$LD50)


#Colored by UV tolerance

#for carbon
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$LD50, size=1)) +
  geom_segment(data=carb.arrows[carb.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=carb.arrows[carb.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=carbon), inherit.aes=FALSE) +
  labs(color = "LD50 (sec)") + scale_color_viridis(discrete=FALSE, direction=-1) + guides(size=FALSE)


#for chemical
ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$LD50, size=1)) +
  geom_segment(data=chem.arrows[chem.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=chem.arrows[chem.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=chem), inherit.aes=FALSE) +
  labs(color = "LD50 (sec)") + scale_color_viridis(discrete=FALSE, direction=-1) + guides(size=FALSE)


```

```r
#5-8 combined all traits NMDS



#this sheet has carbon and chemical assays combined
dat.wide <- read.csv("CC2018_Wide form biolog data_3.csv")


names(dat.wide)[names(dat.wide) == "X3.Methyl.Glucose"] <- "3 Methyl Glucose"

names(dat.wide)[names(dat.wide) == "X1..NaCl"] <- "1% NaCl"

names(dat.wide)[names(dat.wide) == "X4..NaCl"] <- "4% NaCl"

names(dat.wide)[names(dat.wide) == "X8..NaCl"] <- "8% NaCl"

dat.wide$X4..NaCl

#preserve text columns in own object
meta.dat <- dat.wide[,1:7]


#create matrix from cols 7:77 (only assay values)
dat.comm <- dat.wide[,8:100]



#run NMDS
dat.nmds <- metaMDS(dat.comm,
          k = 2,
          maxit = 999,
          trymax = 250,
          wascores = TRUE, na.rm=TRUE)

#make data frame from the NMDS points
points <- as.data.frame(dat.nmds$points)

dat.nmdv <- envfit(dat.nmds$points, dat.comm) #create ordination
dat.arrows <- dat.nmdv$vectors$arrows #pull out arrows
dat.arrows$carbon <- rownames(dat.arrows) #name all vectors for the assay type
dat.arrows <- as.data.frame(dat.nmdv$vectors$arrows*sqrt(dat.nmdv$vectors$r)) #use only really significant arrows
dat.arrows$p.values <- dat.nmdv$vectors$pvals #add in column for p vals
dat.arrows$carbon <- rownames(dat.arrows)
points$family <- as.character(meta.dat$Strain.family) #add in family metadata
points$Strain.ID <- as.character(meta.dat$Strain.ID) #add in strain id meta
points$host <- as.character(meta.dat$Flower) #add in host meta
points$trans <- as.character(meta.dat$Transect) #add in transect meta


#make scatter plot add arrows where p<0.01 by strain family
plot <- ggscatter(data = points, x = "MDS1", y = "MDS2") +
  geom_point(data=points, aes(color=points$family, size=12, shape=points$host)) +
  geom_segment(data=dat.arrows[dat.arrows$p.values < 0.01,],
               mapping = aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),alpha = 0.4, inherit.aes = FALSE) + geom_text_repel(data=dat.arrows[dat.arrows$p.values < 0.01,],
            aes(x=MDS1,y=MDS2,label=carbon, size=12), inherit.aes=FALSE) +
  guides(size=FALSE, color = guide_legend(override.aes = list(size = 5)), shape = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Strain Family", shape="Host") +
   theme(legend.position = "left", legend.text=element_text(size=12),legend.title=element_text(size=16))

plot

ggsave("plot.tiff", dpi=300, dev='png', height=6, width=8, units="in")

```

