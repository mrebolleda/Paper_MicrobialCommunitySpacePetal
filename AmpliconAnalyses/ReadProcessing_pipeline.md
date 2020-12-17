---
title: "Dada_diversity"
author: "Maria Rebolleda-Gomez"
date: "2019/01/29"
output: html_document
---

```r
library("vegan")
library("phyloseq")
library("ggplot2")
library("tidyr")
library(plyr)
library(VennDiagram)
library(picante)
library(ape)
library(DECIPHER)
library(phangorn)
```

# Part 1: Dada 2 pipeline
Data was first demultiplex with idemp (https://github.com/yhwu/idemp). Amplicon data is available in the genebank short read repository. 

### ASV table and analyses
Check quality:
```r
library(dada2)
path <- #Set path to reads
fnFs <- sort(list.files(path, pattern="L001_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="L001_R2", full.names = TRUE))

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
```

Get sample names
```r
setwd() # Set to working directory
metadata=read.table("Input/190628_Ashman_16S_KEH_190626_edits_CC2018.txt",header = T, sep = "\t")
sample.names=substr(fnFs,regexpr("16S_", fnFs, fixed=TRUE)[1]+4,regexpr("_S", fnFs, fixed=TRUE)[1]-1)
```

Filter and trim minimizing errors but mantaining most reads
```r
filt_path <- file.path("path_to_filtered", "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,240),
              maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
              compress=TRUE)

```

Learn error rate. 
```r
errF <- learnErrors(filtFs, multithread=TRUE, randomize=T)
errR <- learnErrors(filtRs, multithread=TRUE,  randomize=T)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

Dereplication to obtain "unique sequences"
```r
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

```

Sample inference: infer the core sequence variants using the error models learned from the data.

```r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR,  multithread=TRUE)
```

Merge pair ends
```r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
```


Make a sequence table (something like a higher resolution OTU table)
```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Check sequence lengths (must be around 250bp)
hist(nchar(getSequences(seqtab)))
#Percentage of sequences longer than 256 --> 0%
sum(nchar(getSequences(seqtab))>310)/dim(seqtab)[2]*100
##Percentage of sequences shorter than 248 --> 0.6% !!!!
sum(nchar(getSequences(seqtab))<260)/dim(seqtab)[2]*100

#Remove sequences much shorter or longer than expected
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(260,310)]
```

Remove chimeras
```r
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```


Get counts of reads at different points in the process. 
```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
matplot(t(track), type = "l")

track=as.data.frame(track)
track$PercentLoss=(track$input-track$nonchim)/track$input*100
track[track$PercentLoss>=20,]
```

Assign OTUs
```r
#taxa <- assignTaxonomy(seqtab.nochim, "raw/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "raw/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

Load as a phylloseq object to remove chloroplast
```r
metadata$SampleID=as.character(gsub("\\.", "_", metadata$SampleID))
metadata=metadata[match(rownames(seqtab.nochim),(metadata$SampleID)),]
rownames(metadata)=rownames(seqtab.nochim)

library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

#Checks if there are any mitochondria or chloroplast reads, and if there are it filters them out. 
#Caution: It also filters out unknown orders or families "NAs"
if (sum(taxa.print[,4]=="Chloroplast"|taxa.print[,5]=="Mitochondria",na.rm=T)){
  ps <- subset_taxa(ps,Order!="Chloroplast"&Family!="Mitochondria")} else{
    message("No chloroplast or mitochondria found")
  }


```

Save files
```r
write.table(seqtab.nochim,"/Output/OTU_nochim_all_dada2.txt")
write.table(taxa,"/Output/taxa_dada2_Sylva")
write.csv(metadata, "metadata.csv")

```

# Part 2: ESV table filtering and formating
For this part it is not necessary to download the reads. Data files are available in this repository. 
```r
otu_tab=read.table("/Output/OTU_nochim_all_dada2.txt")
taxa_tab=read.table("/Output/taxa_dada2_Sylva")
mymap <- read.csv("metadata.csv")
```

Construct tree for dada2 data
```r
seqs <- colnames(otu_tab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
```


```r

colHETU=c("#fcbba1","#ef3b2c","#a50f15")
colVEAL=c("#a2c8ec","#5f9ed1","#014a70")
rarecurve(otu_tab, label=F, col=c(rep(colHETU,4),rep(colVEAL,4)))
abline(v=1500, lty=2)
```

## Data Filtering 
Remove low abundance reads and median standarize samples. 

```r
rownames(mymap)=mymap$SampleID
ps=phyloseq(otu_table(otu_tab, taxa_are_rows=FALSE), 
               sample_data(mymap), 
               tax_table(taxa),
            phy_tree(fitGTR$tree))

#saveRDS(ps,"phyloseq_20191024.RDS")

#Calculate median number of reads/sample as well as min 
sample.reads <- rowSums(otu_table(ps))
median(sample.reads)
#12122
min(sample.reads[sample.reads>0])
#1237

#Remove samples with very low read count. 
plot(sample.reads, sample.reads, xlab="Number of reads", ylab="Number of reads")
abline(h=5000, lty=2)
#It looks like there is a jump between ~2400 to ~9100 reads per sample, and of the three samples with less than ~2400 one has 0 reads. 

#Check which are the samples with low reads
sample.reads[sample.reads<5000]
#HETU 6 B
#VEAL 3 B
#VEAL 5 T

#Remove very low samples
ps_pruned = prune_samples(sample_sums(ps)>=2500, ps)

##Sanity check
which(rowSums(otu_table(ps))<2500)
which(rowSums(otu_table(ps_pruned))<2500)

#Remove rare OTUs 
relab=otu_table(ps_pruned)/rowSums(otu_table(ps_pruned)) 
keep=colMeans(relab)>1e-5
names(keep)=NULL
ps_clean=prune_taxa(keep, ps_pruned)
metadata=as.data.frame(sample_data(ps_clean))

##Sanity check
ncol(otu_table(ps_pruned))-sum(colMeans(relab)<=1e-5)==ncol(otu_table(ps_clean))

#Median standarize samples
med = median(rowSums(otu_table(ps_clean)))
standf = function(x, m=med) round(m * (x / sum(x)))
ps_clean_std = transform_sample_counts(ps_clean, standf)#

#saveRDS(ps_clean_std,"phyloseq_20191025_medianstd.RDS")

```


# Part 3: Calculate beta and alpha diversity
```r

#This function calculates distance matrices with four different diversity indeces and
#outputs a list were each element of the list is one of the diversity matrices. 

dist=function(ps){
  #Check if OTUs are rows or columns- for vegan has to be cols. 
  if (dim(otu_table(ps))[1]>dim(otu_table(ps))[2]){
    warning("There are more rows than columns, check that samples are row and not columns")
  }
  
  mat=as.matrix(otu_table(ps)@.Data)
  colnames(mat)=NULL
  #Get Sorensen and Bray-Curtis distances
  distances=list(dist_soren=vegdist(mat,method="bray",binary=T),
               dist_bray=vegdist(mat,method="bray",binary=F))
  
  #Get unifrac distances:
          distances$unifrac=phyloseq::distance(ps,method="unifrac")
  distances$wunifrac=phyloseq::distance(ps,method="wunifrac",type="samples")
  distances
}

distances=dist(ps_clean_std)

#Function to do PCoA and calculate precentage of variation explained by each axis. 
PCoA_all=function(x){
  pcoa.x=cmdscale(x,eig=T)
  #add percentage of variation explained by each of the first two axis
  pcoa.x$v=c(round(pcoa.x$eig[1]/sum(pcoa.x$eig)*100),
              round(pcoa.x$eig[2]/sum(pcoa.x$eig)*100))
  pcoa.x
}

ord=lapply(distances,PCoA_all)
str(ord)

#Check that order of points is the same for the different measures
head(rownames(ord$dist_soren$points))
head(rownames(ord$unifrac$points))

#Remove from metadata empty samples
metadata_sub=metadata[metadata$SampleID%in%labels(distances$dist_soren),]
metadata_sub=metadata_sub[order(metadata_sub$SampleID),]

#Plot ordinations
shapes=c(6,2,16)
cols=c("#cb181d","#006ba4")
par(oma=c(1, 1, 1, 8), xpd=TRUE)
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))

for (i in 1:4){
    plot(ord[[i]]$points,pch=shapes[metadata_sub$Treatment],col=cols[metadata_sub$Geno],
       xlab=paste("PCoA1(",ord[[i]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[i]]$v[2],"%)",sep=""), main=names(ord)[i])
}

#Add legend to the 4-way plot
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0.58,0.85,legend=c("HETU","VEAL"),pch=16,col=cols,bty ="n",cex=1.2)
legend(0.58,0.55, legend=levels(metadata_sub$Treatment),pch=c(6,2,16),bty ="n",cex=1.2)

dev.off()
par(mar=c(4,4,1,8))
plot(ord[[3]]$points,pch=shapes[metadata_sub$Treatment],col=cols[metadata_sub$Geno],
       xlab=paste("PCoA1(",ord[[3]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[3]]$v[2],"%)",sep=""))
par(xpd=NA)
legend(0.41,0.32,legend=c("HETU","VEAL"),pch=16,col=cols,bty ="n",cex=1.2)
legend(0.41,0.2, legend=levels(metadata_sub$Treatment),pch=c(6,2,16),bty ="n",cex=1.2)

#Permanovas
permanovas=vector("list",4L)
for (i in 1:4){
  permanovas[[i]]=adonis(distances[[i]]~metadata_sub$Geno*metadata_sub$Treatment)
}

```

Richness (given the sampling scheme and thus, not independet from eveness)
```r
presence=function(x){
  x!=0
}
pres_abs=apply(as.matrix(otu_table(ps_clean_std)),c(1,2),presence)
richness=rowSums(pres_abs)
richness_data=data.frame(richness,species=c(rep("HETU",11),rep("VEAL",10)))

rich_summary=ddply(richness_data,.(species),
      summarise,Mean=mean(richness),stdev=sd(richness))
t.test(richness~species,data=richness_data)
```

