---
title: "AlignSangerIllumina"
author: "Maria Rebolleda-Gomez"
date: "2020/08/20"
output: html_document
---

# Align reads and make tree
1. Load packages. 
```r

# Library
library(tidyverse)
library(data.table)
#library(DECIPHER)
library(ape)
library(Biostrings)
library(msa)
library(ggtree)
library(ggExtra)
library(phangorn)
library(phyloseq)
library(castor)
library(picante)
```

2. Load data
```r
ESV16=fread("ESVs.csv")
Isolates=fread("Data_Isolates.csv")
```

3. Format data  for alignment
```r
#Subsample taxa by abundance and only new ASV
ESVsub=ESV16[mean_relab>0.00001,]%>%unique(by = "ASV")

#Open all sequences as biostring object
sequences=DNAStringSet(c(ESVsub$ASV,Isolates$Sequence))

#Trim sequences to similar lenghth
#sequences[width(sequences) >= 1000]=subseq(sequences[width(sequences) >= 1000], start=690, end=1020)

#Label sequences
names(sequences)=c(paste(ESVsub$Genus,1:nrow(ESVsub),sep="_"),
                   paste(Isolates$Species_BLAST,1:nrow(Isolates),sep="_"))

#Two sequences are weird and are messing up the tree
WeirdSequences=c("Pseudomonas caricapapayae_28",'Stenotrophomonas tumulicola_89')
sequences=sequences[!names(sequences)%in%WeirdSequences]
```

4. Align sequences 
```r
#Align sequences with ClustalW default 
alignment=msa(sequences, method="ClustalW", type="dna", order="input")


msaPrettyPrint(alignment, output="pdf", y=c(600,1100),paperWidth=8.5, paperHeight=11, subset=Pseudomonas, alFile ="alignmentPseudo.fasta",file="alignmentPseudo.pdf",showNumbering="none", showLogo="none",consensusColor="ColdHot", showLegend=FALSE, showNames="left",shadingMode="similar")


#Save alignment
phy_mult = as.phyDat(alignment)
write.phyDat(phy_mult,file='~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/Alignment.phy')

DNAStr <- as(alignment, "DNAStringSet")
writeXStringSet(DNAStr, file="~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/alignment_ClustalW.fasta")
```

5. Manually trim alignment with TrimAL. *I tried the automatic methods, but they trimmed too many potentially relevant gaps* 
```bash
trimal -in alignment_ClustalW.fasta -out alignment_trimmed.fasta -htmlout alignment_trimmed.html -select { 0-764,1080-1453 }

```

TrimAL cuts the species names, to make phylogenetic tree I added the names again. 
```r
alignment_trimmed=readDNAStringSet("~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/alignment_trimmed.fasta")
names(alignment_trimmed)=names(sequences)
writeXStringSet(alignment_trimmed, file="~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/alignment_trimmed_newnames.fasta")
```

6. Make phylogenetic tree. Made maximum likelihood tree with tree-IQ. 
```bash
iqtree -s alignment_trimmed_newnames.fasta -nt AUTO -o Chloroflexi_sp._90
```
*tree-IQ* removes repeated sequences. I re-did the tree building with the alignment without repeats. 
```r
iqtree -s alignment_trimmed_newnames.fasta.uniqueseq.phy -nt AUTO -o Chloroflexi_sp._90 -pre uniqueSeqs_out
```


7. Load tree and plot
  7.1. Format metdata for tree. 
```r
Isolates=Isolates[-c(28,89)]

data_tree=data.frame(host_spp=c(ESVsub$spp,Isolates$Host),
                     cult=c(rep(FALSE,nrow(ESVsub)),rep(TRUE,nrow(Isolates))),
                     stringsAsFactors = F)
rownames(data_tree)=gsub(" ", "_", names(sequences))

#Create discrete categories for relative abundance 

breaks=function(x,n){
  l=length(x)
  groupsize=(l/n)%>%round
  interv=seq(1, l,by=groupsize)
  ord=order(x)
  v=vector("numeric", l)
  for (i in 1:n){
    v[ord[interv[i]:interv[i+1]]]=i
  }
  v
}

data_tree$relab=NA
data_tree$relab[!data_tree$cult]=breaks(ESVsub$mean_relab,5)

#Load tree file 
tree=read.tree("~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/Output/uniqueSeqs_out.treefile")

#Remove taxa that is not in the tree (repeated)
data_tree=data_tree[rownames(data_tree)%in%tree$tip.label,]

#We need to add BOTH for the taxa present in both.
#In the ESV data unique ESV selected only from HETU for those that were on both. 
Bothspp=ESV16[ESV16$ASV%in%ESVsub$ASV&spp=="VEAL",]
Bothspp[,sumrel:=sum(Relative_abundance),by="ASV"]
Bothspp=unique(Bothspp,by="ASV")[sumrel!=0,]

vealids=names(sequences)[sequences%in%Bothspp$ASV]
data_tree$host_spp[data_tree$host_spp=="HETU"&rownames(data_tree)%in%vealids]="BOTH"


#For the isolates 
dt=data.table(spp=rownames(data_tree)[data_tree$cult],
              host=data_tree$host_spp[data_tree$cult])
dt[, c("Genus", "sp","id") := tstrsplit(spp, "_", fixed=TRUE)]
dt[,sp2:=paste(Genus, sp, sep=" ")]

dt[,is.HETU:=sp2%in%Isolates[Host=="HETU",Species_BLAST]]
dt[,is.VEAL:=sp2%in%Isolates[Host=="VEAL",Species_BLAST]]

dt[,both:=is.HETU+is.VEAL]
dt$host[dt$both==2]="BOTH"

data_tree$host_spp[rownames(data_tree)%in%dt$spp]=dt$host

#Finally change the name of long species
long=grep("Allorhizobium",rownames(data_tree))
rownames(data_tree)[long]=paste(rep("Rhizobium-Agrobacterium",3),80:82,sep="_")

tree$tip.label[grep("Allorhizobium",tree$tip.label)]=paste(rep("Rhizobium-Agrobacterium",3),80:82,sep="_")


```

  7.2. Plot tree
```r

t=ggtree(tree, layout = 'circular')+
  geom_tiplab2(size=2.6, align=T, color="gray40", linetype="F1")
  #geom_tippoint(aes(color=cult))


heatmap.color=c("gray80","gray60","gray40","gray20","black","#4f372d","#eb6841","#cc2a36","#edc951","#00a0b0")
gheatmap(t,data_tree, width = 0.2, offset=1.8, 
         colnames=F)+
  scale_fill_manual(values=heatmap.color)

```

Remove gaps but mantained trimmed sequneces
```bash
sed 's/[-]//g' alignment_trimmed_newnames.fasta > alignment_trimmed_newnames_nogap.fasta
```

# Pirwise alignment
Open and perform pairwise alignment of pairwise sequences
```r
trimmed=readDNAStringSet("alignment_trimmed_newnames_nogap.fasta")

## R function for pairwise alignment, getting percent identity and number of mistmatch
align_to_all <- function(seq1, seq2, type = "global"){
  # Make pairwise alignment
  algn <- Biostrings::pairwiseAlignment(seq1, seq2, type = type)
  
  #Get percent identity
  id <- pid(algn)
  
  #Get the number of mistmatches
  n_mis <- nmismatch(algn)
  
  #Get alignment lenghth and compare to small sequence 
  al_len <- nchar(algn)
  l1 <- nchar(seq1)
  l2 <- nchar(seq2)
  dif_len <- min(l1,l2)%>%-al_len
  
  c(identity=id, nmism=n_mis, alignment_length=al_len, dif_lengths=dif_len)
}


#For each ESV sequence, find the closest isolate and save percentage similarity
l_seqs1 <- nrow(Isolates)-1
l_seqs2 <- nrow(ESVsub)

tmp <- matrix(NA, ncol=4, nrow=l_seqs1)
similarity <- matrix(NA, ncol=5, nrow=l_seqs2)

for (j in 1:l_seqs2){
  for (i in (l_seqs2+1):(l_seqs2+l_seqs1)){
  tmp[i-l_seqs2,] <- align_to_all(trimmed[j],trimmed[i])
  }
  
#Get data from the maximum similarity alignment
index <- which(tmp[,1]==max(tmp[,1],na.rm=T))[1]
similarity[j,] <- c(tmp[index,],index)

}

similarity=as.data.table(similarity)
names(similarity)=c("identity","nmism","alignment_length","delta_length","index_no")
similarity$mean_relab=ESVsub$mean_relab
similarity$host_spp=data_tree$host_spp[1:l_seqs2]

fwrite(similarity,"../Output/SequenceSimilarity.csv")

cols2=c("#4f372d","#cc2a36","#00a0b0")

median=median(similarity$identity)
# 97.66276
mean=median(similarity$identity)
#95.77618

p=ggplot(similarity, aes(y=mean_relab,x=identity,colour=host_spp))+
  geom_point()+
  scale_color_manual(values=cols2)+
  xlab("% sequence identity")+
  ylab("Mean relative abundance")+
  geom_vline(xintercept=97, linetype=2, color="gray60")+
  geom_vline(xintercept=median, linetype=2, color="red")+
  theme_bw()
  
  
ggMarginal(p, type = "histogram", margins="x", groupFill = F)


t.test(similarity$mean_relab[similarity$identity>median],
       similarity$mean_relab[similarity$identity<median])



```

# Unifrac calculation and simulations

### Observed Unifrac
Format data and calculate unifrac between cultivated and sequencing
```r
#Tree loads as unrooted, so for unifrac I re-assigned the root
tree=root(tree, outgroup = "Chloroflexi_sp._90", resolve.root = T)

#Create presence absence table with all
pres_abs=data.frame(cultivables=as.numeric(data_tree$cult),
                    sequencing=as.numeric(!data_tree$cult))[-183,]

rownames(pres_abs)=rownames(data_tree)[-183]

ps=phyloseq(otu_table(pres_abs,taxa_are_rows=T),phy_tree(tree))
unifrac=phyloseq::distance(ps,method="unifrac")

```
### Null-models
Case 1. Null expectation mantaining the same sample sizes

```r
#3. Take random samples - randomize either by species or by cult-no cult 
random=vector("numeric", 1000)

for (i in 1:1000){
  tmp_otu=pres_abs
  rownames(tmp_otu)=sample(rownames(pres_abs))
  ps_tmp=phyloseq(otu_table(tmp_otu,taxa_are_rows=T),phy_tree(tree))
  random[i]=phyloseq::distance(ps_tmp,method="unifrac")[1]
}


ggplot(data=as.data.frame(random), aes(x=random))+
  geom_histogram(bins = 30, aes(y=(..count..)/sum(..count..)))+
  geom_vline(xintercept = unifrac[1], linetype=2, color="red")+
  xlab("Unifrac distance")+
  ylab("Frequency")+
  theme_bw()
  
#Rough frequentist version of p-value. 
#Given large sample sizes it should be eq. to geting the are under curve after onserved
pvalue=sum(random>=unifrac[1])/1000

```


Case 2. Null expectation with even sample sizes
```r
# Create function to sample random presence absence
sample_rdm_rows=function(all, n_samp, rdm=T, samp=NULL , prob=NULL){
  if (rdm){
    rows_pres <- sample(1:length(all), size=n_samp, prob=prob)
  } else {
    rows_pres <- sample(samp, size=n_samp, prob=prob)
  }
  y <- all
  y[rows_pres] <- 1
  y
}

#Create matrix empty with all 0s
nrows <- nrow(pres_abs)
empty <- matrix(0L, ncol=100, nrow=nrows)

# Apply sampling function to each column. 
# With this arguments it samples from all the branches with equal prob, 
# and independently of their cultivation status. 
null <- apply(empty, 2, sample_rdm_rows, n_samp=sum(pres_abs$cultivables))
rownames(null) <- rownames(pres_abs)

#Get unifrac
ps_null <- phyloseq(otu_table(null, taxa_are_rows=T), phy_tree(tree))
unifrac_null <- phyloseq::distance(ps_null,method="unifrac") %>%
  as.vector

ggplot(data=as.data.frame(unifrac_null), aes(x=unifrac_null))+
  geom_histogram()+
  geom_vline(xintercept = unifrac[1], linetype=2, color="red")+
  xlab("Unifrac distance")+
  theme_bw()

pvalue2=sum(unifrac_null>=unifrac[1])/length(unifrac_null)
  
```

Case 3. Predicted from sampling same number of ESVs as there are isolates (no abundance)
```r
empty <- matrix(0L, ncol=1000, nrow=nrows)

# With this arguments it samples from only ESVs with equal prob, 
esv_eq <- apply(empty, 2, sample_rdm_rows, n_samp=sum(pres_abs$cultivables), rdm=F, samp=which(pres_abs$sequencing==1))

#Calculate unifrac for each sample with cultivables
unifrac_cult=function(x,cult){
  pres_abs_tmp <- c(cult, x)%>%
    matrix(ncol=2)
  rownames(pres_abs_tmp) <- rownames(pres_abs)
  ps <- phyloseq(otu_table(pres_abs_tmp, taxa_are_rows=T), phy_tree(tree))
  uni <- phyloseq::distance(ps,method="unifrac")
  uni[1]
}

unifrac_samesize=apply(esv_eq, 2, unifrac_cult, cult=pres_abs$cultivables)

uni_data=data.frame(unifrac=c(sample(unifrac_null,1000),unifrac_samesize),
                    model=rep(c("Null","Equal Size"),each=1000))

ggplot(data=uni_data, aes(x=unifrac, fill=model))+
  geom_density(alpha=0.7,aes(y=(..count..)/sum(..count..)))+
  geom_vline(xintercept = unifrac[1], linetype=2, color="red")+
  scale_fill_manual(values=c("#d7191c","#bababa"))+
  xlab("Unifrac distance")+
  ylab("frequency")+
  theme_bw()
  
```

Case 4. Predicted from sampling same number of ESVs as there are isolates (no abundance)
```r
#With these arguments it samples from only ESVs according to thier relative abundances
esv_ab <- apply(empty, 2, sample_rdm_rows, n_samp=sum(pres_abs$cultivables), rdm=F, samp=which(pres_abs$sequencing==1), prob=ESVsub$mean_relab)

unifrac_relab=apply(esv_ab, 2, unifrac_cult, cult=pres_abs$cultivables)

uni_data=data.frame(unifrac=c(sample(unifrac_null,1000),unifrac_samesize,unifrac_relab),
                    model=rep(c("Null","Equal Size", "+rel. abundance"),each=1000))

ggplot(data=uni_data, aes(x=unifrac, fill=model))+
  geom_vline(xintercept = quantile(unifrac_null, probs=c(0.025,0.975)),  color="gray20")+
  geom_density(alpha=0.7)+
  geom_vline(xintercept = unifrac[1], linetype=2, color="red")+
  scale_fill_manual(values=c("#d7191c","#bababa","#fdae61"))+
  xlab("Unifrac distance")+
  theme_bw()



```

# Calculate phylogenetic distance

```r
#Calculate expected phylogenetic distance and variance
E_pd=expected.pd(tree)
V_pd=variance.pd(tree,upper.bound = F)

#Save data
setwd("~/Dropbox/Projects/FloweMicrobeProject/Flower-microbe Color Curiculum Dev/Analyses/AmpliconSequencing/")
write.csv(E_pd, "Output/expected_pd.csv")
write.csv(V_pd, "Output/variance_expectedpd.csv")


PD=pd(t(pres_abs),tree,include.root = F)

ex=E_pd$expected.pd[PD$SR]
sd=sqrt(V_pd$variance.pd[PD$SR])
Z=(PD$PD-ex)/sd
Z

pd_data=data_frame(Z=Z, samples=c("Cultivated","Amplicon sequencing"))

ggplot(pd_data,aes(x=samples,y=Z))+
  geom_rect(xmin=0,xmax=Inf,ymin=-1.96,ymax=1.96, alpha=0.05)+
  geom_hline(aes(yintercept=-1.96),linetype="dashed")+
  geom_hline(aes(yintercept=1.96),linetype="dashed")+
  geom_hline(aes(yintercept=0))+
  geom_point(size=3,shape=21, aes(fill=samples))+
  ylab("Phylogenetic distance (z score)")+
  scale_fill_manual(values=c("#eb6841","#edc951"))+
  theme_bw()

```





