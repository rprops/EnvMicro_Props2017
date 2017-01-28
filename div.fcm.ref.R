library("Phenoflow")
library("dplyr")
library("phyloseq")
library("data.table")
library("qdap")
library("vegan")
library("ggplot2")

### Calculate 16S diversity for reference samples
phy.ref <- import_mothur(mothur_shared_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list_CYANO.shared")
tax.ref <- data.frame(fread("16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons_CYANO.taxonomy"))
tax.ref <- genX(tax.ref, "(", ")"); tax.ref  <- data.frame(tax.ref,stringsAsFactors =FALSE)
tax.ref <- tax_table(as.matrix(tax.ref)); rownames(tax.ref) <- tax.ref[,1]
phy.ref <- phyloseq(phy.ref, tax.ref)

### Only consider samples with > 10000 reads
phy.ref <- prune_samples(sample_sums(phy.ref)>10000, phy.ref)


# div.ref <- Diversity_16S(phy.ref, R=100, brea=FALSE, thresh=500)
div.ref <- data.frame(div.ref)


### Dissimilarity matrix for beta diversity
phy.ref <- prune_taxa(taxa_sums(phy.ref) > 0, phy.ref)
dist.ref <- vegdist(t(otu_table(phy.ref)), upper=TRUE)
dist.ref <- as.matrix(dist.ref)
cols <- c()
for(i in 1:ncol(dist.ref)) cols <- c(cols, rep(colnames(dist.ref)[i], nrow(dist.ref)))
dist.ref <- data.frame(Sample1 = row.names(dist.ref), Sample2 = cols, distance = c(dist.ref))
dist.ref <- data.frame(Sample_seq_comb = paste(dist.ref[,1], dist.ref[,2], sep="-"), distance.seq = dist.ref[,3])
  
### Export data
write.csv2(dist.ref, "dist.16S.ref.csv")
write.csv2(div.ref, "otu.diversity16S.ref.csv")

### Calculate FCM diversity for reference samples

path = "./data_reference/FCM_CW"

flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

### Randomly resample to the lowest sample size
# flowData_transformed <- FCS_resample(flowData_transformed, replace=FALSE)

### Calculate average fingerprint over replicates with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128,
                    bw=0.01,normalize=function(x) x)
INDICES <- gsub(flowCore::sampleNames(flowData_transformed), pattern="_rep.*", replacement="")
fbasis <- by(fbasis@basis, INDICES = INDICES, FUN = colMeans)
fbasis <- do.call(rbind, fbasis)
fbasis <- round(fbasis, 3)

### Calculate pair-wise distances
dist.fbasis <- vegdist(fbasis, upper=TRUE)
dist.fbasis <- as.matrix(dist.fbasis)
cols <- c()
for(i in 1:ncol(dist.fbasis)) cols <- c(cols, rep(colnames(dist.fbasis)[i], nrow(dist.fbasis)))
dist.fbasis <- data.frame(Sample1 = row.names(dist.fbasis), Sample2 = cols, distance.fcm = c(dist.fbasis))

### Calculate ecological parameters from normalized fingerprint 
### Densities will be normalized to the interval [0,1]
### n = number of replicates
### d = rounding factor
# Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)
Diversity.ref <- Diversity_rf(flowData_transformed, d=3, param=param, R=3)

### Rename for further processing
div.fcm.ref <- Diversity.ref 

### Now select the ones for which parallel sequencing data is available
meta.seq <- read.csv2("metadata_ref.csv")
meta.seq <- meta.seq[meta.seq$Sequenced=="yes",]

div.fcm.ref <- div.fcm.ref[grep(pattern = paste(meta.seq$Sample.name,collapse="|"),  x = as.character(div.fcm.ref$Sample_names)),]
dist.fbasis <- dist.fbasis[grep(pattern = paste(meta.seq$Sample.name,collapse="|"),  x = as.character(dist.fbasis$Sample1)),]
dist.fbasis <- droplevels(dist.fbasis)

### Now calculate the averages of these technical replicates
groupLevels <- factor(gsub(div.fcm.ref$Sample_names, pattern="_rep.*", replacement="", fixed=FALSE))
means <- do.call(rbind,by(div.fcm.ref[,2:4], INDICES = groupLevels, 
                FUN = colMeans))
errors <- do.call(rbind,by(div.fcm.ref[,5:7], INDICES = groupLevels, 
             FUN = function(x) sqrt(colSums(x^2))/ncol(x)))
div.fcm.ref.merged <- data.frame(cbind(means,errors), sample_fcm = rownames(means))

### And match them to the corresponding sequencing data
div.seq.ref <- read.csv2("otu.diversity16S.ref.csv")
lb <- read.csv2("labels_ref.csv")
div.fcm.ref.merged <- inner_join(div.fcm.ref.merged, lb, by=c("sample_fcm"="Sample_fcm"))
dist.fbasis <- inner_join(dist.fbasis, lb, by=c("Sample1"="Sample_fcm"))
dist.fbasis <- inner_join(dist.fbasis, lb, by=c("Sample2"="Sample_fcm"))
dist.fbasis <- data.frame(dist.fbasis[,1:3], Sample_seq_comb = paste(dist.fbasis[,4], dist.fbasis[,5], sep="-"))


  
### Adjust colnames
colnames(div.fcm.ref.merged) <- c("D0.fcm","D1.fcm","D2.fcm","sd.D0.fcm","sd.D1.fcm","sd.D2.fcm","Sample_fcm","Sample_seq")

### Merge
div.total.ref <- inner_join(div.fcm.ref.merged, div.seq.ref, by=c("Sample_seq"="Sample"))
dist.total.ref <- inner_join(dist.fbasis, dist.ref, by=c("Sample_seq_comb"="Sample_seq_comb"))
dist.total.ref <- dist.total.ref[dist.total.ref$distance.fcm != 0 | dist.total.ref$distance.seq != 0,]

### Tidy up a bit

div.total.ref <- data.frame(div.total.ref[,7:8], div.total.ref[,1:6], div.total.ref[,9:18])

### Write csv

# write.csv2(div.total.ref, "div.ref.merged.csv")

ggplot(data=dist.total.ref, aes(x=distance.fcm, y=distance.seq))+
  geom_point()+
  geom_smooth()

