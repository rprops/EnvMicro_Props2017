library('qdap')
library('data.table')
library("phyloseq")
library("dplyr")
library("Phenoflow")
library("ggplot2")
library("gridExtra")
source("functions.R")

### Import otu data
physeq.otu <- import_mothur(mothur_shared_file ="16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared" ,
                            mothur_constaxonomy_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")
otu_table(physeq.otu) <- t(otu_table(physeq.otu))
physeq.otu <- prune_samples(sample_sums(physeq.otu)>10000, physeq.otu)
physeq.otu <- prune_taxa(taxa_sums(physeq.otu)>0, physeq.otu)

### Select samples for which you have fcm data
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu), pattern=".renamed", replacement="")
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu), pattern=".renamed", replacement="")
meta.seq <- read.csv2("metadata.csv")

### make sure to have run Comparison.R before for data.total.final
physeq.otu <- prune_samples(sample_names(physeq.otu) %in% data.total.final$Sample, physeq.otu)

### Annotate with metadata
rownames(data.total.final) <- data.total.final$Sample
sample_data(physeq.otu) <- data.total.final

### Rescale
physeq.otu <- scale_reads(physeq.otu)
physeq.otu <- transform_sample_counts(physeq.otu, function(x) x/sum(x))

### Run beta diversity analysis on 16s data
pcoa <- ordinate(
  physeq = physeq.otu, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

pcoa.df <- data.frame(pcoa$vectors, sample_data(physeq.otu))
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

beta.pcoa <- ggplot(data=pcoa.df, aes(x=Axis.1, y=Axis.2, colour=Season, shape=Lake))+
  geom_point(alpha=0.7, size=7, aes(fill=Season))+
  scale_shape_manual(values=c(21,24))+
  # scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  # guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  scale_colour_brewer(palette="Accent")+
  scale_fill_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"), fill="Season", shape="Environment",
       title="Taxonomic beta diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))
print(beta.pcoa)

### Start beta diversity analysis on FCM data
path = "data_reference/FCM_MI"
flowData_transformed  <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData_transformed,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
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

### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

### Run beta diversity
pos <- gsub(rownames(fbasis@basis),pattern="_rep.*", replacement="")
beta.div <- beta_div_fcm(fbasis, INDICES=pos, ord.type="PCoA")
df.beta.fcm <- data.frame(Sample = rownames(beta.div$points), beta.div$points)
df.beta.fcm$Sample <- gsub(df.beta.fcm$Sample, pattern="MI5", replacement = "M15")
df.beta.fcm <- inner_join(df.beta.fcm, data.total.final, by=c("Sample"="Sample_fcm"))
var2 <- round(beta.div$eig/sum(beta.div$eig)*100, 1)
df.beta.fcm <- droplevels(df.beta.fcm)


### Plot beta diversity
beta.pcoa.fcm <- ggplot(data=df.beta.fcm, aes(x=X1, y=-X2, shape=Lake, colour=Season))+
  scale_shape_manual(values=c(21,24))+
  geom_point(size=7)+
  geom_point(size=7, aes(fill = Season), alpha=0.7)+
  # scale_size(range=c(4,10), breaks=c(0,0.5,1,1.5,2,2.5,3))+ 
  # guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_bw()+
  scale_fill_brewer(palette="Accent")+
  scale_colour_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var2[1], "%)"), y = paste0("Axis2 (",var2[2], "%)"), fill="Season", shape = "Environment",
       title="Phenotypic beta diversity")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))

print(beta.pcoa.fcm)
grid.arrange(beta.pcoa, beta.pcoa.fcm, ncol=2)


### All together
png("beta-div_fcm_seq.png",width=1.6*7*1.65,height=5*1.5,res=500,units="in")
grid.arrange(beta.pcoa, beta.pcoa.fcm, ncol=2)
dev.off()

### Procrustes analysis
fbasis1 <- fbasis
fbasis1@basis <- fbasis1@basis/apply(fbasis1@basis, 1, max)
fbasis1@basis <- round(fbasis1@basis, 3)
x <- by(fbasis1@basis, INDICES = pos, FUN = colMeans)
x <- do.call(rbind, x)
rownames(x) <- gsub(rownames(x), pattern = "MI5", replacement = "M15")
x <- x[(rownames(x) %in% data.total.final$Sample_fcm), ]

#Rename and order rows of fcm/seq data
tmp <- data.frame(sample_fcm=rownames(x))
tmp <- left_join(tmp, data.total.final, by=c("sample_fcm"="Sample_fcm"))
rownames(x) <- tmp$Sample; remove(tmp)
x <- x[order(rownames(x)),]
x.data <- data.total.final[data.total.final$Sample %in% rownames(x),]
x.data <- x.data[order(x.data$Sample),]
otu_table(physeq.otu) <- otu_table(physeq.otu)[order(rownames(otu_table(physeq.otu))),]

# Run PcoA
dist.fcm <- vegan::vegdist(x)
pcoa.fcm <- cmdscale(dist.fcm)
dist.seq <- vegan::vegdist(otu_table(physeq.otu))
pcoa.seq <- cmdscale(dist.seq)

# Run procrustes + permutation
# dist.proc <- vegan::procrustes(dist.seq, dist.fcm)
dist.prot <- vegan::protest(pcoa.seq,pcoa.fcm)
plot(dist.prot)
summary(dist.prot)

# Run PERMANOVA to evaluate if conclusions are similar between FCM/seq
disper.fcm <- betadisper(dist.fcm, group=x.data$Season)
print(disper.fcm)
plot(disper.fcm)
disper.seq <- betadisper(dist.seq, group=x.data$Season)
print(disper.seq)
plot(disper.seq)

# Permutations are constrained within each sampling year 
perm <- how(nperm = 999)
setBlocks(perm) <- with(x.data, Year)
permanova.fcm <- adonis(dist.fcm~Season*Lake, data=x.data, permutations = perm)
permanova.seq <- adonis(dist.seq~Season*Lake, data=data.frame(sample_data(physeq.otu)), permutations = perm)

