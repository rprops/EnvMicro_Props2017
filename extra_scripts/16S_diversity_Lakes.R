library("data.table")
library("Phenoflow")
library('qdap')
library("phyloseq")
source("functions.R")
### Import data
phy.otu <- import_mothur(mothur_shared_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared",
                         mothur_constaxonomy_file = "16S/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")

### Only look at samples with > 10000 reads
phy.otu <- prune_samples(sample_sums(phy.otu)>10000, phy.otu)

### Replace .renamed for MI13 data
sample_names(phy.otu) <- gsub(sample_names(phy.otu), pattern=".renamed",replacement="",fixed=TRUE)

### Replace "-" from phy.otu sample names
sample_names(phy.otu) <- gsub(sample_names(phy.otu), pattern="-", replacement="")

### Select only the useful samples (F/W/P)
sl <- read.csv2("Lakes_metadata")
sl_F <- sl$Sample_16S[sl$Lake!="Inland" & sl$Platform=="Accuri"]
sl_W <- gsub(sl_F, pattern="F", replace="W")
sl_W[grep(sl_W, pattern="*.W.*")] <- paste(sl_W[grep(sl_W, pattern="*.W.*")],"D",sep="")
sl_W <- sl_W[-c(1:69)]
sl_P <- gsub(sl_F, pattern="F", replace="P")
sl_P <- sl_P[-c(1:69)]

# Free living bacteria
phy.otu.F <- prune_samples(samples=sample_names(phy.otu) %in% sl_F, x=phy.otu)
# Particle associated bacteria
phy.otu.P <- prune_samples(samples=sample_names(phy.otu) %in% sl_P, x=phy.otu)
# Whole fraction
phy.otu.W <- prune_samples(samples=sample_names(phy.otu) %in% sl_W, x=phy.otu)

### Rescale if desired to 10000
# phy.otu.F <- scale_reads(phy.otu.F, n = 10000)
# phy.otu.P <- scale_reads(phy.otu.P, n = 10000)
# phy.otu.W <- scale_reads(phy.otu.W, n = 10000)

### Calculate diversities
div.F <- Diversity_16S(phy.otu.F, R=100, brea=FALSE, thresh=500)
div.F <- data.frame(Sample = rownames(div.F), div.F)
div.W <- Diversity_16S(phy.otu.W, R=100, brea=FALSE, thresh=500)
div.W <- data.frame(div.W)
div.P <- Diversity_16S(phy.otu.P, R=100, brea=FALSE, thresh=500)
div.P <- data.frame(div.P)

### Import FCM ref
write.csv2(div.F, "Lakes_diversity16S_F.csv", row.names = FALSE)
