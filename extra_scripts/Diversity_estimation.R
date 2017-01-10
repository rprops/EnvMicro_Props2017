require('qdap')
require('data.table')
require("phyloseq")
require("dplyr")
# require("Phenoflow")
library("iNEXT")
library("ggplot2")
### Import MED table

#div.otu <- Diversity_16S(physeq.otu, brea=FALSE, R=100, thresh=500)

### Import otu data
physeq.otu <- import_mothur(mothur_shared_file ="total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared" ,
                            mothur_constaxonomy_file = "total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")
otu_table(physeq.otu) <- t(otu_table(physeq.otu))

physeq.otu <- prune_samples(sample_sums(physeq.otu)>10000, physeq.otu)

df <- data.frame(otu_table(physeq.otu))

D <- estimateD(df, datatype = "abundance", base = "size", level = NULL, conf = 0.95)
