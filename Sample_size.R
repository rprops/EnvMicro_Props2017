library("Phenoflow")
library("gridExtra")
library("grid")

### Load data
path = "data_mussel"
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

### Preprocess data according to standard protocol
flowData_transformed <- transform(flowData, `FL1-H` = asinh(`FL1-H`), `SSC-H` = asinh(`SSC-H`), 
                                  `FL3-H` = asinh(`FL3-H`), `FSC-H` = asinh(`FSC-H`))
param = c("FL1-H", "FL3-H", "SSC-H", "FSC-H")
flowData_transformed = flowData_transformed[, param]
remove(flowData)

### Create a PolygonGate for denoising the dataset Define coordinates for
### gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.5,8.5,15,15,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H", "FL3-H")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Total Cells")

### Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data = flowData_transformed[1], filter = polyGate1, 
       scales = list(y = list(limits = c(0, 14)), x = list(limits = c(6, 16))), 
       axis = axis.default, nbin = 125, par.strip.text = list(col = "white", 
                                                              font = 2, cex = 2), smooth = FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 
                                                                     2, max), use.exprs = TRUE)
max = max(summary[, 1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed, `FL1-H` = mytrans(`FL1-H`), 
                                  `FL3-H` = mytrans(`FL3-H`), `SSC-H` = mytrans(`SSC-H`), `FSC-H` = mytrans(`FSC-H`))

### Subsample at various depths and calculate diversity metrics with 100
### bootstraps Notice: this will use some CPU/RAM
for (i in c(10, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 2000, 2500, 
            3000, 5000, 10000, 15000, 20000, 30000)) {
  for (j in 1:100) {
    fs1 <- FCS_resample(flowData_transformed[10], replace = TRUE, sample = i)
    fp <- flowBasis(fs1, param, nbin = 128, bw = 0.01, normalize = function(x) x)
    div.tmp <- Diversity(fp, d = 3, R = 100)
    div.tmp <- cbind(div.tmp, size = i)
    if (j == 1) 
      results <- div.tmp else results <- rbind(results, div.tmp)
  }
  if (i == 10) 
    results.tot <- results else results.tot <- rbind(results.tot, results)
}

### Create plots
D0 <- ggplot(data = results.tot, aes(x = factor(size), y = D0)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + labs(x = "Sample size (nr. of cells)", 
                                                                            y = expression('Phenotypic diversity - D'[0])) + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                                                             hjust = 1), 
                                                                                                                                                axis.text=element_text(size=11),
                                                                                                                                                axis.title=element_text(size=15))
D1 <- ggplot(data = results.tot, aes(x = factor(size), y = D1)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + labs(x = "Sample size (nr. of cells)", 
                                                                            y = expression('Phenotypic diversity - D'[1])) + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                                                             hjust = 1), 
                                                                                                                                                axis.text=element_text(size=11),
                                                                                                                                                axis.title=element_text(size=15))

D2 <- ggplot(data = results.tot, aes(x = factor(size), y = D2)) + # geom_jitter(alpha=0.7, size=1)+
  geom_boxplot(alpha = 0.2, color = "blue", fill = "blue", size = 1) + labs(x = "Sample size (nr. of cells)", 
                                                                            y = expression('Phenotypic diversity - D'[2])) + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                                                             hjust = 1), 
                                                                                                                                                axis.text=element_text(size=11),
                                                                                                                                                axis.title=element_text(size=15))
# Sample used for this assessment
flowCore::sampleNames(flowData_transformed[10])

png(file = "sample_size_effect.png", width = 12, height = 6, res = 500, 
    units = "in", pointsize = 10)
grid.arrange(D0, D1, D2, ncol = 3, top = textGrob("", 
                                                  gp = gpar(fontsize = 20, font = 3)))
dev.off()
