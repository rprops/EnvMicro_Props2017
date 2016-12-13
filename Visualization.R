library("flowViz")
library("flowCore")
library("Phenoflow")
library("gridExtra")
library("ggplot2")
library("Phenoflow")

my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col="black", fill="lightblue", alpha=0.2,border=NA,lwd=2),
  panel.background=list(col="lightgray"),
  background=list(col="white"))


path = "data"
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

flowData_transformed <- flowData_transformed[which(flowCore::sampleNames(flowData_transformed) 
                                                   %in% c("Q3_T_t0_rep1.fcs","Q3_T_t1.5_rep1.fcs","Q3_T_t3_rep1.fcs",
                                                         "C1_t0_rep1.fcs","C1_t1.5_rep1.fcs","C1_t3_rep1.fcs"))]

MyText <- c("Control - 0h", "Control - 1.5h", "Control - 3h",
            "Feeding - 0h", "Feeding - 1.5h", "Feeding - 3h")

# colramp <- colorRampPalette(IDPcolorRamp(21))
require('RColorBrewer')
require("IDPmisc")
png(file="Fig1_scatters_200_log.png", units="in", width=12, height=10,  pointsize=8,   res=500)
print(xyplot(`FL3-H`~`FL1-H`, data=flowData_transformed,index.cond=list(c(1:6)),
             xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
             smooth=TRUE, xlim=c(0.5,1),ylim=c(0.1,1),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=1, by=.1),cex=1),
                         y=list(at=seq(from=0, to=1, by=.2),cex=1)), layout=c(3,2),
             strip=strip.custom(factor.levels=MyText),
             margin=TRUE,
             binTrans="log"
      )
)
dev.off()

comp_t0 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t0.fcs" | pos== "C2_t0.fcs"),
                        comp1=(pos=="Q1_T_t0.fcs" | pos== "Q2_T_t0.fcs" | pos== "Q3_T_t0.fcs"), 
                        thresh=0.04)
comp_t1.5 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t1.5.fcs" | pos== "C2_t1.5.fcs"),
                          comp1=(pos=="Q1_T_t1.5.fcs" | pos== "Q2_T_t1.5.fcs" | pos== "Q3_T_t1.5.fcs"), 
                          thresh=0.04)
comp_t3 <- fp_contrasts(x=fbasis, comp2=(pos=="C1_t3.fcs" | pos== "C2_t3.fcs"),
                        comp1=(pos=="Q1_T_t3.fcs" | pos== "Q2_T_t3.fcs" | pos== "Q3_T_t3.fcs"), 
                        thresh=0.04)


### Merge data frames

comp_total <- rbind(comp_t0, comp_t1.5, comp_t3)
comp_total <- data.frame(comp_total, Timepoint=as.factor(c(rep("0h Feeding", nrow(comp_t0)), 
                                                 rep("1.5h Feeding", nrow(comp_t1.5)),
                                                 rep("3h Feeding", nrow(comp_t3)))))

vtot <- ggplot(comp_total, aes(`FL1.H`, `FL3.H`, z = Density))+
  geom_tile(aes(fill=Density)) + 
  geom_point(colour="gray", alpha=0.4)+
  scale_fill_distiller(palette="RdBu", na.value="white") + 
  stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  geom_contour(color = "white", alpha = 1)+
  facet_grid(~Timepoint)+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)")+
  theme(axis.title=element_text(size=16,face="bold"), strip.text.x=element_text(size=16,face="bold"),
        legend.title=element_text(size=12),legend.text=element_text(size=12),
        axis.text = element_text(size=12))
  
vtot

png(file = "contrasts_0.04.png", width = 12, height = 6, res = 500, 
    units = "in", pointsize = 10)
print(vtot)
dev.off()