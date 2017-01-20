library("flowViz")
library("flowCore")
library("Phenoflow")
library("gridExtra")
library("ggplot2")
library("Phenoflow")
library("egg")
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
# flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
max <- 14.99341
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

# Make xyplot of gating strategy (6 samples)
require('RColorBrewer')
png(file="FigS1_scatters_200_log.png", units="in", width=12, height=10,  pointsize=8,   res=500)
print(xyplot(`FL3-H`~`FL1-H`, data=flowData_transformed,index.cond=list(c(1:6)),
             filter=flist,
             xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
             smooth=FALSE, xlim=c(0.5,1),ylim=c(0.1,1),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=1, by=.1),cex=1),
                         y=list(at=seq(from=0, to=1, by=.2),cex=1)), layout=c(3,2),
             strip=strip.custom(factor.levels=MyText),
             margin=TRUE,
             binTrans="log"
      )
)
dev.off()

# Make xyplot of gating strategy (1 sample)
png(file="Fig1_scatters_400_log.png", units="in", width=8, height=8,  pointsize=8,   res=500)
print(xyplot(`FL3-H`~`FL1-H`, data=flowData_transformed,index.cond=list(c(1)),
             filter=flist[1],
             xbins=400,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
             smooth=FALSE, xlim=c(0.5,1),ylim=c(0.1,1),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=1, by=.1),cex=1),
                         y=list(at=seq(from=0, to=1, by=.2),cex=1)), layout=c(1,1),
             margin=TRUE,
             binTrans="log",
             strip=strip.custom(factor.levels="")
)
)
dev.off()

### Make same biplots but with LNA/HNA filters
my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col=c("black"),fill="lightblue", alpha=0.5,lwd=4),
  panel.background=list(col="lightgray"),
  background=list(col="white"))

sqrcut1 <- matrix(c(asinh(12500),asinh(12500),15,15,3,9.55,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(12500),asinh(12500),3,8,9.55,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")

filters <- filters(list(rGate_LNA, rGate_HNA))
filter1 <- flowCore::filter(filter=rGate_LNA, flowData_transformed)
filter2 <- flowCore::filter(filter=rGate_HNA, flowData_transformed)
flist2 <- list(filter1,filter2)
names(flist2) <- c("LNA","HNA")
flist <- list(filters , filters, filters, filters, filters, filters)
names(flist) <- flowCore::sampleNames(flowData_transformed)


png(file="Fig1B_HNA.LNA_200_log.png", units="in", width=12, height=10,  pointsize=8,   res=500)
print(xyplot(`FL3-H`~`FL1-H`, data=flowData_transformed,index.cond=list(c(1:6)),
             filter=flist,
             xbins=200,nbin=128, par.strip.text=list(col="black", font=3,cex=1.85), 
             smooth=FALSE, xlim=c(0.5,1),ylim=c(0.1,1),xlab=list(label="Green fluorescence intensity (FL1-H)",cex=2),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=2),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=1, by=.1),cex=1),
                         y=list(at=seq(from=0, to=1, by=.2),cex=1)), layout=c(3,2),
             strip=strip.custom(factor.levels=MyText),
             margin=FALSE,
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

### For HNA dynamics (from diversity_analysis.R script)
p12 <- ggplot(data=results, aes(x=factor(Time), y=HNA.cells, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=5,alpha=0.9)+
  scale_fill_manual(values=myColours[c(1,2)])+
  theme_bw()+
  labs(y=expression("HNA cells µL"^"-1"), x="Time (h)", title="C", fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        legend.direction = "horizontal",legend.position = "bottom"
        )+ 
  # guides(fill=FALSE)+
  geom_errorbar(aes(ymin=HNA.cells-sd.HNA.cells, ymax=HNA.cells+sd.HNA.cells), width=0.075)+
  ylim(200,525)

p13 <- ggplot(data=results, aes(x=factor(Time), y=LNA.cells, fill=Treatment)) + 
  # geom_boxplot(alpha=0.9)+
  geom_point(shape=21, size=5,alpha=0.9, position=position_dodge(width=0.3))+
  scale_fill_manual(values=myColours[c(1,2)])+
  # geom_smooth(formula=y ~ x, color="black")+
  # geom_boxplot(mapping=factor(Time),alpha=0.4,outlier.shape=NA)+
  theme_bw()+
  labs(y=expression("LNA cells µL"^"-1"), x="Time (h)", title="B", fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.direction = "horizontal",legend.position = "bottom")+ 
  # guides(fill=FALSE)+
  geom_errorbar(aes(ymin=LNA.cells-sd.LNA.cells, ymax=LNA.cells+sd.LNA.cells), width=0.075
                , position=position_dodge(width=0.3))+
  ylim(1000,1325)



vtot <- ggplot(comp_total, aes(`FL1.H`, `FL3.H`, z = Density))+
  geom_tile(aes(fill=Density)) + 
  geom_point(colour="gray", alpha=0.4)+
  scale_fill_distiller(palette="RdBu", na.value="white") + 
  stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  geom_contour(color = "white", alpha = 1)+
  facet_grid(~Timepoint)+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)",title="A")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )
  
vtot

png(file = "contrasts_0.04.png", width = 12, height = 6, res = 500, 
    units = "in", pointsize = 10)
print(vtot)
dev.off()

png(file = "contrasts_0.04_combined_cells_grid_C_dodge.png", width = 12, height = 12, res = 500, 
    units = "in", pointsize = 10)
# grid.arrange(vtot,p12,p13, nrow=3)
g1 <- ggplotGrob(p13);
g2 <- ggplotGrob(p12b);
gv <- ggplotGrob(vtot);
fg1 <- gtable_frame(g1)
fg2 <- gtable_frame(g2)
fg12 <- gtable_frame(cbind(fg1,fg2))
fg3 <- gtable_frame(gv)
grid.newpage()
combined <- rbind(fg3, fg12)
grid.draw(combined)
dev.off()

### Average HNA removal
t0 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"])
t3 <- mean(results$HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"])
100*(t0-t3)/t0

### Error on average HNA densities
e0 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="0"]^2))/3
e3 <- sqrt(sum(results$sd.HNA.cells[results$Treatment=="Feeding"][results$Time[results$Treatment=="Feeding"]=="3"]^2))/3

### Error on HNA removal ratio
100*sqrt((sqrt(e0^2 + e3^2)/(t0-t3))^2 + (e0/t0)^2)

