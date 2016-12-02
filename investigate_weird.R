### Positions of deviating replicates
faulty <- c("Q3_T_t0.5_rep1.fcs","Q1_T_t0_rep1.fcs","Q1_T_t0.5_rep1.fcs","C1_t1.5_rep1.fcs")

faulty <- c("Q3_T_t0.5_rep1.fcs","Q3_T_t0.5_rep2.fcs","Q3_T_t0.5_rep3.fcs")


pos <- c()
for (i in 1:length(faulty)) pos[i] <- which(sampleNames(flowData_transformed)==faulty[i])

ftest <- flowData_transformed[pos]
### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=ftest[1:3], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)

### 