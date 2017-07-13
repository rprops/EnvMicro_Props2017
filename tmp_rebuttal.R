data.total.final.b <- data.total.final[data.total.final$Lake == "Lake Michigan" | data.total.final$Lake == "Muskegon Lake", ]
data.total.final.b <- droplevels(data.total.final.b)


lmGrid <- expand.grid(intercept = TRUE)
R.cv.D2.b <- train(log2(D2)~log2(D2.fcm)*Lake, data=data.total.final.b, method ='rlm', trControl = trainControl(method ="repeatedcv", repeats = 100), tuneGrid = lmGrid)
R.cv.D2.b$results

# Slope is not significantly different between Muskegon Lake and Lake Michigan
# But intercept is.
# As a result the trends will be preserved between these two ecosystems.
Anova(rlm(log2(D2)~log2(D2.fcm)*Lake, data=data.total.final.b))

my_grob = grobTree(textGrob(bquote(r^2 == .(paste(round(R.cv.D2.b$results$Rsquared, 2)))), x=0.8,  y=0.16, hjust=0,
                            gp=gpar(col="black", fontsize=20, fontface="italic")))


p5.b <- ggplot(data=data.total.final.b,aes(x=D2.fcm,y=D2, fill=Lake))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21, size=6, alpha=0.6, aes(fill=Lake))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[2]),x=expression('Phenotypic diversity - D'[2]), fill="")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(5, 60, 10),minor_breaks = NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1750, 2900, 125),minor_breaks = NULL, limits = c(NA, 2900)) +
  geom_smooth(method="rlm",color="black",formula=y~x)+
  annotation_custom(my_grob)+
  guides(fill=FALSE)

pdf("Fig1_revisedBb.pdf", width = 7, height = 6)
print(p5.b)
dev.off()