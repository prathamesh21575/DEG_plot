Data=read.csv("C:\\Users\\prath\\OneDrive\\Documents\\data3.csv",
              row.names = 1)
Data = log2(Data+ 1)
group = c('DOX','DOX', 'DOX', 'CTRL', 'CTRL', 'DOX', 'CTRL', 'CTRL')
#pvalue meanko meanwt log2fc
p_values = vector()
log2_fold_changes = vector()
mat=matrix(NA,ncol=4,nrow=nrow(Data))
colnames(mat)=c('p_values','meanDOX','meanCTRL','log2fc')
for (i in 1:nrow(Data)) 
{
  gene=as.numeric(Data[i,])
  df=cbind.data.frame(gene,group)
  T = t.test(gene~group, data=df, paired=F, alternative = 'two.sided')
  mat[i,'p_values']=T$p.value
  mat[i,c('meanDOX','meanCTRL')]=T$estimate
  mat[i,'log2fc']=log2(mat[i,'meanDOX']/mat[i,'meanCTRL'])
}
x=mat[,'log2fc']
y=-log10(mat[,'p_values'])
plot(x,y)

library(tidyverse) 
library(RColorBrewer)
library(ggrepel)
library(ggplot2)

mdf=as.data.frame(mat)

significance_threshold = 0.05

mdf$diffexpressed = "NO"
mdf$diffexpressed[mdf$log2fc > 0.6 & mdf$p_values < 0.05] <- "UP"
mdf$diffexpressed[mdf$log2fc < -0.6 & mdf$p_values < 0.05] <- "DOWN"
head(mdf[order(mdf$meanDOX) & mdf$diffexpressed == 'DOWN', ])
mdf$delabel = ifelse(mdf$p_values %in% head(mdf[order(mdf$meanDOX), "pvalue"], 30), mdf$p_values, NA)

myvolcanoplot = ggplot(data = mdf, aes(x , y , col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("green", "black", "#bb0c00"),  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 6), xlim = c(-6, 6)) + 
  labs(color = 'legends',  
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  ggtitle('volcano plot') +  
  geom_text_repel(max.overlaps = Inf)  

myvolcanoplot
dev.off()
