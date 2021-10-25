#NTRprism heatmap plotting 
#Nicolas Altemose, 2021

library(ggplot2)
library(reshape2)

args=commandArgs(TRUE)
perlfile=args[2]
i=args[3]
spanX = as.integer(args[4])
binsize = as.integer(args[5])
prefix = args[6]

data = read.table(perlfile)
heat = as.matrix(data[,5:dim(data)[2]])
heat = heat/rowSums(heat)
heat = heat[,1:(spanX/binsize)]
csums = colSums(heat)

spanY = dim(heat)[1]



colnames(heat) = seq(1,spanX,binsize)
rownames(heat) = dim(heat)[1]:1
gmat = melt(t(heat[spanY:1,]),value.name="density")

p2=ggplot(gmat, aes(Var1, Var2, fill= density)) + 
  geom_raster(hjust=1) +
  scale_fill_gradient(low=rgb(25/255,25/255,59/255), high="yellow")+theme_minimal()+theme(axis.text.y=element_blank())+
  xlab("interval length (bp)")+ylab(paste("top",spanY,"kmers, descending freq."))+ggtitle(i)

ggsave(paste0(prefix,".",i,".heatmap.png"), p2, device="png", width=5, height=3)
