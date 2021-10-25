#NTRprism spectrum plotting
#Nicolas Altemose, 2021

library(ggplot2)

args=commandArgs(TRUE)
perlfile=args[2]
i=args[3] #region
span = as.integer(args[4])
prefix=args[5]

data = read.table(perlfile)
heat = as.matrix(data[,5:dim(data)[2]])
heat = heat/rowSums(heat)
heat = heat[,1:span]
csums = colSums(heat)

topten = order(csums[1:(length(csums)-1)],decreasing=T)[1:10]
print(i)
labelstring=""
for(j in 1:10){
	print(paste0(j,": ",topten[j]))
	labelstring=paste0(labelstring,"\n",topten[j])
}

coldata = data.frame("Interval"=1:(length(csums)-1),"Density"=csums[1:(length(csums)-1)]/dim(data)[1])
plot.new()
p2 = ggplot(coldata, aes(Interval, Density)) + geom_line() +theme_minimal() + coord_cartesian(xlim=c(0,span),clip="off") + ggtitle(paste(i))+
		theme(plot.margin = unit(c(1, 5, 1, 1), "lines")) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
		annotate("text", label = labelstring, x = span+span/20, y = 0.2*csums[topten[1]]/dim(data)[1], size = 3, colour = "black",hjust = 0,vjust=0)
ggsave(paste0(prefix,".",i,".spectrum.png"), p2, device="png", width=8, height=3)
