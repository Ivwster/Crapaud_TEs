#Script to create the violin plots for supplementary figure 4

library(ggplot2)
library(ggpubr)

setwd("~/Crapaud_documents/R")




combined_stats <- read.table("All_SSNstats_and_GC.txt",header=TRUE,sep='\t')


mean(combined_stats[combined_stats$Type=="Outside",]$Degree)
mean(combined_stats[combined_stats$Type=="Inside",]$Degree)



hist(log(combined_stats[combined_stats$Type=="Outside",]$Degree))
hist(combined_stats[combined_stats$Type=="Inside",]$Degree)

qqnorm(log(combined_stats[combined_stats$Type=="Outside",]$Degree))
qqline(log(combined_stats[combined_stats$Type=="Outside",]$Degree))

wilcox.test(combined_stats[combined_stats$Type=="Outside",]$NumberOfUndirectedEdges, combined_stats[combined_stats$Type=="Inside",]$NumberOfUndirectedEdges, paired = FALSE)

pA <- ggviolin(combined_stats[combined_stats$Type != "unknown",], x="Type",y="NumberOfUndirectedEdges",fill="Type", ylab = "Number of Undirected Edges",
                xlab=FALSE
                
  
)+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17)+font("ylab", size = 18)+theme(legend.position = "none")

pB <- ggviolin(combined_stats[combined_stats$Type!="unknown",], x="Type",y="GC", fill = "Type",xlab=FALSE, ylab="GC (%)")+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17)+font("ylab", size = 18)+theme(legend.position = "none")

#+geom_jitter(aes(color = Full.solos_active))

ggarrange(pB,pA,labels = c("A","B"))

