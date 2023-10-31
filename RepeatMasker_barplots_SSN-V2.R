## R-script file for supplementary figure 6

library(ggplot2)
library(dplyr)
library(showtext)

setwd("~/Crapaud_documents/RepeatMasker_SSN-V2")

solos <- read.table("RMcountV3_Subfamilies_V2_noLTR11or12_counts_250filter.out", header=F)
colnames(solos) <- c("variant","strain","grp","nr")
solos$variant <- as.character(solos$variant)

solos$variant <- factor(solos$variant, levels=c("LTR1","LTR2","LTR3","LTR4","LTR5","LTR6","LTR7","LTR8","LTR9","LTR10","LTR11","LTR12","LTR13","LTR14"))
solos$strain <- factor(solos$strain, levels = c("PaWa63p","CBS237.71m","PcWa139m","CBS124.78p","CBS411.78m","CBS415.72m","CBS112042p"))



## The stacked barplots are plotted five at a time and then merged in inkscape 

ggplot(data=solos[solos$variant %in% c("LTR1","LTR2","LTR3","LTR4","LTR5"),], aes(x= strain,y=nr,fill=grp)) +
  geom_bar(stat="identity", position = "stack") +
  facet_grid(~ variant) +
  scale_y_continuous(limits= c(0,120), breaks = c(0,50,100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 1, size = 18),axis.title.x = element_blank())

ggplot(data=solos[solos$variant %in% c("LTR6","LTR7","LTR8","LTR9","LTR10"),], aes(x= strain,y=nr,fill=grp)) +
  geom_bar(stat="identity", position = "stack") +
  facet_grid(~ variant) +
  scale_y_continuous(limits= c(0,120), breaks = c(0,50,100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 1, size = 18),axis.title.x = element_blank())


## The LTR1 plot is removed from this and the rest are moved to the left for the final figure

ggplot(data=solos[solos$variant %in% c("LTR11","LTR12","LTR13","LTR14","LTR1"),], aes(x= strain,y=nr,fill=grp)) +
  geom_bar(stat="identity", position = "stack") +
  facet_grid(~ variant) +
  scale_y_continuous(limits= c(0,120), breaks = c(0,50,100))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 1, size = 18),axis.title.x = element_blank())







LTR1 <- solos[solos$variant == "LTR1",]
LTR2 <- solos[solos$variant == "LTR2",]
LTR3 <- solos[solos$variant == "LTR3",]
LTR4 <- solos[solos$variant == "LTR4",]
LTR5 <- solos[solos$variant == "LTR5",]
LTR6 <- solos[solos$variant == "LTR6",]
LTR7 <- solos[solos$variant == "LTR7",]
LTR8 <- solos[solos$variant == "LTR8",]
LTR9 <- solos[solos$variant == "LTR9",]
LTR10 <- solos[solos$variant == "LTR10",]
LTR11 <- solos[solos$variant == "LTR11",]
LTR12 <- solos[solos$variant == "LTR12",]
LTR13 <- solos[solos$variant == "LTR13",]
LTR14 <- solos[solos$variant == "LTR14",]

LTR1_cont <- xtabs(nr~strain+grp, data=LTR1)
LTR2_cont <- xtabs(nr~strain+grp, data=LTR2)
LTR3_cont <- xtabs(nr~strain+grp, data=LTR3)
LTR4_cont <- xtabs(nr~strain+grp, data=LTR4)
LTR5_cont <- xtabs(nr~strain+grp, data=LTR5)
LTR6_cont <- xtabs(nr~strain+grp, data=LTR6)
LTR7_cont <- xtabs(nr~strain+grp, data=LTR7)
LTR8_cont <- xtabs(nr~strain+grp, data=LTR8)
LTR9_cont <- xtabs(nr~strain+grp, data=LTR9)
LTR10_cont <- xtabs(nr~strain+grp, data=LTR10)
LTR11_cont <- xtabs(nr~strain+grp, data=LTR11)
LTR12_cont <- xtabs(nr~strain+grp, data=LTR12)
LTR13_cont <- xtabs(nr~strain+grp, data=LTR13)
LTR14_cont <- xtabs(nr~strain+grp, data=LTR14)

chisq.test(LTR3_cont)
