library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(ggpmisc)
library(pwr)
setwd("~/Crapaud_documents/R")


windowsFonts("Arial" = windowsFont("Arial"))



TE <- read.table("all_podospora_finalrepeatcounts.out", header=F)

colnames(TE) <- c("species","rep","superfamily","abundance")

TE$abundance_log <- log(TE$abundance)
TE$rep_corrected <- str_remove(TE$rep,"_i")

TE$species <- as.factor(TE$species)

theme_set(theme_bw())
options(ggrepel.max.overlaps = Inf)

positions <- c("PaWa63p","CBS237.71m","PcWa139m","CBS124.78p","CBS411.78m","CBS415.72m","CBS112042p")
name_list <- c("P.anserina", "P. pauciseta", "P. comata", "P. pseudoanserina", "P. pseudopauciseta", "P. pseudocomata", "P. bellae-mahoneyi")

Genome <- c("PaWa63p","CBS237.71m","PcWa139m","CBS124.78p","CBS411.78m","CBS415.72m","CBS112042p")
genome_size <- c(36002190,35576233,34605593,34930426,36170437,35065207,34469758)
Total_TE_bp <- c(sum(TE[TE$species=="PaWa63p",]$abundance),sum(TE[TE$species=="CBS237.71m",]$abundance),sum(TE[TE$species=="PcWa139m",]$abundance),sum(TE[TE$species=="CBS124.78p",]$abundance),sum(TE[TE$species=="CBS411.78m",]$abundance),sum(TE[TE$species=="CBS415.72m",]$abundance),sum(TE[TE$species=="CBS112042p",]$abundance))
Crapaud_bp <- c(579827,672789,314911,620742,906421,128887,554248)
Grenouille_bp <- c(426455,440010,49348,83307,654666,22923,92095)
species <- c("P. anserina","P. pauciseta", "P. comata", "P. pseudoanserina", "P. pseudopauciseta","P. pseudocomata","P. bellae-mahoneyi")

stats <- data.frame(Genome,genome_size,Total_TE_bp,Crapaud_bp,Grenouille_bp,species)

stats$genome_size <- stats$genome_size/1000000
colnames(stats) <- c("genome","Genome_size_Mbp","Total_TE_bp","Crapaud_bp","Grenouille_bp","species")


#### Code for Figure 1 ####

pA <- ggstripchart(TE, x = "species", y = "abundance", color = "superfamily", label = "rep",
                   size = 5,
                   xlab=FALSE,
                   ylab="Abundance (bp)",
                   repel = TRUE,
                   font.label = c(20,"bold"),
                   label.select = list(criteria = "`y` > 50000"),
                   jitter = 0.03,
                   order = positions
                   )+rotate_x_text(45)+ theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17, margin=margin(0,0,0.5,0, "cm"))+font("ylab", size = 18)

pB <- ggscatter(stats, x= "Total_TE_bp", y="Genome_size_Mbp",label = "species",
                      size = 4,
                      add = "reg.line",
                      add.params = list(color ="#D55E00"),
                      ylab="Genome size (Mbp)",
                      xlab ="Total TE (bp)",
                      conf.int = FALSE,
                      repel = TRUE,
                      font.family = "Arial",
                      font.label = c(16,"italic"),
                      #cor.coef = TRUE,
                      #cor.coeff.args = list(method="pearson", label.sep = '\n', size = 5)
)+theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17, margin=margin(0,0,0.5,0, "cm"))+font("ylab", size = 18)+font("xlab", size = 18)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=6)


## pC (non-crapaud bp vs genome size) was not used for final figure
#pC <- ggscatter(stats, x= "Non-Crapaud_bp", y="Genome_size_Mbp",label = "species",
              # size = 4,
               #add = "reg.line",
               #add.params = list(color ="#D55E00"),
               #conf.int = FALSE,
               #repel = TRUE,
              # font.family = "Arial",
              # font.label = c(16,"italic"),
               #cor.coef = TRUE,
               #cor.coeff.args = list(method="pearson", label.sep = '\n', size = 5)
#)+theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17, margin=margin(0,0,0.5,0, "cm"))+font("ylab", size = 18)+font("xlab", size = 18)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=6)

pD <- ggscatter(stats, x= "Crapaud_bp", y="Genome_size_Mbp",label = "species",
               size = 4,
               add = "reg.line",
               add.params = list(color ="#D55E00"),
               conf.int = FALSE,
               repel = TRUE,
               ylab = "Genome size (Mbp)",
               xlab = "crapaud (bp)",
               font.family = "Arial",
               font.label = c(16,"italic"),
               #cor.coef = TRUE,
               #cor.coeff.args = list(method="pearson", label.sep = '\n', size = 5)
)+theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 17,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 17, margin=margin(0,0,0.5,0, "cm"))+font("ylab", size = 18)+font("xlab", size = 18)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=6)



#### Code for Supplementary figure 1 #####



d <- ggscatter(stats, x= "Grenouille_bp", y="Genome_size_Mbp",label = "species",
               size = 3,
               add = "reg.line",
               add.params = list(color ="#D55E00"),
               conf.int = FALSE,
               repel = TRUE,
               font.family = "Arial",
               font.label = c(13,"italic"),
               xlab = "Grenouille (bp)",
               ylab = "Genome Size (Mbp)",
               cor.coef = TRUE,
               cor.coeff.args = list(method="pearson", label.sep = '\n', size = 5)
)+theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 20,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 20, margin=margin(0,0,0.5,0, "cm"))+font("ylab", size = 18)+font("xlab", size = 18)




#### Code for figure 3B ###

gc <- read.table("Full_solos_GC_variants_v2.txt")
colnames(gc) <- c("Sequence", "Subfamily", "GC %", "Type")


viol <- ggviolin(gc, x="Type", y="GC %",
                 size = 1,
                 fill = "gray50",
                 xlab = FALSE)+ theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 27,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 30, margin=margin(0,0,1,0, "cm"))+font("ylab", size = 30)

viol2 <- ggviolin(gc, x="Subfamily", y="GC %",
              size = 1,
              fill = "Subfamily",
              palette = c("#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#809693","#4FC601","#3B5DFF","#4A3B53","#CCCCCC"),
              breaks =3,
              order = c("LTR1","LTR2","LTR3","LTR4","LTR5","LTR6","LTR7","LTR8","LTR9","LTR10","LTR11","LTR12","LTR13","LTR14","LTR15","LTR16","Unknown"),xlab = FALSE)+rotate_x_text(45)+ theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))+font("y.text", size = 27,margin=margin(0,0,0,0.5,"cm"))+font("x.text", size = 30, margin=margin(0,0,1,0, "cm"))+font("ylab", size = 30)+scale_y_continuous(breaks = seq(0,60, by =20))+rremove("legend")


ggarrange(viol,viol2, nrow = 1, ncol = 2,align = c("h"))




## For Supplementary figure 11 ## 

greno <- read.table("Grenouille_LTR_all_final_proper_GC_16-6-2023.txt", header = F)
names(greno) <- c("copy", "GC %")
greno_full_list = c("_R_CBS237.71m_chromosome_5.2_621461-622020","_R_CBS411.78m_chromosome_6_175384-175932","_R_CBS411.78m_chromosome_4_577942-578497","CBS411.78m_chromosome_1_8999726-9000283","CBS411.78m_chromosome_4_4028265-4028820","","CBS411.78m_chromosome_3_3666814-3667369","CBS411.78m_chromosome_6_1531613-1532083","CBS411.78m_chromosome_6_2639187-2639744","CBS411.78m_chromosome_2_2771359-2771916","CBS411.78m_chromosome_3_4129408-4129965","_R_PaWa63p_chromosome_5_1157225-1157783","_R_CBS237.71m_chromosome_1_7138365-7138920","_R_CBS237.71m_chromosome_4_4008859-4009414","_R_CBS237.71m_chromosome_5.2_96770-97327","CBS237.71m_chromosome_4_3968148-3968704","CBS112042p_chromosome_1_6782118-6782674","_R_CBS237.71m_chromosome_7_922580-923108")
greno_full <- greno[greno$copy %in% greno_full_list,]
greno_solo <- subset(greno,!(copy %in% greno_full_list))
greno_full$type <- "Full"
greno_solo$type <- "Solo/Fragment"
greno_2 <- rbind(greno_solo,greno_full)

ggviolin(greno_2, y = "GC %", x="type", size = 1,add = "jitter", xlab = FALSE)+font("y.text",size=20)

