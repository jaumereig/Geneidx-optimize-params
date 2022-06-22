setwd("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/arthropods/foreval_stats")
library("ggplot2")
library("dplyr")
library("ggfortify")
library("gridExtra")
##########
# PART 1 #
##########
# LOAD HEADER FILES IN FOLDER
temp = list.files(pattern="*.header")
#print(list.files(pattern="*.header"))
#for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header = TRUE, sep = '\t'))
# READ AND PROCESS MEANS IN EACH FILE, PICK ONLY THE TOP AVERAGE AT NUCLEOTIDE+EXON+INTRON LEVEL
read_and_process <- function(filename){
  data <- read.csv(filename, header = TRUE, sep = '\t')
  data$nucl.mean <- rowMeans(data[,c("SPn", "SNn")], na.rm=TRUE)
  data$exon.mean <- rowMeans(data[,c("SPe", "SNe")], na.rm=TRUE)
  data$intron.mean <- rowMeans(data[,c("SPi", "SNi")], na.rm=TRUE)
  data$all.mean <- rowMeans(data[,c("nucl.mean", "exon.mean", "intron.mean")], na.rm = TRUE)
  data = data[order(data$all.mean, decreasing = TRUE),]
  return(data[data$all.mean > (max(data$all.mean)*0.97),])
}
# CREATE DATAFRAME WITH TOP RESULT BY SPECIES
top.params.df <- lapply(temp, read_and_process)
top.params.df <- do.call(rbind.data.frame, top.params.df)
drop <- c("program", "SNi", "SPi", "SNichain", "SPichain", "SNt", "SPt", "SNg", "SPg", "Matching_intron_chains", "Missed_introns", "Reference_introns", "Novel_introns", "Predicted_introns", "Missed_loci", "Reference_loci", "Novel_loci", "Predicted_loci")
top.params.df <- top.params.df[,!names(top.params.df) %in% drop]
# RETRIEVE EXCEL WITH SPECIES INFO
selected.arthro <- read.csv2("../stats_arthropods_genomes.tsv", sep = "\t", header = TRUE)
names(selected.arthro)[names(selected.arthro) == "Species"] <- "species"
# MERGE BOTH DATAFRAMES IN ONE
arthro.params.top <- merge(selected.arthro, top.params.df, by = "species")
##########
# PART 2 #
##########
# ALL PARAMS TOGETHER
arthro.params.top <- within(arthro.params.top,  stacked <- paste(HSP, site_factor, exon_factor, exon_weight, sep="_"))
plot(table(arthro.params.top$species))
plot(table(arthro.params.top$stacked)) # 63 out of 81 parameters
tail(sort(table(arthro.params.top$stacked)))
ggplot(arthro.params.top, aes(y = stacked))+
  geom_bar(width = 0.5, colour="white", fill="black")+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 12))+
  ggtitle("Count sets of parameters")+
  xlab("Parameter combinations")+
  ylab("Parameter set")
# DIFFERENTIAL PARAMETER SETS COUNTS BY GC CONTENT AND GENOME SIZE
arthro.params.top <- arthro.params.top%>%mutate(GC_group = case_when(
  GC_content<=median(arthro.best.params$GC_content) ~ 'low_GC',
  GC_content>median(arthro.best.params$GC_content) ~ 'high_GC'
))
arthro.params.top <- arthro.params.top%>%mutate(GenomeSize_group = case_when(
  Genome_size<=median(arthro.best.params$Genome_size) ~ 'low_Gsize',
  Genome_size>median(arthro.best.params$Genome_size) ~ 'high_Gsize'
))
gc <- ggplot(data = arthro.params.top, aes(x = stacked, fill = GC_group))+
  geom_bar(position = "stack")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 60, vjust = 0.9, hjust = 1),
        axis.text.y = element_text(size = 20))+
  ggtitle("Count sets of parameters by GC group")+
  ylab("Parameter counts")+
  xlab("Parameter set")
gs <- ggplot(data = arthro.params.top, aes(x = stacked, fill = GenomeSize_group))+
  geom_bar(position = "stack")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 60, vjust = 0.9, hjust = 1),
        axis.text.y = element_text(size = 20))+
  ggtitle("Count sets of parameters by Genome size")+
  ylab("Parameter counts")+
  xlab("Parameter set")
grid.arrange(gc, gs, nrow = 2, ncol = 1)
table(arthro.params.top$GC_group)
table(arthro.params.top$GenomeSize_group)

gs.nucl.top.aov <- aov(nucl.mean ~ GenomeSize_group, data = arthro.params.top)
summary.aov(gs.nucl.top.aov)
pairwise.t.test(arthro.params.top$nucl.mean, arthro.params.top$GenomeSize_group, p.adjust="none", alternative=c("greater"))
pairwise.t.test(arthro.params.top$exon.mean, arthro.params.top$GenomeSize_group, p.adjust="none", alternative=c("less"))
gs.exon.top.aov <- aov(exon.mean ~ GenomeSize_group, data = arthro.params.top)
summary.aov(gs.exon.top.aov)

gc.nucl.top.aov <- aov(nucl.mean ~ GC_group, data = arthro.params.top)
summary.aov(gc.nucl.top.aov)
pairwise.t.test(arthro.params.top$nucl.mean, arthro.params.top$GC_group, p.adjust="none", alternative=c("greater"))
pairwise.t.test(arthro.params.top$exon.mean, arthro.params.top$GC_group, p.adjust="none", alternative=c("less"))
gc.exon.top.aov <- aov(exon.mean ~ GC_group, data = arthro.params.top)
summary.aov(gc.exon.top.aov)

# DIFFERENTIAL ACCURACY BY GC CONTENT AND GENOME SIZE
p.gc <- ggplot(arthro.params.top, aes(y =nucl.mean, color=GC_group))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))
p1.gc <- ggplot(arthro.params.top, aes(y =exon.mean, color=GC_group))+
  geom_boxplot()+
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))
grid.arrange(p.gc, p1.gc, ncol=1, nrow = 2)

h.gs <- ggplot(arthro.params.top, aes(y =nucl.mean, color=GenomeSize_group))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))
h1.gs <- ggplot(arthro.params.top, aes(y =exon.mean, color=GenomeSize_group))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))
grid.arrange(h.gs, h1.gs, ncol=1, nrow = 2)
