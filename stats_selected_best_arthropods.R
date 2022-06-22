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
  return(data[1,])
}
# CREATE DATAFRAME WITH TOP RESULT BY SPECIES
best.params.df <- lapply(temp, read_and_process)
best.params.df <- do.call(rbind.data.frame, best.params.df)
drop <- c("program", "SNi", "SPi", "SNichain", "SPichain", "SNt", "SPt", "SNg", "SPg", "Matching_intron_chains", "Missed_introns", "Reference_introns", "Novel_introns", "Predicted_introns", "Missed_loci", "Reference_loci", "Novel_loci", "Predicted_loci")
best.params.df <- best.params.df[,!names(best.params.df) %in% drop]
# RETRIEVE EXCEL WITH SPECIES INFO
selected.arthro <- read.csv2("../stats_arthropods_genomes.tsv", sep = "\t", header = TRUE)
names(selected.arthro)[names(selected.arthro) == "Species"] <- "species"
# MERGE BOTH DATAFRAMES IN ONE
arthro.best.params <- merge(selected.arthro, best.params.df, by = "species")
# ALL PARAMS TOGETHER
arthro.best.params <- within(arthro.best.params,  stacked <- paste(HSP, site_factor, exon_factor, exon_weight, sep="_"))
##########
# PART 2 #
##########
# Arthropods: sensitivity and specificity at nucleotide level
art.nucl.avg <- ggplot(arthro.best.params, aes(x =  exon_factor, y =nucl.mean, color=as.factor(HSP)))+
  geom_point(size = 2)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at nucleotide level") + ggtitle("Arthropods: sensitivity and specificity at nucleotide level")
art.nucl.avg
# DIFFERENTIAL ACCURACY BY GC CONTENT
arthro.best.params <- arthro.best.params%>%mutate(GC_group = case_when(
  GC_content<=mean(arthro.best.params$GC_content) ~ 'low_GC',
  GC_content>mean(arthro.best.params$GC_content) ~ 'high_GC'
))
p <- ggplot(arthro.best.params, aes(y =nucl.mean, color=GC_group))+
  geom_boxplot()+
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Arthropods nucleotide accuracy by GC content")
p1 <- ggplot(arthro.best.params, aes(y =exon.mean, color=GC_group))+
  geom_boxplot()+
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Arthropods exon accuracy by GC content")
grid.arrange(p, p1, ncol=1, nrow = 2)
gc.nucl.aov <- aov(nucl.mean ~ GC_group, data = arthro.best.params)
summary.aov(gc.nucl.aov)
pairwise.t.test(arthro.best.params$nucl.mean, arthro.best.params$GC_group, p.adjust="none", alternative=c("two.sided"))
gc.exon.aov <- aov(exon.mean ~ GC_group, data = arthro.best.params)
summary.aov(gc.exon.aov)
# DIFFERENTIAL ACCURACY BY GENOME SIZE
arthro.best.params <- arthro.best.params%>%mutate(Genome_size_group = case_when(
  Genome_size<=mean(arthro.best.params$Genome_size) ~ 'low_Gsize',
  Genome_size>mean(arthro.best.params$Genome_size) ~ 'high_Gsize'
))
h <- ggplot(arthro.best.params, aes(y =nucl.mean, color=Genome_size_group))+
  geom_boxplot()+
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Arthropods exon accuracy by genome size")
h1 <- ggplot(arthro.best.params, aes(y =exon.mean, color=Genome_size_group))+
  geom_boxplot()+
  theme_bw()+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Arthropods exon accuracy by genome size")
grid.arrange(h, h1, ncol=1, nrow = 2)
gs.nucl.aov <- aov(nucl.mean ~ Genome_size_group, data = arthro.best.params)
summary.aov(gs.nucl.aov)
pairwise.t.test(arthro.best.params$nucl.mean, arthro.best.params$Genome_size_group, p.adjust="none", alternative=c("two.sided"))
gs.exon.aov <- aov(exon.mean ~ Genome_size_group, data = arthro.best.params)
summary.aov(gs.exon.aov)
# HISTOGRAM COUNTS OF EACH PARAMETER
par(mfrow=c(2,2))
hist(x = arthro.best.params$exon_factor, main = "Hisogram exon factor", xlab = "Exon factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = arthro.best.params$exon_weight, main = "Histogram exon weigth", xlab = "Exon weight", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = arthro.best.params$site_factor, main = "Histogram site factor", xlab = "Site factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = arthro.best.params$HSP, main = "Histogram HSP factor", xlab = "HSP factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))
# PCA ANALYSIS
pca.params <- c("HSP", "exon_factor", "site_factor", "exon_weight")
pca.params.df <- arthro.best.params[,names(arthro.best.params) %in% pca.params]
pca_res <- prcomp(pca.params.df, scale. = TRUE)
autoplot(pca_res, data = arthro.best.params, colour = "GC_group")
# TABLE WITH COUNTS OF EACH PARAMETER SET
t <- table(arthro.best.params[,c("stacked")])
t <- sort(t, decreasing = TRUE)
# PREVIOUS PLOT IN GGPLOT
t.p <- ggplot(arthro.best.params)+
  geom_bar(aes(y=stacked), width = 0.5, colour="white", fill="black")+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Count sets of parameter")+
  ylab("Parameter combinations")
t.p
# MEAN AT NUCLEOTIDE LEVEL BY PARAMETER SET
params.id <- c(arthro.best.params$stacked)
nucl.mean.param <- (arthro.best.params$nucl.mean)
df.params.nucl <- data.frame(params.id, nucl.mean.param)
df.params.nucl <- aggregate(.~params.id, data = df.params.nucl, mean)
# PLOT MEAN AT NUCLEOTIDE LEVEL BY PARAMETER SET
ggplot(df.params.nucl)+
  geom_bar(aes(x = params.id, y = nucl.mean.param), stat = "identity")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 60, vjust = 0.9, hjust = 1),
        axis.text.y = element_text(size = 20))+
  ggtitle("Nucleotide accuracy by best parameter set")+ xlab("parameter set")+ ylab("Average nucleotide accuracy")
# MEAN AT EXON LEVEL BY PARAMETER SET
exon.mean.param <- (arthro.best.params$exon.mean)
df.params.exon <- data.frame(params.id, exon.mean.param)
df.params.exon <- aggregate(.~params.id, data = df.params.exon, mean)
# PLOT MEAN AT EXON LEVEL BY PARAMETER SET
ggplot(df.params.exon)+
  geom_bar(aes(x = params.id, y = exon.mean.param), stat = "identity")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Exon mean by parameter set")

# PLOT EXON AVERAGE BY PARAMETER SET AND COUNTS OF EACH PARAMETER SET -- MORE PRETTY
exon.test <- tibble(
  params.id <- c(df.params.exon$params.id),
  exon.mean.value <- c(df.params.exon$exon.mean.param),
  t <- c(t)
)

ggplot(arthro.best.params) + 
  geom_bar(mapping = aes(x = stacked), fill = "black") +
  geom_point(data = df.params.exon, aes(x = params.id, y = exon.mean.param), stat = "identity", colour = "blue")+
  #scale_x_date(name = "Parameters set") +
  theme_bw()+
  scale_y_continuous(name = "Parameters occurrences", 
                     sec.axis = sec_axis(~./1, name = "Average sp+sn exon level")) + 
  theme(
    axis.title.y = element_text(color = "black", size = 20),
    axis.title.y.right = element_text(color = "blue", size = 20),
    axis.text.y.right = element_text(color = "blue", size = 20),
    axis.text.y.left = element_text(size = 20),
    axis.ticks.y.right = element_line(color = "blue"),
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 60, vjust = 0.9, hjust = 1))
  
# PLOT NUCLEOTIDE AVERAGE BY PARAMETER SET AND COUNTS OF EACH PARAMETER SET -- MORE PRETTY
nucl.test <- tibble(
  params.id <- c(df.params.nucl$params.id),
  nucl.mean.value <- c(df.params.nucl$nucl.mean.param),
  t <- c(t)
)

ggplot(arthro.best.params) + 
  geom_bar(mapping = aes(x = stacked), fill = "black") +
  geom_point(data = df.params.nucl, aes(x = params.id, y = nucl.mean.param), stat = "identity", colour = "blue")+
  #scale_x_date(name = "Parameters set") +
  theme_bw()+
  scale_y_continuous(name = "Parameters occurrences", 
                     sec.axis = sec_axis(~./1, name = "Average sp+sn nucl level")) + 
  theme(
    axis.title.y = element_text(color = "black", size = 20),
    axis.title.y.right = element_text(color = "blue", size = 20),
    axis.text.y.right = element_text(color = "blue", size = 20),
    axis.text.y.left = element_text(size = 20),
    axis.ticks.y.right = element_line(color = "blue"),
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 60, vjust = 0.9, hjust = 1))
