library("ggplot2")
##################
#LOAD STATS FILES#
##################
anopheles.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
apis.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Apis_mellifera.AmelHAv.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
drosophila.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Drosophila_melanogaster.BDGP6.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
bombus.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Bombus_terrestris.Bter.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
elegans.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
solenopsis.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Solenopsis_invicta.UNIL.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
funestus.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_funestus.AfunF3.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
bemisia.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Bemisia_tabaci.ASIAII5.dna_rm.toplevel.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
glossina.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Glossina_austeni.GausT1.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
danaus.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Danaus_plexippus.Dplex2.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
necator.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Necator_americanus.Namericanus3.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
erecta.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Drosophila_erecta.dere.dna.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
mojavensis.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Drosophila_mojavensis.dmojcaf.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
aedes.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Aedes_aegypti.AaegL5.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.range.header", header = TRUE)
####################
#AVERAGE ALL LEVELS#
####################
# ANOPHELES
anopheles.df$nucl.mean.anopheles.df <- rowMeans(anopheles.df[,c("SPn", "SNn")], na.rm=TRUE)
anopheles.df$exon.mean.anopheles.df <- rowMeans(anopheles.df[,c("SPe", "SNe")], na.rm=TRUE)
anopheles.df$intron.mean.anopheles.df <- rowMeans(anopheles.df[,c("SPi", "SNi")], na.rm=TRUE)
anopheles.df$all.mean.anopheles.df <- rowMeans(anopheles.df[,c("nucl.mean.anopheles.df", "exon.mean.anopheles.df", "intron.mean.anopheles.df")], na.rm = TRUE)
# APIS
apis.df$nucl.mean.apis.df <- rowMeans(apis.df[,c("SPn", "SNn")], na.rm=TRUE)
apis.df$exon.mean.apis.df <- rowMeans(apis.df[,c("SPe", "SNe")], na.rm=TRUE)
apis.df$intron.mean.apis.df <- rowMeans(apis.df[,c("SPi", "SNi")], na.rm=TRUE)
apis.df$all.mean.apis.df <- rowMeans(apis.df[,c("nucl.mean.apis.df", "exon.mean.apis.df", "intron.mean.apis.df")], na.rm = TRUE)
# DROSOPHILA
drosophila.df$nucl.mean.drosophila.df <- rowMeans(drosophila.df[,c("SPn", "SNn")], na.rm=TRUE)
drosophila.df$exon.mean.drosophila.df <- rowMeans(drosophila.df[,c("SPe", "SNe")], na.rm=TRUE)
drosophila.df$intron.mean.drosophila.df <- rowMeans(drosophila.df[,c("SPi", "SNi")], na.rm=TRUE)
drosophila.df$all.mean.drosophila.df <- rowMeans(drosophila.df[,c("nucl.mean.drosophila.df", "exon.mean.drosophila.df", "intron.mean.drosophila.df")], na.rm = TRUE)
#BOMBUS
bombus.df$nucl.mean.bombus.df <- rowMeans(bombus.df[,c("SPn", "SNn")], na.rm=TRUE)
bombus.df$exon.mean.bombus.df <- rowMeans(bombus.df[,c("SPe", "SNe")], na.rm=TRUE)
bombus.df$intron.mean.bombus.df <- rowMeans(bombus.df[,c("SPi", "SNi")], na.rm=TRUE)
bombus.df$all.mean.bombus.df <- rowMeans(bombus.df[,c("nucl.mean.bombus.df", "exon.mean.bombus.df", "intron.mean.bombus.df")], na.rm = TRUE)
# ELEGANS
elegans.df$nucl.mean.elegans.df <- rowMeans(elegans.df[,c("SPn", "SNn")], na.rm=TRUE)
elegans.df$exon.mean.elegans.df <- rowMeans(elegans.df[,c("SPe", "SNe")], na.rm=TRUE)
elegans.df$intron.mean.elegans.df <- rowMeans(elegans.df[,c("SPi", "SNi")], na.rm=TRUE)
elegans.df$all.mean.elegans.df <- rowMeans(elegans.df[,c("nucl.mean.elegans.df", "exon.mean.elegans.df", "intron.mean.elegans.df")], na.rm = TRUE)
# SOLENOPSIS
solenopsis.df$nucl.mean.solenopsis.df <- rowMeans(solenopsis.df[,c("SPn", "SNn")], na.rm=TRUE)
solenopsis.df$exon.mean.solenopsis.df <- rowMeans(solenopsis.df[,c("SPe", "SNe")], na.rm=TRUE)
solenopsis.df$intron.mean.solenopsis.df <- rowMeans(solenopsis.df[,c("SPi", "SNi")], na.rm=TRUE)
solenopsis.df$all.mean.solenopsis.df <- rowMeans(solenopsis.df[,c("nucl.mean.solenopsis.df", "exon.mean.solenopsis.df", "intron.mean.solenopsis.df")], na.rm = TRUE)
# FUNESTUS
funestus.df$nucl.mean.funestus.df <- rowMeans(funestus.df[,c("SPn", "SNn")], na.rm=TRUE)
funestus.df$exon.mean.funestus.df <- rowMeans(funestus.df[,c("SPe", "SNe")], na.rm=TRUE)
funestus.df$intron.mean.funestus.df <- rowMeans(funestus.df[,c("SPi", "SNi")], na.rm=TRUE)
funestus.df$all.mean.funestus.df <- rowMeans(funestus.df[,c("nucl.mean.funestus.df", "exon.mean.funestus.df", "intron.mean.funestus.df")], na.rm = TRUE)
# BEMISIA
bemisia.df$nucl.mean.bemisia.df <- rowMeans(bemisia.df[,c("SPn", "SNn")], na.rm=TRUE)
bemisia.df$exon.mean.bemisia.df <- rowMeans(bemisia.df[,c("SPe", "SNe")], na.rm=TRUE)
bemisia.df$intron.mean.bemisia.df <- rowMeans(bemisia.df[,c("SPi", "SNi")], na.rm=TRUE)
bemisia.df$all.mean.bemisia.df <- rowMeans(bemisia.df[,c("nucl.mean.bemisia.df", "exon.mean.bemisia.df", "intron.mean.bemisia.df")], na.rm = TRUE)
# GLOSSINA
glossina.df$nucl.mean.glossina.df <- rowMeans(glossina.df[,c("SPn", "SNn")], na.rm=TRUE)
glossina.df$exon.mean.glossina.df <- rowMeans(glossina.df[,c("SPe", "SNe")], na.rm=TRUE)
glossina.df$intron.mean.glossina.df <- rowMeans(glossina.df[,c("SPi", "SNi")], na.rm=TRUE)
glossina.df$all.mean.glossina.df <- rowMeans(glossina.df[,c("nucl.mean.glossina.df", "exon.mean.glossina.df", "intron.mean.glossina.df")], na.rm = TRUE)
# DANAUS
danaus.df$nucl.mean.danaus.df <- rowMeans(danaus.df[,c("SPn", "SNn")], na.rm=TRUE)
danaus.df$exon.mean.danaus.df <- rowMeans(danaus.df[,c("SPe", "SNe")], na.rm=TRUE)
danaus.df$intron.mean.danaus.df <- rowMeans(danaus.df[,c("SPi", "SNi")], na.rm=TRUE)
danaus.df$all.mean.danaus.df <- rowMeans(danaus.df[,c("nucl.mean.danaus.df", "exon.mean.danaus.df", "intron.mean.danaus.df")], na.rm = TRUE)
# NECATOR
necator.df$nucl.mean.necator.df <- rowMeans(necator.df[,c("SPn", "SNn")], na.rm=TRUE)
necator.df$exon.mean.necator.df <- rowMeans(necator.df[,c("SPe", "SNe")], na.rm=TRUE)
necator.df$intron.mean.necator.df <- rowMeans(necator.df[,c("SPi", "SNi")], na.rm=TRUE)
necator.df$all.mean.necator.df <- rowMeans(necator.df[,c("nucl.mean.necator.df", "exon.mean.necator.df", "intron.mean.necator.df")], na.rm = TRUE)
# ERECTA
erecta.df$nucl.mean.erecta.df <- rowMeans(erecta.df[,c("SPn", "SNn")], na.rm=TRUE)
erecta.df$exon.mean.erecta.df <- rowMeans(erecta.df[,c("SPe", "SNe")], na.rm=TRUE)
erecta.df$intron.mean.erecta.df <- rowMeans(erecta.df[,c("SPi", "SNi")], na.rm=TRUE)
erecta.df$all.mean.erecta.df <- rowMeans(erecta.df[,c("nucl.mean.erecta.df", "exon.mean.erecta.df", "intron.mean.erecta.df")], na.rm = TRUE)
for (val in erecta.df$HSP) {
  if (val == 'toplevel.0.4') {
    erecta.df$HSP[erecta.df$HSP == val] <- 0.4
  }else if (val == 'toplevel.0.6') {
    erecta.df$HSP[erecta.df$HSP == val] <- 0.6
  }else if (val == 'toplevel.0.8') {
    erecta.df$HSP[erecta.df$HSP == val] <- 0.8
  }else if (val == 'toplevel.1.2') {
    erecta.df$HSP[erecta.df$HSP == val] <- 1.2
  }
}
erecta.df <- erecta.df[erecta.df$Predicted_exons != 0, ]
# MOJAVENSIS
mojavensis.df$nucl.mean.mojavensis.df <- rowMeans(mojavensis.df[,c("SPn", "SNn")], na.rm=TRUE)
mojavensis.df$exon.mean.mojavensis.df <- rowMeans(mojavensis.df[,c("SPe", "SNe")], na.rm=TRUE)
mojavensis.df$intron.mean.mojavensis.df <- rowMeans(mojavensis.df[,c("SPi", "SNi")], na.rm=TRUE)
mojavensis.df$all.mean.mojavensis.df <- rowMeans(mojavensis.df[,c("nucl.mean.mojavensis.df", "exon.mean.mojavensis.df", "intron.mean.mojavensis.df")], na.rm = TRUE)
for (val in mojavensis.df$HSP) {
  if (val == 'toplevel.0.4') {
    mojavensis.df$HSP[mojavensis.df$HSP == val] <- 0.4
  }else if (val == 'toplevel.0.6') {
    mojavensis.df$HSP[mojavensis.df$HSP == val] <- 0.6
  }else if (val == 'toplevel.0.8') {
    mojavensis.df$HSP[mojavensis.df$HSP == val] <- 0.8
  }else if (val == 'toplevel.1.2') {
    mojavensis.df$HSP[mojavensis.df$HSP == val] <- 1.2
  }
}
mojavensis.df <- mojavensis.df[mojavensis.df$Predicted_exons != 0, ]
# AEDES
aedes.df$nucl.mean.aedes.df <- rowMeans(aedes.df[,c("SPn", "SNn")], na.rm=TRUE)
aedes.df$exon.mean.aedes.df <- rowMeans(aedes.df[,c("SPe", "SNe")], na.rm=TRUE)
aedes.df$intron.mean.aedes.df <- rowMeans(aedes.df[,c("SPi", "SNi")], na.rm=TRUE)
aedes.df$all.mean.aedes.df <- rowMeans(aedes.df[,c("nucl.mean.aedes.df", "exon.mean.aedes.df", "intron.mean.aedes.df")], na.rm = TRUE)
for (val in aedes.df$HSP) {
  if (val == 'toplevel.0.4') {
    aedes.df$HSP[aedes.df$HSP == val] <- 0.4
  }else if (val == 'toplevel.0.6') {
    aedes.df$HSP[aedes.df$HSP == val] <- 0.6
  }else if (val == 'toplevel.0.8') {
    aedes.df$HSP[aedes.df$HSP == val] <- 0.8
  }else if (val == 'toplevel.1.2') {
    aedes.df$HSP[aedes.df$HSP == val] <- 1.2
  }
}
aedes.df <- aedes.df[aedes.df$Predicted_exons != 0, ]

#################
#FORMATTING DATA#
#################
# ANOPHELES
anopheles.df$exon_weight = anopheles.df$exon_weight*(-1)
exon_w.anopheles.df <- as.factor(anopheles.df$exon_weight)
site_f.anopheles.df <- as.factor(anopheles.df$site_factor)
exon_f.anopheles.df <- as.factor(anopheles.df$exon_factor)
hsp <- as.factor(anopheles.df$HSP)
# APIS
apis.df$exon_weight = apis.df$exon_weight*(-1)
exon_w.apis.df <- as.factor(apis.df$exon_weight)
site_f.apis.df <- as.factor(apis.df$site_factor)
exon_f.apis.df <- as.factor(apis.df$exon_factor)
hsp <- as.factor(apis.df$HSP)
# DROSOPHILA

#BOMBUS

# ELEGANS

# SOLENOPSIS

############
#MAKE PLOTS#
############
# ANOPHELES
anop.1 <- ggplot(anopheles.df, aes(x = exon_w.anopheles.df, y =nucl.mean.anopheles.df, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at nucleotide level") + ggtitle("Anopheles gambiae: sensitivity and specificity at nucleotide level")
anop.1
anop.2 <- ggplot(anopheles.df, aes(x = exon_w.anopheles.df, y =exon.mean.anopheles.df, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at exon level") + ggtitle("Anopheles gambiae: sensitivity and specificity at exon level")
anop.2
anop.3 <- ggplot(anopheles.df, aes(x = exon_w.anopheles.df, y =Predicted_exons, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Number of exons") + ggtitle("Anopheles gambiae: number of predicted exons")
anop.3

anop.4 <- ggplot(anopheles.df, aes(x = exon_w.anopheles.df, y =Predicted_loci, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Number of loci") + ggtitle("Anopheles gambiae: number of predicted loci")
anop.4
# APIS
apis.1 <- ggplot(apis.df, aes(x = exon_w.apis.df, y =nucl.mean.apis.df, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at nucleotide level") + ggtitle("Apis mellifera: sensitivity and specificity at nucleotide level")
apis.1
apis.2 <- ggplot(apis.df, aes(x = exon_w.apis.df, y =exon.mean.apis.df, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at exon level") + ggtitle("Apis mellifera: sensitivity and specificity at exon level")
apis.2
apis.3 <- ggplot(apis.df, aes(x = exon_w.apis.df, y =Predicted_exons, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Number of exons") + ggtitle("Apis mellifera: number of predicted exons")
apis.3
apis.4 <- ggplot(apis.df, aes(x = exon_w.apis.df, y =Predicted_loci, color=hsp))+
  geom_point(size = 3)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), 
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20),
        axis.text.y = element_text(size = 20))+
  xlab("Exon weight")+ ylab("Number of loci") + ggtitle("Apis mellifera: number of predicted loci")
apis.4
# DROSOPHILA

#BOMBUS

# ELEGANS

# SOLENOPSIS
