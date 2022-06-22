library(ggplot2)
library(dplyr)
##############################
### PARAMS THAT CAN CHANGE ###
###     HSP                ###
###     EXON FACTOR        ###
###     SITE FACTOR        ###
###     NO_SCORE           ###
###     EXON WEIGHT        ###
##############################
###     ANOPHELES          ###
##############################
### PART 1: NO_SCORE = 0 AND EXON_W = -4.5 ARE CONSTANT VALUES FOR A RANGE OF HSP (0.1, 0.2) AND EXON_F, SITE_F (O.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
stats <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.tsv.header", header = TRUE)
attach(stats)
HSP <- as.factor(stats$HSP)
exon_f <- as.factor(stats$exon_factor)
site_f <- as.factor(stats$site_factor)

ggplot(stats, aes(x = site_f, y =SPn, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP)
# Lower HSP and exon factor better specificity at nucleotide level. Highest points reached with middle values of site factor: 0.4
# Data presents big margin between top results and the rest: clear optimal combination for SPn: 0.1, 0.1, 0.4
ggplot(stats, aes(x = site_f, y =SNn, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP) 
# Higher HSP, exon factor and site factor implies a jump on the y axis (rise of sensitivity at nucleotide level). 
# Exon factor .1, .2, .3 are significantly lower. Exon factor at levels .6, .7 show the same values (kind of a plateau maybe)
ggplot(stats, aes(x = site_f, y =SNe, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP)
# HSP 0.1 0.2 yields same results. Increasing site factor increases sensitivity at exon level.
# Exon factor at 0.5 0.6 has the highest results.
ggplot(stats, aes(x = site_f, y =SPe, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP)
# Lower HSP (0.1 > 0.2) and exon factor (0.2 > 0.1) provide the best results given that site factor is at 0.5
ggplot(stats, aes(x = site_f, y =SNi, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP)
# HSP 0.2 better than 0.1 but very similar. Exon factor at 0.5 0.6 present the best results with 0.7 extremely close. Highest site factor the best.
ggplot(stats, aes(x = site_f, y =SPi, color=exon_f ))+
  geom_point()+
  facet_wrap(~ HSP)
# HSP 0.2 better than 0.1 but very similar. Exon factor at 0.2 has the higher percentage combined with site factor 0.4 or 0.5

# PART 1.1: Average percentages of each level (exon, intron, nucleotide) and repeat analysis
stats$nucl.mean <- rowMeans(stats[,c("SPn", "SNn")], na.rm=TRUE)
stats$exon.mean <- rowMeans(stats[,c("SPe", "SNe")], na.rm=TRUE)
stats$intron.mean <- rowMeans(stats[,c("SPi", "SNi")], na.rm=TRUE)
ggplot(stats, aes(x = site_f, y =nucl.mean, color=exon_f ))+ #> max(stats$nucl.mean) [1] 79.35
  geom_point()+
  facet_wrap(~ HSP)
# site= 0.7, hsp= 0.2, exon= 0.3, 0.4
ggplot(stats, aes(x = site_f, y =exon.mean, color=exon_f ))+ # 50.4
  geom_point()+
  facet_wrap(~ HSP)
# site= 0.6, hsp= 0.1, exon= 0.3
ggplot(stats, aes(x = site_f, y =intron.mean, color=exon_f ))+ # 53.45
  geom_point()+
  facet_wrap(~ HSP)
# site= 0.6, 0.7, hsp= 0.2, exon= 0.3, 0.4

# PART 2: NO_SCORE = 0 AND EXON_W = -4.5 AND HSP = 0.7 ARE CONSTANT VALUES FOR A RANGE OF EXON_F, SITE_F (O.1, 0.2, 0.3, 0.4, 0.5, 0.6) 
stats.part.2 <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.part2.header", header = TRUE)
attach(stats.part.2)
HSP.part.2 <- as.factor(stats.part.2$HSP)
exon_f.part.2 <- as.factor(stats.part.2$exon_factor)
site_f.part.2 <- as.factor(stats.part.2$site_factor)
stats.part.2$nucl.mean.part.2 <- rowMeans(stats.part.2[,c("SPn", "SNn")], na.rm=TRUE)
stats.part.2$exon.mean.part.2 <- rowMeans(stats.part.2[,c("SPe", "SNe")], na.rm=TRUE)
stats.part.2$intron.mean.part.2 <- rowMeans(stats.part.2[,c("SPi", "SNi")], na.rm=TRUE)

ggplot(stats.part.2, aes(x = site_f.part.2, y =nucl.mean.part.2, color=exon_f.part.2 ))+ # 80
  geom_point()
# site= 0.6, exon= 0.4

ggplot(stats.part.2, aes(x = site_f.part.2, y =exon.mean.part.2, color=exon_f.part.2 ))+ # 48.4
  geom_point()
# site= 0.6, exon= 0.3, 0.4

ggplot(stats.part.2, aes(x = site_f.part.2, y =intron.mean.part.2, color=exon_f.part.2 ))+ # 52.95
  geom_point()
# site= 0.6, exon= 0.4, 0.3
# Here is confirmed that for any HSP, and keeping exon_w and NO_SCORE constant, the site factor needs to be higher in order to result in better percentages 
# and the exon factor also needs to be around 0.3 0.4. Seems like HSP 0.1 performs slightly better than 0.7 but not a significant difference

# PART 3: NO_SCORE = 0 HSP = 0. 2 ARE CONSTANT VALUES FOR A RANGE OF EXON_F (0.2, 0.3, 0.4, 0.5), SITE_F (0.5, 0.6, 0.7), EXON_W (-0.5, -1.5, -2.5, -3.5, -5.5)
stats.part.3 <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.part3.header", header = TRUE)
attach(stats.part.3)
exon_w.part.3 <- as.factor(stats.part.3$exon_weight)
exon_f.part.3 <- as.factor(stats.part.3$exon_factor)
site_f.part.3 <- as.factor(stats.part.3$site_factor)
stats.part.3$nucl.mean.part.3 <- rowMeans(stats.part.3[,c("SPn", "SNn")], na.rm=TRUE)
stats.part.3$exon.mean.part.3 <- rowMeans(stats.part.3[,c("SPe", "SNe")], na.rm=TRUE)
stats.part.3$intron.mean.part.3 <- rowMeans(stats.part.3[,c("SPi", "SNi")], na.rm=TRUE)

ggplot(stats.part.3, aes(x = exon_f.part.3, y =nucl.mean.part.3, color=exon_w.part.3 ))+
  geom_point()+
  facet_wrap(~ stats.part.3$site_factor)
# exon_w= -3.5 exon_f = 0.3 site_F = 0.6

ggplot(stats.part.3, aes(x = exon_f.part.3, y =exon.mean.part.3, color=exon_w.part.3 ))+
  geom_point()+
  facet_wrap(~ site_factor)
# exon_w= -4.5 exon_f = 0.3 site = 0.6

ggplot(stats.part.3, aes(x = exon_f.part.3, y =intron.mean.part.3, color=exon_w.part.3))+
  geom_point()+
  facet_wrap(~ site_factor)
# exon_w= -3.5 exon_f = 0.3 site= 0.6

# PART 4: EXON_F = 0.3 AND HSP = 0.2 AND SITE_F = 0.6 ARE CONSTANT VALUES FOR A RANGE OF EXON_W (-3.5, -4.5), NO_SCORE (0.0, -0.2, -0.4, -0.8, -1.6, -3.2)
stats.part.4 <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.part4.header", header = TRUE)
attach(stats.part.4)
exon_w.part.4 <- as.factor(stats.part.4$exon_weight)
score.part.4 <- as.factor(stats.part.4$score)
stats.part.4$nucl.mean.part.4 <- rowMeans(stats.part.4[,c("SPn", "SNn")], na.rm=TRUE)
stats.part.4$exon.mean.part.4 <- rowMeans(stats.part.4[,c("SPe", "SNe")], na.rm=TRUE)
stats.part.4$intron.mean.part.4 <- rowMeans(stats.part.4[,c("SPi", "SNi")], na.rm=TRUE)

ggplot(stats.part.4, aes(x = score.part.4, y =nucl.mean.part.4, color=exon_w.part.4))+
  geom_point()
# Score=0 is the best option clearly. With better results at exon_w= -3.5
ggplot(stats.part.4, aes(x = score.part.4, y =exon.mean.part.4, color=exon_w.part.4))+
  geom_point()
# Score=0 is the best option clearly. With better results at exon_w= -4.5
ggplot(stats.part.4, aes(x = score.part.4, y =intron.mean.part.4, color=exon_w.part.4))+
  geom_point()
# Score=0 is the best option clearly. With better results at exon_w= -3.5

# PART 5: NO_SCORE = 0.0 AND SITE_F = 0.6 AND EXON_F = 0.3 ARE CONSTANT VALUES FOR A RANGE OF EXON_W (-3.5 -4.0 -4.5), HSP (0.1 0.2 0.4 0.6 0.8)
stats.part.5 <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.part5.header", header = TRUE)
attach(stats.part.5)
exon_w.part.5 <- as.factor(stats.part.5$exon_weight)
hsp.part.5 <- as.factor(stats.part.5$HSP)
stats.part.5$nucl.mean.part.5 <- rowMeans(stats.part.5[,c("SPn", "SNn")], na.rm=TRUE)
stats.part.5$exon.mean.part.5 <- rowMeans(stats.part.5[,c("SPe", "SNe")], na.rm=TRUE)
stats.part.5$intron.mean.part.5 <- rowMeans(stats.part.5[,c("SPi", "SNi")], na.rm=TRUE)

ggplot(stats.part.5, aes(x = hsp.part.5, y =nucl.mean.part.5, color=exon_w.part.5))+
  geom_point()
# exon_W= -4 or -3.5 and hsp = 0.8
ggplot(stats.part.5, aes(x = hsp.part.5, y =exon.mean.part.5, color=exon_w.part.5))+
  geom_point()
# exon_w = -4.5 and hsp = 0.1
ggplot(stats.part.5, aes(x = hsp.part.5, y =intron.mean.part.5, color=exon_w.part.5))+
  geom_point()
# exon_w = -4 and hsp = 0.2

### TESTING ALL DATASETS IN ONE
total <- rbind(stats, stats.part.2, stats.part.3, stats.part.4, stats.part.5)
total$nucl.mean.total <- rowMeans(total[,c("SPn", "SNn")], na.rm=TRUE)
total$exon.mean.total <- rowMeans(total[,c("SPe", "SNe")], na.rm=TRUE)
total$intron.mean.total <- rowMeans(total[,c("SPi", "SNi")], na.rm=TRUE)
total$all.mean.total <- rowMeans(total[,c("nucl.mean.total", "exon.mean.total", "intron.mean.total")], na.rm = TRUE)
exon_w.total <- as.factor(total$exon_weight)
hsp.total <- as.factor(total$HSP)
site_f.total <- as.factor(total$site_factor)
exon_f.total <- as.factor(total$exon_factor)
score.total <- as.factor(total$score)

ggplot(total, aes(x = exon_w.total, y = all.mean.total, color = hsp.total))+
  geom_point()+
  facet_wrap(~ site_factor+ exon_factor)

subtotal <- rbind(stats, stats.part.2, stats.part.5)
subtotal$nucl.mean.subtotal <- rowMeans(subtotal[,c("SPn", "SNn")], na.rm=TRUE)
subtotal$exon.mean.subtotal <- rowMeans(subtotal[,c("SPe", "SNe")], na.rm=TRUE)
subtotal$intron.mean.subtotal <- rowMeans(subtotal[,c("SPi", "SNi")], na.rm=TRUE)
subtotal$all.mean.subtotal <- rowMeans(subtotal[,c("nucl.mean.subtotal", "exon.mean.subtotal", "intron.mean.subtotal")], na.rm = TRUE)
exon_w.subtotal <- as.factor(subtotal$exon_weight)
hsp <- as.factor(subtotal$HSP)
site_f.subtotal <- as.factor(subtotal$site_factor)
exon_f.subtotal <- as.factor(subtotal$exon_factor)
score.subtotal <- as.factor(subtotal$score)

ggplot(subtotal, aes(x = exon_f.subtotal, y = all.mean.subtotal, color = hsp))+
  geom_point()+
  #geom_label(data = . %>% group_by(site_factor) %>% filter(y == max(all.mean.subtotal)), aes(label = sprintf('%0.2f', y)), hjust = -0.5)+
  theme_bw()+
  facet_wrap(~ site_factor+ exon_weight, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #theme(element_text(co))
  xlab("Exon factor")+ ylab("Average Sn+Sp at all levels") + ggtitle("Anopheles gambiae: sensitivity and specificity at nucleotide+exon+intron level")

### TESTING INTRON OPTIMIZATION -> TOWARDS THE BEST PARAMETER
intron.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Anopheles_gambiae.AgamP4.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.intron_optimization.header", header = TRUE)
intron.df$nucl.mean.intron.df <- rowMeans(intron.df[,c("SPn", "SNn")], na.rm=TRUE)
intron.df$exon.mean.intron.df <- rowMeans(intron.df[,c("SPe", "SNe")], na.rm=TRUE)
intron.df$intron.mean.intron.df <- rowMeans(intron.df[,c("SPi", "SNi")], na.rm=TRUE)
intron.df$all.mean.intron.df <- rowMeans(intron.df[,c("nucl.mean.intron.df", "exon.mean.intron.df", "intron.mean.intron.df")], na.rm = TRUE)

exon_w.intron.df <- as.factor(intron.df$exon_weight)
site_f.intron.df <- as.factor(intron.df$site_factor)
exon_f.intron.df <- as.factor(intron.df$exon_factor)

ggplot(intron.df, aes(x = exon_w.intron.df, y =intron.mean.intron.df, color=exon_f.intron.df))+
  geom_point()+
  facet_wrap(~intron.df$site_f)
  
##############################
###          APIS         ###
##############################
apis.df <- read.table("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/pipeline_benchmarking/Apis_mellifera.AmelHAv.dna_rm.toplevel.hsp.exonf.sitef.score.exonw.header", header = TRUE)
apis.df <- apis.df[-c(49, 50, 69, 70, 73, 74, 97, 100, 103, 106),]
apis.df$nucl.mean.apis.df <- rowMeans(apis.df[,c("SPn", "SNn")], na.rm=TRUE)
apis.df$exon.mean.apis.df <- rowMeans(apis.df[,c("SPe", "SNe")], na.rm=TRUE)
apis.df$intron.mean.apis.df <- rowMeans(apis.df[,c("SPi", "SNi")], na.rm=TRUE)
apis.df$all.mean.apis.df <- rowMeans(apis.df[,c("nucl.mean.apis.df", "exon.mean.apis.df", "intron.mean.apis.df")], na.rm = TRUE)

exon_w.apis.df <- as.factor(apis.df$exon_weight)
site_f.apis.df <- as.factor(apis.df$site_factor)
exon_f.apis.df <- as.factor(apis.df$exon_factor)
hsp <- as.factor(apis.df$HSP)


ggplot(apis.df, aes(x = exon_f.apis.df, y =all.mean.apis.df, color=hsp))+
  geom_point()+
  theme_bw()+
  facet_grid(~exon_weight + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Exon factor")+ ylab("Average Sn+Sp at all levels") + ggtitle("Apis mellifera: sensitivity and specificity at nucleotide+exon+intron level")
# site_f = 0.5 exon_f = 0.25 exon_w = -4.25 hsp = 0.9

