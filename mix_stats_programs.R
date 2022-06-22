library("ggplot2")
library("tidyr")
library("gridExtra")
library("dplyr")
library("scales")
augustus <- data.frame(species  = c("Anopheles_gambiae", "Apis_mellifera",  "Bombus_terrestris", "C_elegans", "Drosophila_melanogaster", "Bemisia_tabaci", "Drosophila_erecta", "Drosophia_mojavensis", "Necator_americanus", "Aedes_aegypti", "Glossina_austeni", "Solenopsis_invicta", "Anopheles_funestus", "Danaus_plexippus"),
                       program = c("Augustus", "Augustus", "Augustus", "Augustus","Augustus", "Augustus","Augustus", "Augustus","Augustus", "Augustus","Augustus", "Augustus","Augustus", "Augustus"),
                       SNn = c(89.1, 29.0, 39.0, 89.9, 91.4, 15.7, 92.8, 87.9, 56.8, 11.1, 40.1, 25.8, 88.9, 53.1),
                       SPn = c(55.8, 42.1, 34.4, 73.5, 72.3, 31,   70.5, 66.6, 47.1, 78.9, 70.4, 39.5, 54,   39.1),
                       SNe = c(48.0, 12.5, 17.8, 66.2, 55.1, 8.5,  54,   47.8, 34.6, 8.4,  5.8,  12.4, 47.6, 33.0),
                       SPe = c(39.3, 36.6, 31.1, 67.1, 53.3, 24.5, 51.8, 48.6, 41.8, 56.6, 33.5, 35.9, 35.4, 46.2),
                       SNi = c(72.0, 13.6, 19.5, 77.5, 80.3, 7.5,  79.1, 69.9, 40.2, 8.4,  8.5,  14.2, 71.5, 34.6),
                       SPi = c(53.2, 39.8, 34.1, 80.2, 72.5, 40.7, 71.3, 66.9, 51.9, 60.1, 43.8, 40.6, 48.2, 47.7),
                       Predicted_exons = c(65588, 23949, 40803, 119119, 56622, 28964, 56616, 54161, 100103, 9295, 14985, 40640, 71825, 68925),
                       Predicted_loci = c(10495, 3337, 6057, 21310, 11056, 3343, 11192, 10453, 40655, 2615, 4050, 4911, 11627, 8401),
                       time_user = c(80414, 85124, 64436, 34039, 45211, 54755, 48506, 71238, 127410, 195615, 95010, 140503, 101015, 42713),
                       genome_size = c(1195047198, 283349351, 225250884, 407584042, 248654244, 100286401, 143726002, 370264922, 244075060, 378101515, 193826310, 152712140, 444543995, 248676414))

genemark <- data.frame(species  = c("Aedes_aegypti", "Anopheles_gambiae", "Apis_mellifera", "Bemisia_tabaci", "Bombus_terristris", "C_elegans", "Drosophila_melanogaster", "Glossina_austeni", "Necator_americanus", "Solenopsis_invicta", "Drosophia_mojavensis", "Drosophila_erecta", "Anopheles_funestus", "Danaus_plexippus"),
                       program = c("GeneMark", "GeneMark", "GeneMark", "GeneMark","GeneMark", "GeneMark","GeneMark", "GeneMark","GeneMark", "GeneMark","GeneMark", "GeneMark","GeneMark", "GeneMark"),
                       SNn = c(81.6, 77.6, 72.6, 19.5, 84.6, 89.7, 85.1, 76.1, 76.4, 64.8, 82.5, 86.0, 86.7, 92.1),
                       SPn = c(72.3, 87.9, 92.8, 21.3, 80.6, 90.2, 95.4, 87.9, 73.7, 31.5, 89.1, 94.7, 52.0, 66.1),
                       SNe = c(54.5, 51.9, 33.1, 9.9,  60.3, 76.7, 66.6, 33.2, 60.4, 36.1, 56.5, 65.2, 62.9, 75.3),
                       SPe = c(38.2, 47.1, 32.7, 9.6,  52.5, 73.3, 63.4, 40.6, 57.8, 12.5, 48.5, 61.3, 31.0, 49.6),
                       SNi = c(52.5, 59.1, 39.6, 9.0,  65.8, 82.2, 73.2, 40.1, 65.5, 39.3, 65.7, 72.3, 69.6, 76.0),
                       SPi = c(40.4, 52.0, 38.7, 9.3,  58.9, 77.0, 67.3, 47.4, 59.2, 14.0, 54.0, 65.7, 34.1, 52.0),
                       Predicted_exons = c(89031, 59114, 70983, 85905, 81780, 126196, 57650, 70406, 126334, 254990, 64177, 57854, 108381, 146547),
                       Predicted_loci = c(27090,  12788, 9285, 17437, 13974, 18280, 12924, 14049, 12862, 48515, 13298, 12751, 25492, 243347),
                       time_user = c(21022, 10085, 6018, 9906, 12119, 3577, 6341, 9845, 9340, 11769, 7269, 18637, 13739, 49444)*4,
                       genome_size = c(1195047198, 283349351, 225250884, 407584042, 248654244, 100286401, 143726002, 370264922, 244075060, 378101515, 193826310, 152712140, 444543995, 248676414))

geneid <- data.frame(species = c("Aedes_aegypti", "Anopheles_funestus", "Anopheles_gambiae", "Apis_mellifera", "Bemisia_tabaci", "Bombus_terristris", "C_elegans", "Drosophila_erecta", "Drosophia_mojavensis", "Drosophila_melanogaster", "Glossina_austeni", "Necator_americanus", "Solenopsis_invicta", "Danaus_plexippus"),
                       program = c("Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid", "Geneid"),
                       SNn = c(69.5, 71.3, 75.8, 63.1, 24.7, 76.5, 80.5, 80.8, 71.4, 79.4, 6.10, 57.5, 71.7, 35.8),
                       SPn = c(86.9, 53.1, 83.7, 84.9, 27.1, 74.2, 89.1, 94.2, 91.9, 95.7, 65.9, 55.7, 52.1, 78.1),
                       SNe = c(52.0, 48.1, 54.7, 31.9, 17.3, 54.6, 63.7, 61.2, 50.4, 62.2, 3.70, 35.3, 47.0, 33.1),
                       SPe = c(44.0, 24.5, 43.1, 34.3, 14.5, 40.7, 67.8, 56.2, 48.8, 61.7, 10.5, 34.2, 29.4, 43.4),
                       SNi = c(51.3, 50.5, 58.3, 32.3, 16.3, 54.5, 63.6, 65.2, 54.1, 65.8, 2.10, 31.7, 47.4, 28.2),
                       SPi = c(43.1, 26.2, 44.9, 40.6, 18.9, 48.3, 70.3, 62.6, 55.3, 64.3, 5.90, 38.7, 37.9, 37.1),
                       Predicted_exons = c(73775, 104675, 68164, 65241, 99189, 95595, 113365, 59400, 56874, 55244, 40754, 124604, 141470, 73818),
                       Predicted_loci = c(16957,  264901, 15325, 17365, 37919, 27085, 21761, 16532, 15961, 13135, 7618, 40690, 49402, 10237),
                       time_user = c(881, 476, 347, 228, 297, 324, 104, 179, 401, 222, 374, 216, 413, 314),
                       genome_size = c(1195047198, 283349351, 225250884, 407584042, 248654244, 100286401, 143726002, 370264922, 244075060, 378101515, 193826310, 152712140, 444543995, 248676414))

geneid.BLASTx <- data.frame(species = c("Anopheles_gambiae", "Apis_mellifera", "Bombus_terristris", "C_elegans", "Drosophila_melanogaster", "Solenopsis_invicta", "Bemisia_tabaci", "Anopheles_funestus", "Glossina_austeni", "Danaus_plexippus", "Necator_americanus", "Aedes_aegypti", "Drosophila_erecta", "Drosophia_mojavensis"),
                            program = c("Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx", "Geneid+BLASTx"),
                            SNn = c(  68.7, 54.1,  65.1, 82.9, 78.2, 59.5, 55.8, 76,    53.6, 87.5, 53.0, 73.7, 86.7, 72.7),
                            SPn = c(  94.7, 92.4,  88.8, 89.1, 96,   74.3, 60.5, 91.32, 89.1, 77.5, 74.5, 88.3, 84.8, 92.4),
                            SNe = c(  46.5, 25.6,  40.9, 64.9, 56.4, 36.6, 35.6, 49.9,  19.8, 63.5, 34.9, 47.9, 60.9, 50.9),
                            SPe = c(  55.1, 40.9,  56.5, 62.0, 63.1, 47.6, 39.0, 55.1,  41.2, 61.7, 41.7, 57.7, 63.5, 52.4),
                            SNi = c(  49.9, 25.7,  41.7, 71.1, 61.9, 36.3, 33.9, 53.4,  22.1, 62.8, 34.5, 46.0, 64.2, 56.5),
                            SPi = c(  58.1, 44.6,  60.9, 67.1, 67.6, 52.5, 40.5, 57.7,  46.1, 63.6, 41.6, 56.2, 68.2, 56.8),
                            Predicted_exons = c(  45254, 43909, 51546, 126353, 48945,   67885, 75903, 48351, 41432, 25523, 101069, 51826, 52135, 53587),
                            Predicted_loci = c(  10280, 9254,  10002, 19058,   11403,   16847, 16554, 10800,  9697, 8610, 15825, 12779, 13582, 11976),
                            time_user = c(103,193, 198, 146, 160, 315, 460, 263, 288, 233, 709, 835, 260, 400)*4,
                            genome_size = c(225250884, 407584042, 100286401, 143726002, 378101515, 444543995, 248654244, 283349351, 193826310, 248676414, 152712140, 1195047198, 370264922, 244075060))


programs.mix <- rbind(augustus, genemark, geneid, geneid.BLASTx)
#cols.name <- c("SNn", "SPn", "SNe", "SPe", "SNi", "SPi", "Predicted_exons", "Predicted_loci", "time_user", "genome_size")
#programs.mix[cols.name] <- sapply(programs.mix[cols.name],as.numeric)

programs.mix$nucl.mean.programs.mix <- rowMeans(programs.mix[,c("SPn", "SNn")], na.rm=TRUE)
programs.mix$exon.mean.programs.mix <- rowMeans(programs.mix[,c("SPe", "SNe")], na.rm=TRUE)
programs.mix$intron.mean.programs.mix <- rowMeans(programs.mix[,c("SPi", "SNi")], na.rm=TRUE)
programs.mix$all.mean.programs.mix <- rowMeans(programs.mix[,c("nucl.mean.programs.mix", "exon.mean.programs.mix", "intron.mean.programs.mix")], na.rm = TRUE)

# NUCLEOTIDE LEVEL
# SENSITIVITY
nu.se <- ggplot(programs.mix, aes(x=program, y=SNn, fill=program), show.legend = FALSE) +
  #geom_violin(width=1.4, alpha=0.2) +
  geom_boxplot(show.legend = FALSE)+
  geom_jitter(size=2, alpha=2.2, aes(colour=program), show.legend = FALSE) +
  theme(plot.title = element_text(size=40), 
        #legend.position="none",
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        #legend.text = element_text(size = 40),
        #legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  ggtitle("Nucleotide sensitivity") + xlab("") + ylab("Accuracy")
nu.se
# NUCLEOTIDE LEVEL
# SPECIFICITY
nu.sp <- ggplot(programs.mix, aes(x=program, y=SPn, fill=program), show.legend = FALSE) +
  #geom_violin(width=1.4, alpha=0.2) +
  geom_boxplot(show.legend = FALSE)+
  geom_jitter(size=2, alpha=2.2, aes(colour=program), show.legend = FALSE) +
  theme(legend.position="none", plot.title = element_text(size=40), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  ggtitle("Nucleotide specificity") + xlab("") + ylab("Accuracy")
nu.sp
#grid.arrange(nu.se, nu.sp, nrow =1, ncol = 2)
# NUCLEOTIDE LEVEL
# AVERAGE SENSITIVITY+SPECIFICITY
nu.avg <- ggplot(programs.mix, aes(x=program, y=nucl.mean.programs.mix, fill=program)) +
  #geom_violin(width=1.4, alpha=0.2) +
  geom_boxplot()+
  geom_jitter(size=2, alpha=2.2, aes(colour=program)) +
  theme(legend.position="none", plot.title = element_text(size=40), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  ggtitle("Nucleotide average") + xlab("") + ylab("Accuracy")
nu.avg
grid.arrange(nu.se, nu.sp, nu.avg, nrow =1, ncol = 3)
mix.nucl.aov <- aov(nucl.mean.programs.mix ~ program, data = programs.mix)
summary.aov(mix.nucl.aov)
pairwise.t.test(programs.mix$nucl.mean.programs.mix, programs.mix$program, p.adjust="none", alternative=c("greater"))
pairwise.t.test(programs.mix$SNn, programs.mix$program, p.adjust="none", alternative=c("greater"))
pairwise.t.test(programs.mix$SPn, programs.mix$program, p.adjust="none", alternative=c("greater"))
# EXON LEVEL
# SENSITIVITY
ex.se <- ggplot(programs.mix, aes(x=program, y=SNe, fill=program), show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE) +
  #geom_violin(width=1.4, alpha=0.2) +  
  geom_jitter(size=0.6, alpha=0.9, aes(colour=program), show.legend = FALSE) +
  theme(legend.position="none", plot.title = element_text(size=40), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  ggtitle("Exon sensitivity") + xlab("") + ylab("Accuracy")
ex.se
# EXON LEVEL
# SPECIFICITY
ex.sp <- ggplot(programs.mix, aes(x=program, y=SPe, fill=program), show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE) +
  #geom_violin(width=1.4, alpha=0.2) + 
  geom_jitter(size=0.6, alpha=0.9, aes(colour=program), show.legend = FALSE) +
  theme(legend.position="none", plot.title = element_text(size=40), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  ggtitle("Exon specificity") + xlab("") + ylab("Accuracy")
ex.sp
# EXON LEVEL
# AVERAGE SENSITIVITY+SPECIFICITY
ex.avg <- ggplot(programs.mix, aes(x=program, y=exon.mean.programs.mix, fill=program)) +
  geom_boxplot() +
  #geom_violin(width=1.4, alpha=0.2) +
  geom_jitter(size=0.6, alpha=0.9, aes(colour=program)) +
  theme(legend.position="none", plot.title = element_text(size=40), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
  theme_bw() +
  #ylim(0,100) +
  ggtitle("Exon average") + xlab("") + ylab("Accuracy")
ex.avg
grid.arrange(ex.se, ex.sp, ex.avg, nrow =1, ncol = 3)
mix.ex.aov <- aov(exon.mean.programs.mix ~ program, data = programs.mix)
summary.aov(mix.ex.aov)
pairwise.t.test(programs.mix$exon.mean.programs.mix, programs.mix$program, p.adjust="none", alternative=c("greater"))
pairwise.t.test(programs.mix$SNe, programs.mix$program, p.adjust="none", alternative=c("greater"))
pairwise.t.test(programs.mix$SPe, programs.mix$program, p.adjust="none", alternative=c("greater"))
# EXECUTION TIME
ex.time <- ggplot(programs.mix, aes(x=genome_size, y=log(time_user, base = exp(10)), group=program, color=program)) +
  geom_point()+
  geom_smooth(method = 'lm', se=F)+
  theme_bw()+
  theme(plot.title = element_text(size=40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)
        #panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        #panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
  ) +
  ggtitle("Execution time")+xlab("Genome size, bp")+ylab("log(execution time, seconds)")
ex.time
ex.time.aov <- aov(time_user ~ genome_size, data = programs.mix)
summary.aov(mix.ex.aov)
pairwise.t.test(programs.mix$time_user, programs.mix$program, p.adjust="none", alternative=c("two.sided"))
#qqnorm(programs.mix$nucl.mean.programs.mix, pch = 1, frame = FALSE)
#qqline(programs.mix$nucl.mean.programs.mix, col = "steelblue", lwd = 2)


ggplot(programs.mix, aes(x=genome_size, y=time_user, group=program, color=program)) +
  geom_point()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
              labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = c(1e8, 5e8, 10e8),
                     labels = c(expression("1·"~10^8),
                                expression("5·"~10^8),
                                expression("1·"~10^9))) +
  geom_smooth(method='lm', formula = y ~ x, alpha = 0.2) +
  theme_bw() + 
  theme(plot.title = element_text(size=40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  ylab("log(Execution time, seconds)") + 
  xlab("Genome size, bp") +
  ggtitle("Execution time")
