setwd("/home/jaume/Escritorio/Internship_CRG/training/invertebrates/arthropods/foreval_stats")
library("ggplot2")
library("dplyr")
library("ggfortify")
library("gridExtra")
library("AICcmodavg")
##########
# PART 1 #
##########
# LOAD HEADER FILES IN FOLDER
temp = list.files(pattern="*.header")
#print(list.files(pattern="*.header"))
#for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header = TRUE, sep = '\t'))
# READ AND PROCESS MEANS IN EACH FILE, PICK ONLY THE TOP AVERAGE AT NUCLEOTIDE+EXON+INTRON LEVEL
read_and_process_all <- function(filename){
  data <- read.csv(filename, header = TRUE, sep = '\t')
  data$nucl.mean <- rowMeans(data[,c("SPn", "SNn")], na.rm=TRUE)
  data$exon.mean <- rowMeans(data[,c("SPe", "SNe")], na.rm=TRUE)
  data$intron.mean <- rowMeans(data[,c("SPi", "SNi")], na.rm=TRUE)
  data$all.mean <- rowMeans(data[,c("nucl.mean", "exon.mean", "intron.mean")], na.rm = TRUE)
  data = data[order(data$all.mean, decreasing = TRUE),]
}
# CREATE DATAFRAME WITH TOP RESULT BY SPECIES
all.params.df <- lapply(temp, read_and_process_all)
all.params.df <- do.call(rbind.data.frame, all.params.df)
drop <- c("program", "SNi", "SPi", "SNichain", "SPichain", "SNt", "SPt", "SNg", "SPg", "Matching_intron_chains", "Missed_introns", "Reference_introns", "Novel_introns", "Predicted_introns", "Missed_loci", "Reference_loci", "Novel_loci", "Predicted_loci")
all.params.df <- all.params.df[,!names(all.params.df) %in% drop]
# RETRIEVE EXCEL WITH SPECIES INFO
selected.arthro <- read.csv2("../stats_arthropods_genomes.tsv", sep = "\t", header = TRUE)
names(selected.arthro)[names(selected.arthro) == "Species"] <- "species"
# MERGE BOTH DATAFRAMES IN ONE
arthro.all.params <- merge(selected.arthro, all.params.df, by = "species")
# REMOVE ROWS WITH ZERO PREDICTED EXONS
arthro.all.params <- arthro.all.params[arthro.all.params$Predicted_exons != 0, ]
# ALL PARAMS TOGETHER
arthro.all.params <- within(arthro.all.params,  stacked <- paste(HSP, site_factor, exon_factor, exon_weight, sep="_"))
##########
# PART 2 #
##########
# Arthropods: sensitivity and specificity at nucleotide level
nucl.avg.arthro.all <- ggplot(arthro.all.params, aes(x =  exon_factor, y =nucl.mean, color=as.factor(HSP)))+
  geom_point(size = 2)+
  theme_bw()+
  #facet_grid(~exon_factor + site_factor, labeller = label_both)+
  facet_wrap(~exon_factor + site_factor, labeller = label_both)+
  theme(strip.background = element_rect(color="black", fill="#FC4E07", size=1, linetype="solid"), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(vjust = 0.5, hjust=1))+
  xlab("Exon weight")+ ylab("Average Sn+Sp at nucleotide level") + ggtitle("Arthropods: sensitivity and specificity at nucleotide level")
# DIFFERENTIAL ACCURACY BY GC CONTENT
arthro.all.params <- arthro.all.params%>%mutate(GC_group = case_when(
  GC_content<=mean(arthro.all.params$GC_content) ~ 'low_GC',
  GC_content>mean(arthro.all.params$GC_content) ~ 'high_GC'
))
p <- ggplot(arthro.all.params, aes(y =nucl.mean, color=GC_group))+
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
p1 <- ggplot(arthro.all.params, aes(y =exon.mean, color=GC_group))+
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
grid.arrange(p, p1, ncol=1, nrow = 2)
gc.nucl.aov <- aov(nucl.mean ~ GC_group, data = arthro.all.params)
summary.aov(gc.nucl.aov)
pairwise.t.test(arthro.all.params$nucl.mean, arthro.all.params$GC_group, p.adjust="none", alternative=c("greater"))
gc.exon.aov <- aov(exon.mean ~ GC_group, data = arthro.all.params)
summary.aov(gc.exon.aov)
# DIFFERENTIAL ACCURACY BY GENOME SIZE
arthro.all.params <- arthro.all.params%>%mutate(Genome_size_group = case_when(
  Genome_size<=mean(arthro.all.params$Genome_size) ~ 'low_Gsize',
  Genome_size>mean(arthro.all.params$Genome_size) ~ 'high_Gsize'
))
h <- ggplot(arthro.all.params, aes(y =nucl.mean, color=Genome_size_group))+
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
h1 <- ggplot(arthro.all.params, aes(y =exon.mean, color=Genome_size_group))+
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
grid.arrange(h, h1, ncol=1, nrow = 2)
gs.nucl.aov <- aov(nucl.mean ~ Genome_size_group, data = arthro.all.params)
summary.aov(gs.nucl.aov)
pairwise.t.test(arthro.all.params$nucl.mean, arthro.all.params$Genome_size_group, p.adjust="none", alternative=c("greater"))
gs.exon.aov <- aov(exon.mean ~ Genome_size_group, data = arthro.all.params)
summary.aov(gs.exon.aov)
# MEAN AT NUCLEOTIDE LEVEL BY PARAMETER SET
params.selected.id <- c(arthro.all.params$stacked)
nucl.mean.selected.param <- (arthro.all.params$nucl.mean)
df.params.selected.nucl <- data.frame(params.selected.id, nucl.mean.selected.param)
df.params.selected.nucl <- aggregate(.~params.selected.id, data = df.params.selected.nucl, mean)
best5.nucl <- df.params.selected.nucl[order(df.params.selected.nucl$nucl.mean, decreasing = TRUE),][1:10,]
df.params.selected.nucl <- df.params.selected.nucl  %>%
  mutate(max_nucl=ifelse(df.params.selected.nucl$nucl.mean<min(best5.nucl$nucl.mean.selected.param), "1", "0"))
# PLOT MEAN AT NUCLEOTIDE LEVEL BY PARAMETER SET - HIGHLIGHTED TOP 10
ggplot(df.params.selected.nucl)+
  geom_bar(aes(x = params.selected.id, y = nucl.mean.selected.param, fill=max_nucl), stat = "identity")+
  theme_bw()+
  scale_fill_manual( values = c( "0"="lightgreen", "1"="darkgray" ), guide = "none" )+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Average at nucleotide level by parameter set")
# MEAN AT EXON LEVEL BY PARAMETER SET
exon.mean.selected.param <- (arthro.all.params$exon.mean)
df.params.selected.exon <- data.frame(params.selected.id, exon.mean.selected.param)
df.params.selected.exon <- aggregate(.~params.selected.id, data = df.params.selected.exon, mean)
best5.exon <- df.params.selected.exon[order(df.params.selected.exon$exon.mean, decreasing = TRUE),][1:10,]
df.params.selected.exon <- df.params.selected.exon  %>%
  mutate(max_exon=ifelse(df.params.selected.exon$exon.mean<min(best5.exon$exon.mean.selected.param), "1", "0"))
# PLOT MEAN AT EXON LEVEL BY PARAMETER SET - HIGHLIGHTED TOP 5
ggplot(df.params.selected.exon)+
  geom_bar(aes(x = params.selected.id, y = exon.mean.selected.param, fill = max_exon), stat = "identity")+
  theme_bw()+
  scale_fill_manual( values = c( "0"="lightgreen", "1"="darkgray" ), guide = "none" )+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        title = element_text(size = 30),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Average at exon level by parameter set")
# PLOT EXON AND NUCLEOTIDE AVERAGE BY PARAMETER SET AND COUNTS OF EACH PARAMETER SET 
exon.avg.param.set <- ggplot(arthro.all.params)+
  geom_bar(aes(x=stacked), width = 0.5, colour="white", fill="black")+
  geom_point(data = df.params.selected.exon, aes(x = params.selected.id, y = exon.mean.selected.param), stat = "identity", colour = "orange")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1))
nucl.avg.param.set <- ggplot(arthro.all.params)+
  geom_bar(aes(x=stacked), width = 0.5, colour="white", fill="black")+
  geom_point(data = df.params.selected.nucl, aes(x = params.selected.id, y = nucl.mean.selected.param), stat = "identity", colour = "orange")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1))
# PLOT EXON AND NUCLEOTIDE AVERAGE BY PARAMETER SET AND COUNTS OF EACH PARAMETER SET -- MORE PRETTY
exon.test <- tibble(
  params.selected.id <- c(df.params.selected.exon$params.selected.id),
  exon.mean.value <- c(df.params.selected.exon$df.params.selected.nucl),
  t <- c(t)
)
exon.avg.param.set.2 <- ggplot(arthro.all.params) + 
  geom_bar(mapping = aes(x = stacked), fill = "black") +
  geom_point(data = df.params.selected.exon, aes(x = params.selected.id, y = exon.mean.selected.param), stat = "identity", colour = "blue")+
  #scale_x_date(name = "Parameters set") +
  theme_bw()+
  scale_y_continuous(name = "Parameters occurrences", 
                     sec.axis = sec_axis(~./1, name = "Average sp+sn exon level")) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "blue"),
    axis.text.y.right = element_text(color = "blue"),
    axis.ticks.y.right = element_line(color = "blue"),
    axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1))
nucl.test <- tibble(
  params.selected.id <- c(df.params.selected.nucl$params.selected.id),
  nucl.mean.value <- c(df.params.selected.nucl$df.params.selected.nucl),
  t <- c(t)
)
nucl.avg.param.set.2 <- ggplot(arthro.all.params) + 
  geom_bar(mapping = aes(x = stacked), fill = "black") +
  geom_point(data = df.params.selected.nucl, aes(x = params.selected.id, y = nucl.mean.selected.param), stat = "identity", colour = "blue")+
  #scale_x_date(name = "Parameters set") +
  theme_bw()+
  scale_y_continuous(name = "Parameters occurrences", 
                     sec.axis = sec_axis(~./1, name = "Average sp+sn nucl level")) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "blue"),
    axis.text.y.right = element_text(color = "blue"),
    axis.ticks.y.right = element_line(color = "blue"),
    axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1))

# SELECT BEST 40 PARAM SETS INDEPENDENT OF WHICH SPECIES THEY ARE FROM
top.40.params <- arthro.all.params[order(arthro.all.params$all.mean, decreasing = TRUE),][1:40,]
table(top.40.params$species) # 5 species represented
table(top.40.params$stacked) # 20 parameters sets represented
hsp <- as.factor(top.40.params$HSP)
exon_f <- as.factor(top.40.params$exon_factor)
site_f <- as.factor(top.40.params$site_factor)
exon_w <- as.factor(top.40.params$exon_weight)
# HISTOGRAM COUNTS OF EACH PARAMETER
par(mfrow=c(2,2))
hist(x = top.40.params$exon_factor, main = "Histogram top 40 params exon factor", xlab = "Top 40 exon factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = top.40.params$exon_weight, main = "Histogram top 40 params exon weight", xlab = "Top 40 exon weight", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = top.40.params$site_factor, main = "Histogram top 40 params site factor", xlab = "Top 40 site factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(x = top.40.params$HSP, main = "Histogram top 40 params HSP factor", xlab = "Top 40 HSP factor", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))
# TABLE WITH COUNTS OF EACH PARAMETER SET (makes more sense after selecting best results, here it will be ~equal because I take all the data)
k <- table(arthro.all.params[,c("stacked")])
k <- sort(k, decreasing = TRUE)
plot(table(arthro.all.params[,c("stacked")]))
# SAME AS PREVIOUS BUT WITH GGPLOT
k.plot <- ggplot(arthro.all.params)+
  geom_bar(aes(y=stacked), width = 0.5, colour="white", fill="black")+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
  ylab("Parameter combinations")
k.plot.top.40 <- ggplot(top.40.params)+
  geom_bar(aes(y=stacked), width = 0.5, colour="white", fill="black")+
  theme_bw()+
  ylab("Parameter combinations")+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  ggtitle("Count top 40 sets of parameters")+
  ylab("Parameter combinations")

# STATISTICAL TESTS AT NUCLEOTIDE LEVEL
one.way <- aov(nucl.mean ~ stacked, data = arthro.all.params)
summary.aov(one.way) # <2e-16 ***
anova.hsp <- aov(nucl.mean ~ HSP, data = arthro.all.params)
summary.aov(anova.hsp) # 0.0556 .
anova.exonf <- aov(nucl.mean ~ exon_factor, data = arthro.all.params)
summary.aov(anova.exonf) # <2e-16 ***
anova.exonw <- aov(nucl.mean ~ exon_weight, data = arthro.all.params)
summary.aov(anova.exonw) # <2e-16 ***
anova.sitef <- aov(nucl.mean ~ site_factor, data = arthro.all.params)
summary.aov(anova.sitef) # 3.98e-16 ***
anova.all <- aov(formula = nucl.mean ~ as.factor(HSP)+as.factor(exon_factor)+as.factor(exon_weight)+as.factor(site_factor), data = arthro.all.params)
summary(anova.all)
interaction.anova.all <- aov(nucl.mean ~ HSP*exon_factor*exon_weight*site_factor, data = arthro.all.params)
summary(interaction.anova.all)
# Find the best-fit model
model.set <- list(one.way, anova.hsp, anova.exonf, anova.exonw, anova.sitef, interaction.anova.all)
model.names <- c("one.way", "anova.hsp", "anova.exonf", "anova.exonw", "anova.sitef", "interaction.anova.all")
aictab(model.set, modnames = model.names)

#plot(interaction.anova.all)

tukey.anova.all<-TukeyHSD(anova.all)
tukey.anova.all

tukey.plot.aov<-aov(nucl.mean ~ as.factor(exon_factor):as.factor(exon_weight):as.factor(site_factor), data=arthro.all.params)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

mean.arthro.all.params.nucl <- arthro.all.params %>%
  group_by(HSP, exon_factor, exon_weight, site_factor) %>%
  summarise(
    nucl = mean(nucl.mean)
  )
mean.arthro.all.params.nucl <- within(mean.arthro.all.params.nucl,  stacked <- paste(HSP, site_factor, exon_factor, exon_weight, sep="_"))
ggplot(mean.arthro.all.params.nucl, aes(x = exon_factor, y = nucl, group=stacked)) +
  theme_classic() +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~ exon_factor)

ggplot(mean.arthro.all.params.nucl, aes(x = site_factor, y = nucl, group=stacked)) +
  theme_classic()+
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~ site_factor)

ggplot(mean.arthro.all.params.nucl, aes(x = HSP, y = nucl, group=stacked)) +
  theme_classic()+
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~ HSP)

ggplot(mean.arthro.all.params.nucl, aes(x = exon_weight, y = nucl, group=stacked)) +
  theme_classic()+
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~ exon_weight)
random <- arthro.all.params[1:2,]

# STATISTICAL TEST OF ALL RUNS BY SPECIES AT NUCLEOTIDE LEVEL
arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',]
one.way <- aov(nucl.mean ~ stacked, data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary.aov(one.way)
anova.hsp <- aov(nucl.mean ~ HSP, data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary.aov(anova.hsp)
anova.exonf <- aov(nucl.mean ~ exon_factor, data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary.aov(anova.exonf)
anova.exonw <- aov(nucl.mean ~ exon_weight, data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary.aov(anova.exonw)
anova.sitef <- aov(nucl.mean ~ site_factor, data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary.aov(anova.sitef)
anova.all <- aov(formula = nucl.mean ~ as.factor(HSP)+as.factor(exon_factor)+as.factor(exon_weight)+as.factor(site_factor), data = arthro.all.params[arthro.all.params$species == 'Ancistrocerus_nigricornis',])
summary(anova.all)
interaction.anova.all <- aov(nucl.mean ~ HSP*exon_factor*exon_weight*site_factor, data = arthro.all.params)
summary(interaction.anova.all)
# Find the best-fit model
model.set <- list(one.way, anova.hsp, anova.exonf, anova.exonw, anova.sitef, interaction.anova.all)
model.names <- c("one.way", "anova.hsp", "anova.exonf", "anova.exonw", "anova.sitef", "interaction.anova.all")
aictab(model.set, modnames = model.names)
# ANOVA NUCLEOTIDE BY SPECIES
u <- unique(arthro.all.params$species)
for (sp in u) {
  aov.by.species <- aov(nucl.mean ~ HSP*exon_factor*exon_weight*site_factor, data = arthro.all.params[arthro.all.params$species == sp,])
  sum.aov <- summary(aov.by.species)
}

# FIT MODEL NUCLEOTIDE LEVEL
fit.model.nucl <- glm(data = arthro.all.params, formula = nucl.mean ~ species + stacked)
fit.model.nucl$coefficients
fit.model.nucl$aic
coef.params <- fit.model.nucl$coefficients[c(55:134)]
coef.species <- fit.model.nucl$coefficients[c(1:54)]

tail(sort(coef.params))
# FIT MODEL EXON LEVEL
fit.model.exon <- glm(data = arthro.all.params, formula = exon.mean ~ species)
fit.model.exon$aic
species.coef.exon <- fit.model.exon$coefficients
species.coef.exon
df.fit.model <- data.frame(species.coef.nucl, species.coef.exon)
boxplot(df.fit.model)
ggplot(df.fit.model)+
  geom_boxplot()
