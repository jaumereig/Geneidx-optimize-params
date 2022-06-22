rm(list=ls())

# require(shiny)
# require(plotly)
require(reshape)
require(ggpubr)
require(viridis)
require(ggplot2)

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("../plots")



### ALL SPECIES SGP2 ###
stats_dir <- "../data"

#eval_file_name <- paste(stats_dir, "/evaluation_5_with_over200.Nov2021.tsv.header", sep = "")
#eval_file_name <- paste(stats_dir, "/evaluation_5_with_over200.Nov2021.only20.tsv.header", sep = "")
# eval_file_name <- paste(stats_dir, "/evaluation_5_with_over200.tsv.header", sep = "")
# eval_file_name <- paste(stats_dir, "/evaluation_5_with_over200.only20.tsv.header", sep = "")
# eval_file_name <- paste(stats_dir, "/evaluation_6_with_over200.tsv.header", sep = "")
eval_file_name <- paste(stats_dir, "/evaluation_8_with_over200.tsv.header", sep = "")

eval_data <- read.table(eval_file_name, header = T)
eval_data <- unique(eval_data)



num_progs_per_sp = aggregate(eval_data$program, by = list(eval_data[,c("species")]),
                             FUN = length)
colnames(num_progs_per_sp) <- c("species", "counts")
species_to_keep = num_progs_per_sp[num_progs_per_sp$counts >= 12,][,"species"]
eval_data <- eval_data[eval_data$species %in% species_to_keep,]


stats_data <- eval_data[eval_data$program == "Geneid+BLASTx_Hum",]
stats_data <- stats_data[stats_data$UNMATCHED == "over200",]
stats_data$meanN <- (stats_data$SNn + stats_data$SPn)/2
stats_data$meanE <- (stats_data$SNe + stats_data$SPe)/2

aggregate(stats_data$meanN, by = list(stats_data[,c("species")]),
          FUN = max)




all_data <- eval_data

all_data[,c("SNn", "SPn","SNe", "SPe","SNi", "SPi",
            "SNichain", "SPichain",
            "SNt", "SPt","SNg", "SPg")] = all_data[,c("SNn", "SPn",
                                                      "SNe", "SPe","SNi","SPi",
                                                      "SNichain", "SPichain",
                                                      "SNt", "SPt",
                                                      "SNg", "SPg")]/100

# new_column_names <- c("species", "program", "UNMATCHED",
#                       "SN", "SP", "SN", "SP",
#                       "SN", "SP", "SN", "SP",
#                       "SN", "SP", "SN", "SP",
#                       "Matching_intron_chains",
#                       "Matching_transcripts", "Matching_loci", "Missed_exons", "Reference_exons", 
#                       "Novel_exons", "Predicted_exons", "Missed_introns", "Reference_introns", 
#                       "Novel_introns", "Predicted_introns", "Missed_loci", "Reference_loci", 
#                       "Novel_loci", "Predicted_loci")
# 
# colnames(all_data) <- new_column_names

# species <- c("H.sapiens", "M.musculus")
species <- unique(all_data$species)
selected_eval_data <- all_data[all_data$species %in% species,]
melted_eval_data <- melt(selected_eval_data, id.vars = c("species", "program"))


evaluation_data_groups <- c('nucleotide', 'nucleotide',
                            'exon', 'exon', 'intron', 'intron',
                            'intron_chain', 'intron_chain',
                            'transcript', 'transcript', 'gene', 'gene',
                            'matching', 'matching', 'matching',
                            'W_M_exons', 'W_M_exons', 'W_M_exons', 'W_M_exons',
                            'W_M_introns', 'W_M_introns', 'W_M_introns', 'W_M_introns',
                            'W_M_loci', 'W_M_loci', 'W_M_loci', 'W_M_loci')

# repetitions_per_species <- nrow(selected_eval_data[selected_eval_data$species == selected_eval_data$species[1],])
# num_of_species <- length(unique(selected_eval_data$species))
# print(paste(repetitions_per_species, num_of_species))
# repeats = repetitions_per_species * num_of_species
repeats <- nrow(selected_eval_data)
group_labels = rep(evaluation_data_groups, each = repeats)

melted_eval_data$grp <- group_labels

reduced_melted_data <- melted_eval_data[melted_eval_data$grp %in% c(
  # 'gene',
  # 'transcript',
  'exon',
  # 'intron',
  # 'intron_chain',
  'nucleotide'),]


reduced_melted_data[reduced_melted_data$value <= 0.01,]

reduced_melted_data <- reduced_melted_data[reduced_melted_data$value > 0.01,]



p <- ggplot(reduced_melted_data, aes(factor(variable), value))
# pdf("evaluation_5programs_over200_withWG.pdf", height = 6, width = 10)
# pdf("evaluation_5programs_over200_withoutWG.pdf", height = 6, width = 10)
pdf("evaluation_8programs_over200_joined_class.pdf", height = 6, width = 12)
# pdf("evaluation_7programs_over200_withWG_noREF.pdf", height = 6, width = 14)
p + geom_boxplot(aes(color = program)) +
  theme_classic() +
  # facet_wrap(~grp, scales = "free_x") +
  facet_wrap(~grp, scales = "free") +
  ylim(0,1) +
  scale_color_discrete(name = "Program") +
  xlab("Metric") + 
  ylab("Value") +
  # labs(title = paste("Evaluation of 5 gene predictions in", length(species), "vertebrate species")) + 
  labs(title = paste("Evaluation in", length(species), "vertebrate species")) +
  # labs(title = paste("Evaluation in", length(species), "vertebrate species WITHOUT WG")) +
  theme(text = element_text(size = 15))
dev.off()



p <- ggplot(reduced_melted_data, aes(factor(variable), value, color = program))
pdf("evaluation_5programs_over200_jitter.pdf", height = 6, width = 10)
# pdf("evaluation_6programs_over200_withWG_jitter.pdf", height = 6, width = 10)
# pdf("evaluation_7programs_over200_withWG_jitter.pdf", height = 6, width = 10)
p + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.05), alpha = 0.3) +
  theme_classic() +
  # facet_wrap(~grp, scales = "free_x") +
  facet_wrap(~grp, scales = "free") +
  ylim(0,1) +
  scale_color_discrete(name = "Program") +
  xlab("Metric") + 
  ylab("Value") +
  labs(title = paste("Evaluation of 5 gene predictions in", length(species), "vertebrate species")) + 
  theme(text = element_text(size = 15))
dev.off()





