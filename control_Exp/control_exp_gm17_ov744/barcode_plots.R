#!/usr/bin/env Rscript




# This script takes the multi sample / timepoint output from 'Bartender'
# It creates barplots for absolute barcode read numbers and relative barcode read abundance

# Requirements:
#
# Input:
# 1) 'Bartender' output file prefix (e.g. 'all_samples')
# 2) Threshold barcode relative abundance in percent (e.g. '2' or '1.2')
# 3) Metadata file for Samples 
#    (!!! needs to be sorted ascendingly by sample number !!!)
#    (!!! all text values in UPPER case !!!)
#    (!!! column 'Sample' may only contain integers !!!)
#    (!!! column 'Sample' needs to start with 1, 2, 3, ...; use column 'Index' for the S## number given by sequencer !!!)
# 4) Metadata file for Barcodes
# 
# 
# 
# 
# 
# 
# 

suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('dplyr')) install.packages('dplyr'))
suppressMessages(if (!require('RColorBrewer')) install.packages('RColorBrewer'))

write("\n\n----------", stdout())
write("Loading libraries", stdout())
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))

#*********************************************************************************************
#********PREPARATION**************************************************************************
#*********************************************************************************************

#read in script arguments <- this needs to be run!
args = commandArgs(trailingOnly=TRUE)


#for debugging
debug <- 1

if (debug==1) {
write("DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE", stdout())
args[1] <- "control_exp"
args[2] <- "0"
args[3] <- "sample_metadata_complete_corrected.csv"
args[4] <- "barcode_metadata.csv"
}

if (length(args)<4) {
  stop("Specify 'Bartender' output file prefix (e.g. 'all_samples'), threshold barcode relative abundance in percent (e.g. '2'), metadata file for samples, and metadata file for barcodes. 4 Parameters total.", call.=FALSE)
}

#make a results directory
write("\n\n----------", stdout())
write("Creating results directory", stdout())
dir.create("results", showWarnings = FALSE)

#create log files (if debug mode off)
if (debug==0) {
write("\n\n----------", stdout())
write("Processing; all further output in results/logs", stdout())
msglog <- file(file.path("results","message_log.txt",fsep = .Platform$file.sep), open = "a")
sink(file.path("results","analysis_log.txt",fsep = .Platform$file.sep),type='output',append=F)
sink(msglog,type='message',append=F)
}

write("\n\n----------", stdout())
write("Chosen barcode percentage cutoff:", stdout())
#threshold percent value from argument
threshold <- as.numeric(args[2])
threshold

#create filenames of input files and for output files from frefix supplied as script argument

name_clusterfile <- paste(args[1],"_cluster.csv",sep="")
name_qualityfile <- paste(args[1],"_quality.csv",sep="")
name_singlefile <- paste(args[1],"_single.csv",sep="")
name_combinedfile <- paste(args[1],"_combined.csv",sep="")

#make list of all files
table_inputs <- c(name_clusterfile,name_qualityfile)
other_inputs <- c(name_singlefile,name_combinedfile)


#read in table files (cluster and quality)
write("\n\n----------", stdout())
write("Reading in Bartender files", stdout())
for (i in 1:length(table_inputs)) { 
  nam1 <- unlist(strsplit(table_inputs[i],"_"))
  nam2 <- nam1[length(nam1)]
  assign(nam2, read.table(table_inputs[i], header=T, sep=","))
  rm(i,nam1,nam2)
}


#read in other files (single and combined)
for (i in 1:length(other_inputs)) { 
  nam1 <- unlist(strsplit(other_inputs[i],"_"))
  nam2 <- nam1[length(nam1)]
  assign(nam2, read.table(other_inputs[i], header=F, sep=","))
  rm(i,nam1,nam2)
}

#read in metadata files
write("Reading in metadata", stdout())
args[3]
args[4]
meta_samples <- read.table(args[3], header=T, sep=",")
meta_barcodes <- read.table(args[4], header=T, sep=",")

meta_samples$Sample <- meta_samples$Index

#remove naming tools
rm(name_clusterfile,name_qualityfile,name_singlefile,name_combinedfile)

#calculate number of samples in input
write("\n\n----------", stdout())
write("Counting samples", stdout())
z<-as.numeric(length(cluster.csv)-3)
z
sample_column_first <- "time_point_1"
sample_column_last <- paste("time_point",z,sep="_")

#create a list of the known/defined barcodes and genotypes
list_known_barcodes <- unique(as.character(meta_barcodes[,"Center"]))
list_known_libraries <- sort(unique(as.character(meta_barcodes[,"Library"])))

#create a dictionary that assigns a certain color to every genotype to be used in plotting
getPalette = colorRampPalette(brewer.pal(12, "Paired"))         #create color gradient from x colors of set (https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html#rcolorbrewer-package)
colors <- getPalette(12)                                        #retrieve all 12 colors included in the palette above

colors <- split(colors, ceiling(seq_along(colors)/length(list_known_libraries))) #remove colors that are not need, only leaving behind as many as there are genotypes
library_colors <- unlist(colors[1])

#name the vector of colors with the genotypes, basically generating a dictionary genotype -> color
names(library_colors) <- list_known_libraries

#add grey for unassigned
unassigned_color <- "#D3D3D3"
names(unassigned_color) <- "unassigned"
library_colors <- append(library_colors, unassigned_color, after = length(library_colors))



#*********************************************************************************************
#********ANALYSIS*****************************************************************************
#*********************************************************************************************
write("\n\n----------", stdout())
write("Data shuffling", stdout())
#ON ABSOLUTE READ BASIS

#per sample, add up all reads that contain a barcode, make long, and then separate 'time_point_#' to '#'
cluster_reads_long <- cluster.csv %>% select(Center,all_of(sample_column_first):all_of(sample_column_last)) %>% 
                          summarise_at(vars(all_of(sample_column_first):all_of(sample_column_last)), sum, na.rm = TRUE) %>% 
                          gather(key="Sample", value="Barcode",all_of(sample_column_first):all_of(sample_column_last)) %>% 
                          separate(.,Sample,c(NA,NA,"Sample")) %>% 
                          transform(., Sample = as.numeric(Sample))

#obtain details from metadata file, add to dataframe above
cluster_reads_long$TotalReads <- meta_samples$ReadsSum
cluster_reads_long$Unknown <- cluster_reads_long$TotalReads - cluster_reads_long$Barcode
cluster_reads_long$Experiment <- meta_samples$Experiment
cluster_reads_long$Medium <- meta_samples$Medium
cluster_reads_long$Substrate <- meta_samples$Substrate
cluster_reads_long$Replicate <- meta_samples$Replicate


#create a dictionary of barcode centers and metadata 
suppressWarnings(suppressMessages(
cluster_dict <- cluster.csv %>% select(Center) %>% left_join(.,meta_barcodes) %>% select(Center,Library) %>% arrange(Center)
))
#if no genotype is specified for a cluster/barcode, name it "unassigned"
cluster_dict$Library <- ifelse(is.na(cluster_dict$Library), "unassigned", as.character(cluster_dict$Library))

#write out cluster.id dictionary
write_tsv(cluster_dict, file.path("results",paste(args[1],"_cluster_dict.tab",sep="")), col_names=T)





#ON RELATIVE ABUNDANCE BASIS

#calculate relative abundance of barcodes (in percent) from absolute read numbers and add to original cluster.csv
for (i in 4:(length(cluster.csv))) {
sapply(names(cluster.csv)[i], function(x) {
  cluster.csv[paste0("Sample_", i-3, "_pct")] <<- cluster.csv[x] / sum(cluster.csv[x]) * 100
})
}
#store the column names for future reference
sample_pct_column_first <- "Sample_1_pct"
sample_pct_column_last <- paste("Sample", z, "pct", sep="_")

#create list of barcode centers that are have percentage >= threshold in ANY sample
list_cluster_pct_threshold <- as.character (cluster.csv %>% 
                                            select(Center, all_of(sample_pct_column_first):all_of(sample_pct_column_last)) %>% 
                                            gather(key="Sample", value="Percentage",all_of(sample_pct_column_first):all_of(sample_pct_column_last)) %>% 
                                            filter(Percentage >= threshold) %>%
                                            select(1) %>% arrange(Center) %>%
                                            unique() %>% t(.)
                                          )


#delete all rows that have barcode centers that are not listed in the list created above
cluster_pct_threshold <- cluster.csv %>% select(Center,all_of(sample_pct_column_first):all_of(sample_pct_column_last)) %>% filter(.[[1]] %in% list_cluster_pct_threshold)

#expand filtered data set by barcode metadata obtained from cluster_dict
#make long and separate 'sample_#_pct' into '#'
#add sample metadata
suppressWarnings(suppressMessages(
clusters_pct_threshold_meta_long <- cluster_pct_threshold %>% 
                                        left_join(.,cluster_dict) %>% 
                                        gather(key="Sample", value="Percentage",all_of(sample_pct_column_first):all_of(sample_pct_column_last)) %>% 
                                        separate(.,Sample,c(NA,"Sample",NA)) %>% 
                                        transform(., Sample = as.integer(Sample)) %>%
                                        left_join(.,meta_samples, ) 
))

clusters_pct_threshold_meta_long <- clusters_pct_threshold_meta_long %>% filter(Experiment=="GM17_OV744")
clusters_pct_threshold_meta_long$Experiment <- as.factor(clusters_pct_threshold_meta_long$Experiment)

#*********************************************************************************************
#********PLOTTING*****************************************************************************
#*********************************************************************************************
write("\n\n----------", stdout())
write("Plotting", stdout())
cluster_reads_long_barcoded <- cluster_reads_long %>% select(Sample,Barcode,Unknown,Experiment,Medium,Substrate,Replicate) %>% gather(key="BarcodePresence", value="Reads",2:3)

#cluster_reads_long_barcoded$Sample <- factor(cluster_reads_long_barcoded$Sample, levels = unique(cluster_reads_long_barcoded$Sample))
#cluster_reads_long_barcoded$BarcodePresence <- factor(cluster_reads_long_barcoded$BarcodePresence, levels = c("Unknown","Barcode"))
#cluster_reads_long_barcoded$Structure <- factor(cluster_reads_long_barcoded$Structure, levels = c("Inoculum","Primary","Secondary","Tertiary"))
#cluster_reads_long_barcoded$Type <- factor(cluster_reads_long_barcoded$Type, levels = c("Tip","Hair","Segment","Inoculum"))

#PLOT - Total read number, total barcode read number
p <- ggplot(cluster_reads_long_barcoded, aes(x=Sample, y=Reads, fill=BarcodePresence)) +
  scale_fill_manual(values = c("Grey","violetred4")) +
  geom_bar(position="stack", stat="identity") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .50)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_total_reads_absolute.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Total percent of reads that contain any barcode
p <- ggplot(cluster_reads_long_barcoded, aes(x=Sample, y=Reads, fill=BarcodePresence)) +
  scale_fill_manual(values = c("Grey","violetred4")) +
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") + 
  ylab("Reads (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_total_reads_relative.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Total read number, total barcode read number, in GRID
p <- ggplot(cluster_reads_long_barcoded, aes(x=Replicate, y=Reads, fill=BarcodePresence)) +
  scale_fill_manual(values = c("Grey","violetred4")) +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  geom_bar(position="stack", stat="identity") +
  xlab("Sample") + 
  ylab("Total reads (absolute)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_total_reads_absolute_grid.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Total percent of reads that contain any barcode, in GRID
p <- ggplot(cluster_reads_long_barcoded, aes(x=Replicate, y=Reads, fill=BarcodePresence)) +
  scale_fill_manual(values = c("Grey","violetred4")) +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") + 
  ylab("Total reads (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_total_reads_relative_grid.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Total barcode read number mass~type heatmap
#p <- ggplot(cluster_reads_long_barcoded %>% filter(BarcodePresence == "Barcode"), aes(x=Replicate, y=Type)) +
#  facet_grid(Type~Structure, drop=T, scales="free", space="free") +
#  geom_tile(aes(fill=Reads)) +
#  scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="blue") +
#  theme_bw(base_size=6)+
#  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
#  theme(legend.position="right")   
#
#ggsave(filename = file.path("results",paste(args[1],"_total_barcodes_absolute_heatmap.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Barcoded read number according to mass~type, and structure, in GRID (barplot)
#p <- ggplot(cluster_reads_long_barcoded %>% filter(BarcodePresence == "Barcode"), aes(x=Structure, y=Reads, fill=Structure)) +
#  scale_fill_manual(values = c("tomato","olivedrab","lightsalmon")) +
#  facet_grid(Type~Mass, drop=T, scales="free_x", space="free") +
#  geom_bar(position="stack", stat="identity") +
#  xlab("Sample and Mass (mg)") + 
#  ylab("Total barcode reads") +
#  theme_bw(base_size=6)+
#  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#  theme(legend.position="none")   
#
#ggsave(filename = file.path("results",paste(args[1],"_total_barcodes_mass_per_type_and_structure_bar.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Barcoded read number according to mass~type, and structure, in GRID (lineplot)
#p <- ggplot(cluster_reads_long_barcoded %>% filter(BarcodePresence == "Barcode"), aes(x=Mass, y=Reads, fill=Structure, color=Structure)) +
#  scale_fill_manual(values = c("tomato","olivedrab","lightsalmon")) +
#  facet_grid(Type~., drop=T, scales="free_x", space="free") +
#  geom_line() +
#  geom_point(size=1) +
#  xlab("Mass (mg)") + 
#  ylab("Total barcode reads") +
#  theme_bw(base_size=6)+
#  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#  theme(legend.position="right")   
#
#ggsave(filename = file.path("results",paste(args[1],"_total_barcodes_mass_per_type_and_structure_line.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Barcoded read number according to surface~type, and structure, in GRID (lineplot)
#p <- ggplot(cluster_reads_long_barcoded %>% filter(BarcodePresence == "Barcode"), aes(x=Surface, y=Reads, fill=Structure, color=Structure)) +
#  scale_fill_manual(values = c("tomato","olivedrab","lightsalmon")) +
#  facet_grid(Type~., drop=T, scales="free_x", space="free") +
#  geom_line() +
#  geom_point(size=1) +
#  xlab("Surface area (cm^2)") + 
#  ylab("Total barcode reads") +
#  theme_bw(base_size=6)+
#  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#  theme(legend.position="right")   
#
#ggsave(filename = file.path("results",paste(args[1],"_total_barcodes_surface_per_type_and_structure_line.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


###################################

#if a samples has no barcodes at all, it will have NaN in the dataframe, which makes it not being plotted. to remove NaN:
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

clusters_pct_threshold_meta_long[is.nan(clusters_pct_threshold_meta_long)] <- 0

#assign and order factor levels
#clusters_pct_threshold_meta_long$Library <- factor(clusters_pct_threshold_meta_long$Library, levels = sort(unique(clusters_pct_threshold_meta_long$Library)))
#clusters_pct_threshold_meta_long$Sample <- factor(clusters_pct_threshold_meta_long$Sample, levels = unique(clusters_pct_threshold_meta_long$Sample))
#clusters_pct_threshold_meta_long$Structure <- factor(clusters_pct_threshold_meta_long$Structure, levels = c("Inoculum","Primary","Secondary","Tertiary"))
#clusters_pct_threshold_meta_long$Type <- factor(clusters_pct_threshold_meta_long$Type, levels = c("Tip","Hair","Segment","Inoculum"))

plotting_df <- clusters_pct_threshold_meta_long %>% group_by(Library,Sample,Medium,Substrate) %>% summarize(Percentage=sum(Percentage))
write_csv(plotting_df,"plotting_df.csv")

#PLOT - Relativ abundance of all barcodes above threshold (in relation to total barcode read number) 
p <- ggplot(clusters_pct_threshold_meta_long %>% group_by(Library,Sample,Medium,Substrate) %>% 
              summarize(Percentage=sum(Percentage)), aes(x=Sample, y=Percentage, fill=Library)) +
  scale_fill_manual(values = library_colors) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") + 
  ylab("Above threshold barcodes (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_threshold_barcodes_relative.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)

#PLOT - Relativ abundance of all barcodes above threshold (in relation to total barcode read number) in GRID
p <- ggplot(clusters_pct_threshold_meta_long %>% group_by(Library,Replicate,Medium,Substrate) %>% 
              summarize(Percentage=sum(Percentage)), aes(x=Replicate, y=Percentage, fill=Library)) +
  scale_fill_manual(values = library_colors) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  xlab("Sample") + 
  ylab("Above threshold barcodes (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_threshold_barcodes_relative_grid.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)



#PLOT - Pie charts of percentage of all barcodes identified in the meta_barcodes file
#for some reason this doesnt show single replicates
#colourCount = length(list_known_libraries)+1                    #calculate number of colors needed for this plot
#getPalette = colorRampPalette(brewer.pal(8, "Set2"))         #create color gradient from x colors of set (https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html#rcolorbrewer-package)
#
#p <- ggplot(clusters_pct_threshold_meta_long %>% filter(.$Genotype %in% list_known_libraries)  %>% group_by(Genotype,Replicate,Type,Structure) %>% 
#              summarize(Percentage=sum(Percentage)), aes(x="Replicate", y=Percentage, fill=Genotype)) +
#  scale_fill_manual(values = library_colors) +
#  scale_y_continuous(breaks=seq(0,100,20)) +
#  geom_bar(position="fill", stat="identity") +
#  coord_polar(theta = "y") +
#  facet_grid(Type~Structure) +
#  theme_bw(base_size=6)+
#  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#  theme(legend.position="right")   
#
#ggsave(filename = file.path("results",paste(args[1],"_defined_barcodes_relative_pie.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


#PLOT - Bar charts of percentage of all barcodes identified in the meta_barcodes file
p <- ggplot(clusters_pct_threshold_meta_long %>% filter(.$Library %in% list_known_libraries)  %>% group_by(Library,Sample) %>% 
              summarize(Percentage=sum(Percentage)), aes(x=Sample, y=Percentage, fill=Library)) +
  scale_fill_manual(values = library_colors) +
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") + 
  ylab("Defined barcodes (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_defined_barcodes_relative_bar.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)

#PLOT - Bar charts of percentage of all barcodes identified in the meta_barcodes file coordinated in grid Substrate vs Medium
p <- ggplot(clusters_pct_threshold_meta_long %>% filter(.$Library %in% list_known_libraries)  %>% group_by(Library,Replicate,Substrate,Medium) %>% 
              summarize(Percentage=sum(Percentage)), aes(x=Replicate, y=Percentage, fill=Library)) +
  scale_fill_manual(values = library_colors) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  xlab("Replicate") + 
  ylab("Defined barcodes (relative)") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="right")   
if (debug==1) {p}
ggsave(filename = file.path("results",paste(args[1],"_defined_barcodes_relative_bar_grid.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)


###################################
######## DISTRIBUTION #############
###################################

#compare inoculum barcode distribution vs experimental samples

#make a results directory
write("\n\n----------", stdout())
write("Creating directory to compare barcode distributions", stdout())
dir.create("results/distribution", showWarnings = FALSE)

#create boxplots of the distribution of percentages of detected barcodes (!!!remove percentage==0!!!)
p <- ggplot(clusters_pct_threshold_meta_long %>% filter(Percentage != 0), aes(x = Replicate, y = Percentage)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(fill = "#4271AE", outlier.colour = "red", alpha = 0.9) +
  coord_cartesian(ylim=c(0,50)) +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"_barcode_abundance_distribution_boxplots.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)

p <- ggplot(clusters_pct_threshold_meta_long %>% filter(Percentage != 0), aes(x = Replicate, y = Percentage)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(fill = "#4271AE", outlier.colour = "red", alpha = 0.9) +
  coord_cartesian(ylim=c(0,0.5)) +
  facet_grid(Substrate~Medium, drop=T, scales="free_x", space="free") +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"_barcode_abundance_distribution_boxplots_zoomed.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)




#create df to work with (based on reads, not percentages)
histo_df <- cluster.csv %>% select(Center,all_of(sample_column_first):all_of(sample_column_last)) %>% 
                            left_join(.,cluster_dict) %>% 
                            gather(key="Sample", value="Reads",all_of(sample_column_first):all_of(sample_column_last)) %>% 
                            filter(Reads!=0) %>%
                            separate(.,Sample,c(NA,NA,"Sample")) %>% 
                            transform(., Sample = as.integer(Sample)) %>%
                            left_join(.,meta_samples) %>% select(Center,Library,Medium,Substrate,Replicate,Reads)

#where should y-axis be cut off?
histo_cutoff <- 1000

#inoculum
p <- ggplot(histo_df %>% filter(Medium=="Inoc"), aes(x=Reads)) + 
  geom_histogram(color="red", fill="white", binwidth=1) +
  xlab("Reads per Barcode") +
  ylab("Number of Barcodes") +
  coord_cartesian(ylim=c(0, histo_cutoff)) +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"histo_inoculum_full.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)




for (Medium_ in levels(clusters_pct_threshold_meta_long$Medium)[-1]) { 
  
  for (Substrate_ in levels(clusters_pct_threshold_meta_long$Substrate)[-4]) {
    
      for (replicate_ in levels(as.factor(clusters_pct_threshold_meta_long$Replicate))) {
    
      histograms <- ggplot(histo_df %>% filter(Medium==Medium_, Substrate==Substrate_, Replicate==replicate_), aes(x=Reads)) + 
                geom_histogram(color="red", fill="white", binwidth=1) +
                xlab("Reads per Barcode") + 
                ylab("Number of Barcodes") +
                coord_cartesian(ylim=c(0, histo_cutoff)) +
                theme_bw(base_size=6)+
                theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
                theme(legend.position="bottom")  
      ggsave(filename = file.path("results/distribution",paste(paste(args[1],"histo",Medium_,Substrate_,replicate_,sep="_"),".pdf",sep="")), plot = histograms, device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)
  
      }
  }
}
rm(temp,list_Mediums,list_Substrates)


p <- ggplot(histo_df %>% filter(Medium=="MOPS", Substrate=="Glucose", Replicate=="A"), aes(x=Reads)) + 
  geom_histogram(color="red", fill="white", binwidth=1) +
  xlab("Reads per Barcode") + 
  ylab("Number of Barcodes") +
  coord_cartesian(ylim=c(0, histo_cutoff)) +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
#if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"histo_primary_tip_a.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)

p <- ggplot(histo_df %>% filter(Medium=="Primary", Substrate=="Tip", Replicate=="B"), aes(x=Reads)) + 
  geom_histogram(color="red", fill="white", binwidth=1) +
  xlab("Reads per Barcode") + 
  ylab("Number of Barcodes") +
  coord_cartesian(ylim=c(0, histo_cutoff)) +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
#if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"histo_primary_tip_b.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)

p <- ggplot(histo_df %>% filter(Medium=="Primary", Substrate=="Tip", Replicate=="C"), aes(x=Reads)) + 
  geom_histogram(color="red", fill="white", binwidth=1) +
  xlab("Reads per Barcode") + 
  ylab("Number of Barcodes") +
  coord_cartesian(ylim=c(0, histo_cutoff)) +
  theme_bw(base_size=6)+
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(legend.position="bottom")   
#if (debug==1) {p}
ggsave(filename = file.path("results/distribution",paste(args[1],"histo_primary_tip_c.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)



write("\n\n----------", stdout())
write("Analysis completed", stdout())
Sys.Date()
Sys.time()

write("\n\n----------", stdout())
sessionInfo()

sink()

write("\n\n----------", stdout())
write("Analysis completed", stdout())
#*********************************************************************************************
#*********************************************************************************************
#*********************************************************************************************
q()
