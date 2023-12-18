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
#    (!!! all values lower case !!!)
#    (!!! column 'Sample' may only contain integers !!!)
# 4) Metadata file for Barcodes
# 
# 
# 
# 
# 
# 
# 


write("Loading libraries", stdout())

suppressMessages(if (!require('tidyverse')) install.packages('tidyverse'))
suppressMessages(if (!require('dplyr')) install.packages('dplyr'))
suppressMessages(if (!require('RColorBrewer')) install.packages('RColorBrewer'))
suppressMessages(if (!require('vegan')) install.packages('vegan'))
suppressMessages(if (!require('janitor')) install.packages('janitor'))

suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(vegan))
suppressMessages(library(janitor))

#*********************************************************************************************
#********PREPARATION**************************************************************************
#*********************************************************************************************

#read in script arguments
args = commandArgs(trailingOnly=TRUE)

#for debugging
write("DEBUG MODE", stdout())
args[1] <- "lib_assoc"
args[2] <- "0"
args[3] <- "sample_metadata.csv"
#args[4] <- "barcode_metadata.csv"

if (length(args)<3) {
  stop("Specify 'Bartender' output file prefix (e.g. 'all_samples'), threshold barcode relative abundance in percent (e.g. '2'), metadata file for samples, and metadata file for barcodes. 4 Parameters total.", call.=FALSE)
}




#threshold percent value from argument
threshold <- as.numeric(args[2])

#create filenames of input files and for output files from frefix supplied as script argument

name_clusterfile <- paste(args[1],"_cluster.csv",sep="")
name_qualityfile <- paste(args[1],"_quality.csv",sep="")
name_singlefile <- paste(args[1],"_single.csv",sep="")
name_combinedfile <- paste(args[1],"_combined.csv",sep="")

#make list of all files
table_inputs <- c(name_clusterfile,name_qualityfile)
other_inputs <- c(name_singlefile,name_combinedfile)


#read in table files (cluster and quality)
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
meta_samples <- read.table(args[3], header=T, sep=",")
#meta_barcodes <- read.table(args[4], header=T, sep=",")

#remove naming tools
rm(name_clusterfile,name_qualityfile,name_singlefile,name_combinedfile)

#*********************************************************************************************
#********META ANALYSIS************************************************************************
#*********************************************************************************************


#calculate number of samples in input
write("Counting samples", stdout())
z<-as.numeric(length(cluster.csv)-3)
z


#make a results directory
write("Creating results directory", stdout())
dir.create("results", showWarnings = FALSE)


#process cluster.csv to sort barcodes that only appear in one single library into its respective dataframe
lib_gm17_wt_10k_lib <-    cluster.csv %>% filter(time_point_1  > 0 & time_point_2 == 0) %>% select(Center)
lib_gm17_wt_10k_lib$Sample <- 1
lib_gm17_tnsa_10k_lib <-      cluster.csv %>% filter(time_point_1 == 0 & time_point_2  > 0) %>% select(Center)
lib_gm17_tnsa_10k_lib$Sample <- 2



#create and write out a metadata file that contains only barcodes that are unique for one library
meta_barcodes <- rbind(lib_gm17_wt_10k_lib,lib_gm17_tnsa_10k_lib) %>%
                  left_join(.,meta_samples) %>% select(-ReadsSum,-Sample)
rm(lib_gm17_tnsa_10k_lib,lib_gm17_wt_10k_lib)
write_csv(meta_barcodes, "results/barcode_metadata.csv", col_names=T)


#create a dictionary of ALL barcodes, and add metadata 
suppressWarnings(suppressMessages(
barcode_dict <- cluster.csv %>% select(Center) %>% left_join(.,meta_barcodes) %>% arrange(Center)
))
#if no genotype can be assigned because the barcode is present in multiple library samples, overwrite with "promiscuous"
barcode_dict$Library <- ifelse(is.na(barcode_dict$Library), "promiscuous", as.character(barcode_dict$Library))
write_tsv(barcode_dict, file.path("results",paste(args[1],"_barcode_dict.tab",sep="")), col_names=T)


#create a table of all sequenced barcodes, regardless if they are unique or not, to create a list for each library
suppressWarnings(suppressMessages(
assoc_barcodes_not_exclusified <- cluster.csv %>% select(c(2,4:(z+3))) %>% 
  gather(key="Sample", value="Reads",time_point_1:time_point_2) %>% 
  separate(.,Sample,c(NA,NA,"Sample")) %>% 
  transform(., Sample = as.numeric(Sample)) %>% filter(Reads>0) %>% 
  left_join(.,meta_samples) %>% arrange(Center) %>% select(Center,Library)
))


#write out a list of barcodes for every library (unique or not)
dir.create("results/venn", showWarnings = FALSE)
write(assoc_barcodes_not_exclusified %>% filter(Library=="GM17_TnSA_10k") %>% .[,"Center"],   file = "results/venn/GM17_TnSA_10k.list")
write(assoc_barcodes_not_exclusified %>% filter(Library=="GM17_WT_10k") %>% .[,"Center"],     file = "results/venn/GM17_WT_10k.list")

#create a summary file with how many barcodes were identified, found to be unique, and promiscuous for each library
lib_summary <- meta_samples %>% select(Library,LibSizeIntend) %>%
  left_join(assoc_barcodes_not_exclusified %>% count(Library,name="Total")) %>% 
  left_join(meta_barcodes %>% count(Library,name="Unique"))
lib_summary[is.na(lib_summary)] <- 0
lib_summary$Promiscuous <- lib_summary$Total - lib_summary$Unique
write_tsv(lib_summary, file.path("results",paste(args[1],"_lib_summary.tab",sep="")), col_names=T)

#plot rarefaction curve
#prepare input file for vegan::rarecurve()
#make long - seperate time_point_# to # - make integer - join sample data - make wide - convert Library column to row names
rarefac_df <- cluster.csv %>% select(Center,time_point_1:time_point_2) %>% 
                        pivot_longer(cols = time_point_1:time_point_2,names_to = "Sample") %>% 
                        separate(.,Sample,c(NA,NA,"Sample")) %>%
                        transform(., Sample = as.integer(Sample)) %>% 
                        left_join(meta_samples) %>%
                        select(Center,value,Library) %>%
                        pivot_wider(names_from = Center, values_from = value) %>%
                        column_to_rownames(var="Library")

#plot curves
pdf("results/rarefaction_curves_vegan.pdf",paper="a4r", useDingbats = FALSE)
rarecurve(rarefac_df,step = 1000,xlab="Sequencing reads", ylab="Identified Barcodes", label=F, cex = 0.6)
rarecurve(rarefac_df,step = 1000,xlab="Sequencing reads", ylab="Identified Barcodes", label=T, cex = 0.6)
dev.off()

# #*********************************************************************************************
# #********ANALYSIS*****************************************************************************
# #*********************************************************************************************
# 
# write("Plotting", stdout())
# #ON ABSOLUTE READ BASIS
# 
# #per sample, add up all reads that contain a barcode, make long, and then separate 'time_point_#' to '#', add metadata
# cluster_reads_long <- cluster.csv %>% select(c(1,4:(z+3))) %>% 
#                           summarise_at(vars(1:(z+1)), sum, na.rm = TRUE) %>% 
#                           gather(key="Sample", value="Barcode",2:(z+1)) %>% 
#                           separate(.,Sample,c(NA,NA,"Sample")) %>% 
#                           transform(., Sample = as.numeric(Sample)) %>% left_join(.,meta_samples)
# 
# #calculate how many reads do not contain barcodes
# cluster_reads_long$Unknown <- cluster_reads_long$ReadsSum - cluster_reads_long$Barcode
# 
# #ON RELATIVE ABUNDANCE BASIS
# 
# #calculate relative abundance of barcodes (in percent) from absolute read numbers and add to original cluster.csv
# for (i in 4:(length(cluster.csv))) {
# sapply(names(cluster.csv)[i], function(x) {
#   cluster.csv[paste0("Sample_", i-3, "_pct")] <<- cluster.csv[x] / sum(cluster.csv[x]) * 100
# })
# }
# 
# #create list of cluster.ids that are have percentage >= threshold in ANY sample
# list_cluster_pct_threshold <- as.character (cluster.csv %>% 
#                                             select(c(2,(z+4):(3+(2*z)))) %>% 
#                                             gather(key="Sample", value="Percentage",2:(z+1)) %>% 
#                                             filter(Percentage >= threshold) %>%
#                                             select(1) %>% arrange(Center) %>%
#                                             unique() %>% t(.)
#                                           )
# 
# 
# #delete all rows that have cluster.ids that are not listed in the list created above
# cluster_pct_threshold <- cluster.csv %>% select(c(2,(z+4):(3+(2*z)))) %>% filter(.[[1]] %in% list_cluster_pct_threshold)
# 
# #expand filtered data set by barcode metadata obtained from barcode_dict
# #make long and separate 'sample_3_pct' into '#'
# #add sample metadata
# 
# suppressWarnings(suppressMessages(
# clusters_pct_threshold_meta_long <- cluster_pct_threshold %>% 
#                                         left_join(.,barcode_dict) %>% 
#                                         gather(key="Sample", value="Percentage",2:(z+1)) %>% 
#                                         separate(.,Sample,c(NA,"Sample",NA)) %>% 
#                                         transform(., Sample = as.numeric(Sample)) %>%
#                                         left_join(.,meta_samples) 
# ))
# 
# 
# options(scipen=999)
# #*********************************************************************************************
# #********PLOTTING*****************************************************************************
# #*********************************************************************************************
# 
# cluster_reads_long_barcoded <- cluster_reads_long %>% select(Sample,Library,Species,Genotype,LibSizeIntend,Replicate,ReadsSum,Barcode,Unknown) %>% gather(key="BarcodePresence", value="Reads",Barcode:Unknown)
# 
# cluster_reads_long_barcoded$Sample <- factor(cluster_reads_long_barcoded$Sample, levels = unique(cluster_reads_long_barcoded$Sample))
# cluster_reads_long_barcoded$BarcodePresence <- factor(cluster_reads_long_barcoded$BarcodePresence, levels = c("Unknown","Barcode"))
# #cluster_reads_long_barcoded$Structure <- factor(cluster_reads_long_barcoded$Structure, levels = c("primary","secondary","tertiary"))
# #cluster_reads_long_barcoded$Type <- factor(cluster_reads_long_barcoded$Type, levels = c("tips","hair","segment"))
# 
# 
# #PLOT - Total read number, total barcode read number, for 10k libraries, in GRID
# p <- ggplot(cluster_reads_long_barcoded %>% filter(LibSizeIntend==10000), aes(x=Species, y=Reads, fill=BarcodePresence)) +
#   scale_fill_manual(values = c("Grey","violetred4")) +
#   facet_grid(Genotype~Species, drop=T, scales="free_x", space="free") +
#   geom_bar(position="stack", stat="identity") +
#   xlab("Sample") + 
#   ylab("Total reads (absolute)") +
#   theme_bw(base_size=6)+
#   theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#   theme(legend.position="bottom")   
# p
# ggsave(filename = file.path("results",paste(args[1],"10k_libs_absolute.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)
# 
# 
# #PLOT - Percentage barcode reads of all reads, for 10k libraries, in GRID
# p <- ggplot(cluster_reads_long_barcoded %>% filter(LibSizeIntend==10000), aes(x=Species, y=Reads, fill=BarcodePresence)) +
#   scale_fill_manual(values = c("Grey","violetred4")) +
#   facet_grid(Genotype~Species, drop=T, scales="free_x", space="free") +
#   geom_bar(position="fill", stat="identity") +
#   xlab("Sample") + 
#   ylab("Total reads (absolute)") +
#   theme_bw(base_size=6)+
#   theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#   theme(legend.position="bottom")   
# p
# ggsave(filename = file.path("results",paste(args[1],"10k_libs_perc.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)
# 
# 
# 
# 
# 
# #PLOT - distribution of barcodes in each library
# barcode_dist <- clusters_pct_threshold_meta_long[order(clusters_pct_threshold_meta_long$Genotype,clusters_pct_threshold_meta_long$Center),] %>% filter(Genotype!="no or multiple assignments")
# barcode_dist$Center <- as.factor(barcode_dist$Center)
# barcode_dist$newID <- seq.int(nrow(barcode_dist))
# barcode_dist$Center <- reorder(barcode_dist$Center, barcode_dist$newID)
# 
# p2 <- ggplot(barcode_dist, aes(x=Library, y=Percentage, fill=Center))+
#   geom_bar(position="fill", stat="identity")+
#   theme_bw(base_size=6)+
#   theme(legend.position = "none")
# 
# ggsave(filename = file.path("results",paste(args[1],"_barcode_dist.pdf",sep="")), plot = p2, device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)
# rm(p2)
# 
# histo_df <- barcode_dist %>% filter(Percentage!=0)
# histo_gm17tn <- histo_df %>% filter(Library=="GM17_TnSA_10k")
# histo_gm17wt <- histo_df %>% filter(Library=="GM17_WT_10k")
# histo_kt2440tn <- histo_df %>% filter(Library=="KT2440_TnSA_10k")
# histo_kt2440wt <- histo_df %>% filter(Library=="KT2440_WT_10k")
# pdf("results/histograms.pdf")
# hist(histo_df$Percentage, xlim=c(0,0.1), breaks=100, col="violetred4", xlab = "Abundance of Barcode", ylab = "Frequency")
# hist(histo_gm17tn$Percentage, xlim=c(0,0.1), breaks=100, col="violetred4", xlab = "Abundance of Barcode", ylab = "Frequency")
# hist(histo_gm17wt$Percentage, xlim=c(0,0.1), breaks=100, col="violetred4", xlab = "Abundance of Barcode", ylab = "Frequency")
# hist(histo_kt2440tn$Percentage, xlim=c(0,0.1), breaks=100, col="violetred4", xlab = "Abundance of Barcode", ylab = "Frequency")
# hist(histo_kt2440wt$Percentage, xlim=c(0,0.1), breaks=100, col="violetred4", xlab = "Abundance of Barcode", ylab = "Frequency")
# dev.off()
# 
# ###################################
# 
# #clusters_pct_threshold_meta_long$Name <- factor(clusters_pct_threshold_meta_long$Name, levels = unique(clusters_pct_threshold_meta_long$Name))
# clusters_pct_threshold_meta_long$Sample <- factor(clusters_pct_threshold_meta_long$Sample, levels = unique(clusters_pct_threshold_meta_long$Sample))
# 
# #PLOT - Relativ abundance of all barcodes above threshold (in relation to total barcode read number) 
# colourCount = 5                    #calculate number of colors needed for this plot
# getPalette = colorRampPalette(brewer.pal(11, "Spectral"))         #create color gradient from x colors of set (https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html#rcolorbrewer-package)
# 
# p <- ggplot(clusters_pct_threshold_meta_long, aes(x=Sample, y=Percentage, fill=Genotype)) +
#   scale_fill_manual(values=getPalette(colourCount)) +
#   scale_y_continuous(breaks=seq(0,100,20)) +
#   geom_bar(position="stack", stat="identity") +
#   xlab("Sample") + 
#   ylab("Above threshold barcodes (relative)") +
#   theme_bw(base_size=6)+
#   theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
#   theme(legend.position="bottom")   
# 
# p
# ggsave(filename = file.path("results",paste(args[1],"_threshold_barcodes_relative.pdf",sep="")), plot = last_plot(), device = "pdf" ,width = 297, height = 210, units = "mm",useDingbats = FALSE)
# 
# 
# 
# #*********************************************************************************************
# #*********************************************************************************************
# #*********************************************************************************************
q()