########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("phyloseq")
  library("ggplot2")
}))

########################################################
#
#       ARGS
#
########################################################

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Folder with the phyloseq.Rdata file.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output directory file).n", call.=FALSE)
}
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input).n", call.=FALSE)
}

output_path = opt$output
input_path = opt$input

if (!dir.exists(output_path)) {
  # If it does not exist, create the folder
  dir.create(output_path, recursive = TRUE)  # Use recursive = TRUE to create any necessary parent directories
  cat("Folder created:", output_path, "\n")
} else {
  cat("Folder already exists:", output_path, "\n")
}

##################################################
#
#		LOADING FILES
#
##################################################

load(paste(input_path, "/phyloseq.Rdata", sep = ""))

##################################################
#
#		Figures
#
##################################################

#potentiellement filtrer sur les 20 especes les plus abondances, faire une figure sur un niveau atxo plus haut pour faciliter la lecture
png(filename = paste(output_path,"/barplot.png", sep = ""), width = 1000, height = 1000)
plot_bar(physeq, fill = "Genus")
dev.off()

png(filename = paste(output_path,"/phylogic_tree.png", sep = ""))
plot_tree(physeq, label.tips="taxa_names", ladderize="left", plot.margin=0.3)
dev.off()

png(filename = paste(output_path,"/heatmap.png", sep = ""))
plot_heatmap(physeq, title = "Heatmap of samples based on the OTU/ASV composition.")
dev.off()

#estimate_richness(physeq, split = TRUE, measures = NULL)
#eventuellement donner un argument pour la liste des mesures de diversite et en donner le choix
png(filename = paste(output_path,"/richness.png", sep = "", width=1600, height=900))
plot_richness(physeq, x = "samples", color = NULL, shape = NULL,
title = "Plot of the richness per sample.", scales = "free_y", nrow = 1, shsi = NULL,
measures = NULL, sortby = NULL)
dev.off()
