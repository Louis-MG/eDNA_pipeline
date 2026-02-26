########################################################
#
#       PACKAGES
#
########################################################

suppressWarnings(suppressPackageStartupMessages({
  library("optparse")
  library("tidyr")
  library("dplyr")
  library("phyloseq")
}))

########################################################
#
#       ARGS
#
########################################################

option_list = list(
  make_option(c("-p", "--phyloseq"), type="character", default=NULL,
              help="Phyloseq object from Lotus3.", metavar="character"),
  make_option(c("-r", "--rdp"), type="character", default=NULL,
              help="RDP hiera object with classification up to the species rank.", metavar="character"),  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output directory file).n", call.=FALSE)
}
if (is.null(opt$phyloseq)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (phyloseq).n", call.=FALSE)
}
if (is.null(opt$rdp)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rdp).n", call.=FALSE)
}

rdp_path=opt$rdp
folder_path = opt$output
phyloseq_path = opt$phyloseq

if (!dir.exists(folder_path)) {
  # If it does not exist, create the folder
  dir.create(output_path, recursive = TRUE)  # Use recursive = TRUE to create any necessary parent directories
  cat("Folder created:", folder_path, "\n")
} else {
  cat("Folder already exists:", folder_path, "\n")
}

##################################################
#
#		LOADING FILES
#
##################################################

data <- load(phyloseq_path)

classifier_output <- read.table(rdp_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(classifier_output) <- c("seq_id", "Root", "Kingdom", "rank2", "confidence2", "Phylum", "rank3", "confidence3", "Class", "rank4", "confidence4", "Order", "rank5", "confidence5", "Family", "rank6", "confidence6", "Genus", "rank7", "confidence7", "Species", "rank8", "confidence8")
tax_table <- classifier_output[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
tax_table[is.na(tax_table)] <- "?"

#remplacer la taxa_table du pyloseq par la nouvelle de rdp
tax_table <- as.matrix(tax_table)
rownames(tax_table) <- classifier_output$seq_id
ps_tax_table <- tax_table(tax_table)

tax_table(physeq) <- ps_tax_table

##################################################
#
#               SAVING PHYLOSEQ
#
##################################################

# writing phyloseq
save(physeq, file = paste(folder_path, "phyloseq.Rdata", sep = "/"))

##################################################
#
#               phyloseq2table
#
##################################################

melted_df <- psmelt(physeq)

# Create a table with sample names and the list of species present in each
species_per_sample <-
  melted_df %>%
  filter(Abundance > 0) %>%  # Only keep species with abundance > 0
  select(Sample, Species) %>%  # Adjust column names if needed
  group_by(Sample) %>%
  summarise(Species_Present = list(unique(Species))) %>%
  unnest(cols = c(Species_Present))

# Pivot to a wide format for a cleaner table
species_table <-
  species_per_sample %>%
  group_by(Sample) %>%
  summarise(Species_List = paste(Species_Present, collapse = ", "))

# Print the table

output_path = paste(folder_path, "sample2species.tsv", sep = "/")
write.table(species_table, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
