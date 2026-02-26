# eDNA_pipeline

Bash pipeline for eDNA metabarcoding studies.
Structure of the git repository:
 - main.sh : main script doing most of the job
 - scripts/ : folder containing helping R scripts
 - env/ : folder containing the yaml files for reproductible mamba/conda environments
 - 

# Structure of the pipeline

The pipeline is composed of the following steps:
 - QC and adapter removal
 - ASV resolution
 - taxonomic assignment (currently looking at [midori II](https://www.reference-midori.info/) )
 - [confidence scoring of taxonomic assignments](https://onlinelibrary.wiley.com/doi/full/10.1002/edn3.70077)
 - report

# Output

The output is composed of:
 - intermediate files
 - OTU.fa (contains ASVs)
 - rdp_midori_"$marker".txt taxonomic assignment
 - a report
 - logs for each step
