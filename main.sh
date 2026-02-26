#!/bin/bash

set -e

Help() {
echo "
This script is the main for the eDNA pipeline. Its arguments are :
	-i --input <PATH> to the input folder.
	-m --map <PATH> to the mapping file, LotuS3 style.
	-o --output <PATH> to the output folder.
	-t --threads <INT> number of threads to use.
	-1 -2 <STR> strings for forward and reverse determinant.
	-p --tech <454/miSeq/hiSeq/PacBio> read technology.
	-n --nanopore <FLAG> \u2691 #comment on fait ca
	-g --marker <CO1/CO2/CO3/Cytb/A6/A8/ND1/ND2/ND3/ND4/ND4L/ND5/ND6/srRNA(12S)/lrRNA(18S)> gene marker used.
	-h --help displays this help message and exits.
"
}

if [ $# -eq 0 ]
then
	Help
	exit 0
fi

declare -A markers=(
    ["CO1"]=1 ["CO2"]=1 ["CO3"]=1 ["Cytb"]=1 ["A6"]=1
    ["A8"]=1 ["ND1"]=1 ["ND2"]=1 ["ND3"]=1 ["ND4"]=1
    ["ND4L"]=1 ["ND5"]=1 ["ND6"]=1 ["srRNA"]=1 ["lrRNA"]=1
)

while [ $# -gt 0 ]
do
	case $1 in
	-i | --input) input="$2"
	shift 2;;
	-m | --map) map="$2"
	shift 2;;
	-o | --output) output="$2"
        shift 2;;
	-t | --threads) threads="$2"
	shift 2;;
	-1) R1="$2"
	shift 2;;
	-2) R2="$2"
	shift 2;;
	-g | --marker) marker="$2"
	shift 2;;
	-p | --tech) tech="$2"
	shift 2;;
	-h | --help) Help; exit 0;;
	-* | --*) unknown="$1"; echo -e "ERROR: unknown argument: $unknown. Exiting."; exit 1;;
	*) shift ;;
	esac
done

map=$(readlink -f "$map")
input=$(readlink -f "$input")
output=$(readlink -f "$output")

if [ -d "$output" ];
then
        echo -e "WARNING: output directory ${output} already exist."
else
        mkdir "$output"
fi


#####################################
#
#       INPUT CHECKS
#
#####################################

if [ ! -d "$input" ];
then
	echo -e "ERROR: input folder ${input} does not exist. Exiting"
	exit 1
elif [ ! -f "$map" ];
then
	echo -e "ERROR: input file ${map} does not exist. Exiting"
	exit 1
elif [ ! "${markers[$marker]}" ];
then
	echo -e "ERROR: marker gene not in the list of markers avalable.Exiting."
	exit 1
fi


#####################################
#
#	QC
#
#####################################

# RUN 1

mamba run -n fastqc fastqc -o "$output"/fastqc -t "$threads" $(readlink -f "$input"/*.fastq*)

mamba run -n fastqc multiqc -o "$output"/multiqc "$output"/fastqc/*


#####################################
#
#       LotuS3/DADA2
#
#####################################

pipeline_dir=$(readlink -f "$0")
pipeline_dir=$(dirname "$pipeline_dir")

#TODO: pertinence de faire un arbre taxo alors qu on a un gros melange despeces qui peuvent etre bien distantes ?
"$pipeline_dir"/tools/LotuS3/lotus3 -i "$input" -o "$output" -t "$threads" -m "$map" -p "$tech" -buildPhylo 0 -CL dada2 -tax_group eukarya -amplicon_type "$amplicon_type"

#####################################
#
#       TAXONOMY
#
#####################################

java -Xmx15g -jar "$pipeline_dir"/tools/rdp_classifier_2.14/dist/classifier.jar classify -t "$pipeline_dir"/tools/rdp_midori_"$marker"/rRNAClassifier.properties -o "$output"/rdp_midori_"$marker"_hier.tsv --format allrank "$output"/OTU.fna &> "$output"/rdp.log
#Rscript to phyloseq
mamb run -n eDNA_R_packages Rscript scripts/phyloseq_to_sample2species_report.R -p "$output"/phyloseq.Rdata -o "$output" -r "$output"/rdp_midori_"$marker"_hier.tsv

#####################################
#
#       CONFIDENCE SCORING
#
#####################################

#TODO: dl le script des chercheurs et essayer de faire les metadatas
mamba run -n eDNA_R_packages Rscript "$pipeline_dir"/scripts/Pipeline_taxo_quality_V5.R
# les fichiers complexes a produire vont l'etre par les auteurs du script donc ca sera facile a faire.
# Soyons patients !

#####################################
#
#       REPORT
#
#####################################


#TODO: utiliser un script phyloseq pour faire les alpha pcoa ou autres

#struct: echantillon metadata liste_especes_trouvees note_sur_5
