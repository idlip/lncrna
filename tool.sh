#!/usr/bin/env bash

# Script to identify lncRNA
LNCRNA_VERSION=0.1

# set -e ## (to exit if any command does not run)
# set -x ## to debug and know command is run during exec

# colors echo -e {{{
c_red="\033[1;31m"
c_green="\033[1;32m"
c_yellow="\033[1;33m"
c_blue="\033[1;34m"
c_magenta="\033[1;35m"
c_cyan="\033[1;36m"
c_reset="\033[0m"
c_bold="\033[1m"
#}}}

# DBG_MSG="echo 'Notice here'"

# print help message for the tool
help_info() {
	printf "
    Usage:
    ${c_green}%s${c_reset} [options] [fastq dir]
    ${c_green}%s${c_reset} -t fastp -a bwa -i rawData

    Options:
      ${c_yellow}-c, --config ${c_reset}
        Configuration file path for the tool
      ${c_yellow}-i, --input ${c_reset}
        Input Fastq files directory
      ${c_yellow}-t, --trim ${c_reset}
        Pre-processing tool to trim bad quality reads                  (options: fastp or trimgalore)
      ${c_yellow}-r, --refgenome  ${c_reset}
        Reference genome file path
      ${c_yellow}-a, --aligner  ${c_reset}
        Aligner tool to be used                                        (options: bwa, minimap2, hisat2, bowtie2, star)
      ${c_yellow}-n, --annotate  ${c_reset}
        lncRNA annotation GTF file
      ${c_yellow}-@, --threads  ${c_reset}                             (Type: Integer)
        Threads for certain tools to be used for efficient computation
      ${c_yellow}-p, --paired  ${c_reset}
        Indicate that raw sequence files are paired end reads
      ${c_yellow}-v, --version${c_reset}
        Print version of the tool

        " "${0##*/}" "${0##*/}" # to printf the bash file name
	exit 0
}

###########################
##### Config template #####
###########################
function _init_config() {
	cat <<-EOF
		hello
		world

	EOF
	exit 0
}

#####################################################################
##### Variables/constant for defining reusable code dynamically #####
#####################################################################
: "${LNCRNA_CONF:=$PWD/lncrna.conf}"
function _parse_conf() {
	if [ -f "$LNCRNA_CONF" ]; then
		source "$LNCRNA_CONF"
	else
		echo -e "${c_red}Error:${c_reset} No config file found"
		exit 1
	fi
}

#################################
##### Check for Executables #####
#################################
_check_pkgs() {
	_dependency=(
		"$ALIGNER_TOOL" "$TRIMMER_TOOL"
	)
	for pkg in "${_dependency[@]}"; do
		if [ ! -x "$(command -v "$pkg")" ]; then
			echo -e "Package ${c_yellow}$pkg${c_reset} not installed."
			echo "Please install it or use the preferred Nix package method."
			exit 0
		fi
	done
}

###############################
##### Check for all files #####
###############################
_check_files() {
	_data_files=(
		"$REFERENCE_GENOME" "$LOGIT_FILE" "$HEXAMER_FILE"
	)
	for file in "${_data_files[@]}"; do
		[ -f "$file" ] ||
			(echo -e "${c_cyan}Error:${c_reset} $file file not found" &&
				exit 1)
	done

}

##########################################################################
##### Function to check for gz|fastq files in given input directory #####
##########################################################################
function check_fastq() {
	echo -e "${c_cyan}Check:${c_reset} Validating the Input directory"
	# RAWDATA_FILES="$(find "$RAWDATA_DIR" -type f \( -name "*.gz" -o -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" \) | sort -V)"
	# RAWDATA_FILES=($(ls -v "${RAWDATA_DIR}"/*.gz | uniq))
	mapfile -t RAWDATA_FILES < <(find "$RAWDATA_DIR" -type f \( -name "*.gz" -o -name "*.fastq" -o -name "*.fq" \) | sort -Vk1)
	### so many ways, depends on best way. But all will work enough!! ###
	# mapfile -t RAWDATA_FILES < <((ls -v "${RAWDATA_DIR}"/*.gz | uniq) 2> /dev/null)
	if [ "$(printf "%s" "${RAWDATA_FILES[@]}" | wc -c)" -gt 2 ]; then # if [$(find "$RAWDATA_DIR" -type f \( -name "*.gz" -o -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" \) | wc -l) ]
		printf "files are legit\n"
	else
		echo -e "${c_red}${c_bold}Error: ${c_reset}The input directory does not have any Fastq files"
		exit 0
	fi
	# mkdir -p $QC_DIR trimmed_data $ALIGNED_DIR/{sam,bam}
}

function version_info() {
	printf "${0##*/} ${c_magenta}%s\n" $LNCRNA_VERSION
	exit 0
}

##### prefer nix or snakemake or conda or ...? Soemthing universal and dependable
# # check distro used and determine the install command to be used
# _check_distro() {
#     _distro="$(grep -i "^ID=" /etc/os-release | cut -d'=' -f2)"
#     case $_distro in
#         "debian" | "ubuntu" | "linuxmint" | "pop" | "kali")
#             install_pkg="sudo apt install"
#             ;;
#         "arch")
#             install_pkg="sudo pacman -S"
#             ;;
#         "centos" | "fedora")
#             install_pkg="sudo dnf...."
#             ;;
#         "nixos")
#             install_pkg="echo really want to"
#             ;;
#     esac
# }
# _check_distro

function parse_fastq() {
	# to parse and fetch ref + annotation file
	# need to use api to get info?
	# As on 2024-08-07, we need manual input of paired or single end, rest taken care!
	echo ""
}

################################################################################
######### List the Fastq files in fashion of sinlge or paired manner ###########
################################################################################
function print_fastq_files() {
	if [ "$FASTQ_READ" == "paired" ]; then
		echo "Reading files in pairs: "
		for ((i = 0; i < ${#RAWDATA_FILES[@]}; i += 2)); do
			echo "${RAWDATA_FILES[i]} ${RAWDATA_FILES[i + 1]}"
			basename "${RAWDATA_FILES[i]}" | cut -d'.' -f1
		done
	else
		echo "Reading single end files"
		for file in "${RAWDATA_FILES[@]}"; do
			echo "$file"
			basename "$file" | cut -d'.' -f1
		done
	fi

}

######################################################################
######## Pre-processing the raw reads to trim bad quality ############
######################################################################
function pre_process() {
	echo -e "${c_cyan}PreProcess:${c_reset} Trimming the fastq reads using $TRIMMER_TOOL"
	if [ "$FASTQ_READ" == "paired" ]; then
		for ((i = 0; i < ${#RAWDATA_FILES[@]}; i += 2)); do
			TRIMOUT1="$TRIM_DIR/$(basename "${RAWDATA_FILES[i]}")"
			TRIMOUT2="$TRIM_DIR/$(basename "${RAWDATA_FILES[i + 1]}")"
			REPORT="$QC_DIR/$(basename "${file}" | cut -d'.' -f1).html"
			case "$TRIMMER_TOOL" in
			"fastp")
				fastp -i "${RAWDATA_FILES[i]}" -I "${RAWDATA_FILES[i + 1]}" -o "$TRIMOUT1" -O "$TRIMOUT2" -h "$REPORT" "$TRIMMER_OPTS"
				;;
			"trimgalore")
				trimgalore -q 20 --stringency 3 --gzip --length 20 --paired "${RAWDATA_FILES[i]}" "${RAWDATA_FILES[i + 1]}" -o $TRIM_DIR/ "$TRIMMER_OPTS"
				;;
			esac
		done
	else
		for file in "${RAWDATA_FILES[@]}"; do
			TRIMOUT="$TRIM_DIR/$(basename "${file}")"
			REPORT="$QC_DIR/$(basename "${file}" | cut -d'.' -f1).html"
			case "$TRIMMER_TOOL" in
			"fastp")
				fastp -i "$file" -o "$TRIMOUT" -h "$REPORT" "$TRIMMER_OPTS"
				;;
			"trimgalore")
				trimgalore -q 20 --stringency 3 --gzip --length 20 "$file" -o "$TRIMOUT" "$TRIMMER_OPTS"
				;;
			esac
		done
	fi
}

######################################################################################
##### Index the reference genome using a tool of choice for paired or single end #####
######################################################################################
function index_refgen() {
	echo -e "${c_cyan}Indexing:${c_reset} Indexing the reference genome ($REFERENCE_GENOME) usign $ALIGNER_TOOL"
	REF_DIR=$(dirname "$REFERENCE_GENOME")
	REF_NAME="echo $REFERENCE_GENOME | cut -d'.' -f1"
	case "$ALIGNER_TOOL" in
	"minimap2")
		minimap2 -d "$REF_NAME".mmi "$REFERENCE_GENOME" "$ALIGNER_OPTS"
		;;
	"hisat2")
		hisat2-build "$REFERENCE_GENOME" "$REF_DIR" "$ALIGNER_OPTS"
		;;
	"bowtie2")
		bowtie2-build "$REFERENCE_GENOME" "$REF_DIR" "$ALIGNER_OPTS"
		;;
	"bwa")
		bwa index "$REFERENCE_GENOME" "$ALIGNER_OPTS"
		;;
	# "star" | "STAR")

	"*")
		echo "Sorry, we do not support this aligner as of now. You can open an issue for adding it."
		exit 0
		;;
	esac
}

##############################################################################
##### Align the sequence against the reference genome using the same tool #####
##############################################################################
function aligner() {
	echo -e "${c_cyan}Aligning:${c_reset} Mapping the reads onto reference genome ($REFERENCE_GENOME) using $ALIGNER_TOOL"
	if [ "$FASTQ_READ" == "paired" ]; then
		for ((i = 0; i < ${#RAWDATA_FILES[@]}; i += 2)); do
			SAMOUT="$ALIGNED_DIR/sam/$(basename "${RAWDATA_FILES[i]}" | cut -d'.' -f1).sam"
			case "$ALIGNER_TOOL" in
			"minimap2")
				minimap2 -a "$REF_NAME".mmi "${RAWDATA_FILES[i]} ${RAWDATA_FILES[i + 1]}" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"hisat2")
				hisat2 -x "$REFERENCE_GENOME" -1 "${RAWDATA_FILES[i]}" -2 "${RAWDATA_FILES[i + 1]}" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bowtie2")
				bowtie2 -x "$REFERENCE_GENOME" -1 "${RAWDATA_FILES[i]}" -2 "${RAWDATA_FILES[i + 1]}" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bwa" | *)
				bwa mem "$REFERENCE_GENOME" "${RAWDATA_FILES[i]} ${RAWDATA_FILES[i + 1]}" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			esac

		done
	else
		for file in "${RAWDATA_FILES[@]}"; do
			SAMOUT="$ALIGNED_DIR/sam/$(basename "$file" | cut -d'.' -f1).sam"
			case "$ALIGNER_TOOL" in
			"minimap2")
				minimap2 -a "$REF_NAME".mmi "$file" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"hisat2")
				hisat2 -x "$REFERENCE_GENOME" "$file" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bowtie2")
				bowtie2 -x "$REFERENCE_GENOME" -U "$file" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bwa" | *)
				bwa mem "$REFERENCE_GENOME" "$file" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			esac

		done
	fi
	# echo "Indexing somehting..."
	# echo "Aligning all files from $input_fastq_dir onto ref.fna"
}

############################################################
##### Convert sam to bam and save space & be efficient #####
############################################################
function sam2bam() {
	echo -e "${c_cyan}BAM:${c_reset} Converting, sorting and Indexing SAM file into BAM file format"
	for file in "$ALIGNED_DIR"/sam/*; do
		BAM_OUT="$ALIGNED_DIR/bam/$(basename "$file" .sam).bam"
		samtools view -Sb -o "$BAM_OUT" "$file"
		samtools sort -o "$BAM_OUT" "$file"
		samtools index "$BAM_OUT"
	done
}

#########################################################################
##### Quantify the gene using featureCounts via annotation of exons #####
#########################################################################
function ann_genes() {
	echo -e "${c_cyan}Quantification:${c_reset} Quantify the lncRNA from annotation file ($ANN_FILE)"
	if [ "$FASTQ_READ" == "paired" ]; then
		featureCounts -p --countReadPairs -t exon -g gene_id -a "$ANN_FILE" -o counts.txt "$ALIGNED_DIR"/bam/*
	else
		featureCounts -t exon -g gene_id -a "$ANN_FILE" -o counts.txt "$ALIGNED_DIR"/bam/*
	fi
}

function diffge() {
	echo -e "${c_cyan}DGE:${c_reset} Deriving top differential gene expression using R DESeq2"
	Rscript --vanilla deseq.R
}

# better to merge bam files or bed files?
# caveat limit with bam would be of memory usage of streaming files

##############################################################
##### Loop all bam files and convert them into bed files #####
##############################################################
function bam2bed() {
	for file in "$ALIGNED_DIR"/bam/*; do
		BED_OUT="$ALIGNED_DIR/bed/$(basename "$file" .bam).bed"
		bedtools bamtobed -bed12 -i "$file" >"$BED_OUT"
	done
}

############################################################
##### Merge the many bed files using sort unix command #####
############################################################
function mergebed() {
	# method 1
	# filearg=""
	# for file in $(ls -v $ALIGNED_DIR/bed/*); do
	#     # to chain the input argument together for all files. bedtools does not take * as expansion
	#     filearg="$filearg -i $file"
	# done
	# cat $ALIGNED_DIR/bed/*.bed | bedtools sort -i stdin | bedtools merge -i stdin > "$BEDFILE"

	# method 2
	sort -u "$ALIGNED_DIR"/bed/*.bed -o "$BEDFILE"
}

##############################################################################################
##### Calculate the coding potential of merged bed file against built model of reference #####
##############################################################################################
function codingpotential() {
	cpat -r "$REFERENCE_GENOME" -g "$BEDFILE" -d "$LOGIT_FILE" -x "$HEXAMER_FILE" -o "$CPAT_OUTPUT_PREFIX" "$CPAT_OPTS"
	Rscript "ncp-filter.R" "results/final-out.ORF_prob.tsv" #FIXME
}

#####################################################################
##### Filter and extract the long non-coding ORF sequence files #####
#####################################################################
function extract_lncrna() {
	seqkit seq -m 200 "results/final-out.ORF_seqs.fa" -o "results/long_seqs.fa"
	seqkit grep -r -f "results/ncp-gene-list.tsv" "results/long_seqs.fa" -o "results/lncrna.fa"
}

# conditional loop over the arguments to move forward
while [ $# -gt 0 ]; do
	case "$1" in
	-i | --input)
		RAWDATA_DIR="$2"
		;;
	-h | --help)
		help_info
		;;
	--init-config)
		_init_config >"$LNCRNA_CONF"
		;;
	-t | --trim)
		TRIMMER_TOOL="$2"
		;;
	-r | --refgenome)
		REFERENCE_GENOME="$2"
		;;
	-a | --aligner)
		ALIGNER_TOOL="$2"
		;;
	-n | --annotate)
		ANN_FILE=$2
		;;
	-v | --version)
		version_info
		;;
	-p | --paired)
		FASTQ_READ="paired"
		;;
	-@ | --threads)
		THREADS="$2"
		;;
	-c | --config)
		LNCRNA_CONF="$2"
		;;
		# *)
		#     help_info
		#     ;;
	esac
	shift
done

_parse_conf
_check_pkgs
check_fastq

# aligner
print_fastq_files
# echo $ALIGNER_TOOL
