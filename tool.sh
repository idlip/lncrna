#!/usr/bin/env bash

# Script to identify lncRNA
LNCRNA_VERSION=0.1

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
    ${c_green}%s${c_reset} -t fastp '-w 8 -z -V' -a bwa '-t 6 -k 33' -i rawData  (with additional parameters for tool)


    Options:
      ${c_yellow}-i, --input${c_reset}
        Input Fastq files directory
      ${c_yellow}-t, --trim 'additional commands' ${c_reset}
        Pre-processing tool to trim bad quality reads                  (options: fastp or trimgalore)
      ${c_yellow}-r, --refgenome  ${c_reset}
        Reference genome file  (recommended to keep in directory)
      ${c_yellow}-a, --aligner 'additional commands'  ${c_reset}
        Aligner tool to be used                                        (options: bwa, minimap2, hisat2, bowtie2, star)
      ${c_yellow}-n, --annotate  ${c_reset}
        lncRNA annotation GTF file



      ${c_yellow}-v, --version${c_reset}
        Print version of the tool

        " "${0##*/}" "${0##*/}" "${0##*/}" # to printf the bash file name
	exit 0
}

function check_fastq() {
	echo -e "${c_cyan}Check:${c_reset} Validating the Input directory"
	# INPUT_FASTQ_FILES="$(find "$INPUT_FASTQ_DIR" -type f \( -name "*.gz" -o -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" \) | sort -V)"
	# INPUT_FASTQ_FILES=($(ls -v "${INPUT_FASTQ_DIR}"/*.gz | uniq))
	mapfile -t INPUT_FASTQ_FILES < <(find "$INPUT_FASTQ_DIR" -type f \( -name "*.gz" -o -name "*.fastq" -o -name "*.fq" \) | sort -Vk1)
	# mapfile -t INPUT_FASTQ_FILES < <((ls -v "${INPUT_FASTQ_DIR}"/*.gz | uniq) 2> /dev/null)
	if [ "$(printf "%s" "${INPUT_FASTQ_FILES[@]}" | wc -c)" -gt 2 ]; then # if [$(find "$INPUT_FASTQ_DIR" -type f \( -name "*.gz" -o -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" \) | wc -l) ]
		printf "files are legit\n"
	else
		echo -e "${c_red}${c_bold}Error: ${c_reset}The input directory does not have any Fastq files"
		exit 0
	fi
	# mkdir -p qc_reports trimmed_data aligned_output/{sam,bam}
}

function version_info() {
	printf "${0##*/} ${c_magenta}%s\n" $LNCRNA_VERSION
	exit 0
}

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

# array of dependency tools required
_dependency=(
	"fastp"
	"bwa"
	"bowtie2"
	"samtools"
	"lolcat"
)

# check dependency package and install it
# _check_pkgs() {
#     for pkg in "${_dependency[@]}"; do
#         if [ ! -x "$(command -v "$pkg")" ]; then
#             echo "Package $pkg not installed, Installing it now."
#             echo "Please provide the sudo password to run '$install_pkg $pkg' for installing the package:"
#             $install_pkg "$pkg"
#         fi
#     done
#     echo "All packages are ready!"
# }
# _check_pkgs

function parse_fastq() {
	# to parse and fetch ref + annotation file
	# need to use api to get info?
	echo hi
}

function print_fastq_files() {
	if [ "$FASTQ_READ" == "paired" ]; then
		echo "Reading files in pairs: "
		for ((i = 0; i < ${#INPUT_FASTQ_FILES[@]}; i += 2)); do
			echo "${INPUT_FASTQ_FILES[i]} ${INPUT_FASTQ_FILES[i + 1]}"
			basename "${INPUT_FASTQ_FILES[i]}" | cut -d'.' -f1
		done
	else
		echo "Reading single end files"
		for file in "${INPUT_FASTQ_FILES[@]}"; do
			echo "$file"
			basename "$file" | cut -d'.' -f1
		done
	fi

}

function pre_process() {
	echo -e "${c_cyan}PreProcess:${c_reset} Trimming the fastq reads using $TRIMMER_TOOL"
	if [ "$FASTQ_READ" == "paired" ]; then
		for ((i = 0; i < ${#INPUT_FASTQ_FILES[@]}; i += 2)); do
			TRIMOUT1="trimmed_output/$(basename "${INPUT_FASTQ_FILES[i]}")"
			TRIMOUT2="trimmed_output/$(basename "${INPUT_FASTQ_FILES[i + 1]}")"
			REPORT="qc_reports/$(basename "${file}" | cut -d'.' -f1).html"
			case "$TRIMMER_TOOL" in
			"fastp")
				fastp -i "${INPUT_FASTQ_FILES[i]}" -I "${INPUT_FASTQ_FILES[i + 1]}" -o "$TRIMOUT1" -O "$TRIMOUT2" -h "$REPORT" "$TRIMMER_OPTS"
				;;
			"trimgalore")
				trimgalore -q 20 --stringency 3 --gzip --length 20 --paired "${INPUT_FASTQ_FILES[i]}" "${INPUT_FASTQ_FILES[i + 1]}" -o trimmed_output/ "$TRIMMER_OPTS"
				;;
			esac
		done
	else
		for file in "${INPUT_FASTQ_FILES[@]}"; do
			TRIMOUT="trimmed_output/$(basename "${file}")"
			REPORT="qc_reports/$(basename "${file}" | cut -d'.' -f1).html"
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

function index_refgen() {
  echo -e "${c_cyan}Indexing:${c_reset} Indexing the reference genome ($REF_GENOME) usign $ALIGNER_TOOL"
	REF_DIR=$(dirname "$REF_GENOME")
	REF_NAME="echo $REF_GENOME | cut -d'.' -f1"
	case "$ALIGNER_TOOL" in
	"minimap2")
		minimap2 -d "$REF_NAME".mmi "$REF_GENOME" "$ALIGNER_OPTS"
		;;
	"hisat2")
		hisat2-build "$REF_GENOME" "$REF_DIR" "$ALIGNER_OPTS"
		;;
	"bowtie2")
		bowtie2-build "$REF_GENOME" "$REF_DIR" "$ALIGNER_OPTS"
		;;
	"bwa")
		bwa index "$REF_GENOME" "$ALIGNER_OPTS"
		;;
	# "star" | "STAR")

	"*")
		echo "Sorry, we do not support this aligner as of now. You can open an issue for adding it."
		exit 0
		;;
	esac
}

function aligner() {
  echo -e "${c_cyan}Aligning:${c_reset} Mapping the reads onto reference genome ($REF_GENOME) using $ALIGNER_TOOL"
	if [ "$FASTQ_READ" == "paired" ]; then
		for ((i = 0; i < ${#INPUT_FASTQ_FILES[@]}; i += 2)); do
			SAMOUT="aligned_output/sam/$(basename "${INPUT_FASTQ_FILES[i]}" | cut -d'.' -f1).sam"
			case "$ALIGNER_TOOL" in
			"minimap2")
				minimap2 -a "$REF_NAME".mmi "${INPUT_FASTQ_FILES[i]} ${INPUT_FASTQ_FILES[i + 1]}" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"hisat2")
				hisat2 -x "$REF_GENOME" -1 "${INPUT_FASTQ_FILES[i]}" -2 "${INPUT_FASTQ_FILES[i + 1]}" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bowtie2")
				bowtie2 -x "$REF_GENOME" -1 "${INPUT_FASTQ_FILES[i]}" -2 "${INPUT_FASTQ_FILES[i + 1]}" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bwa" | *)
				bwa mem "$REF_GENOME" "${INPUT_FASTQ_FILES[i]} ${INPUT_FASTQ_FILES[i + 1]}" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			esac

		done
	else
		for file in "${INPUT_FASTQ_FILES[@]}"; do
			SAMOUT="aligned_output/sam/$(basename "$file" | cut -d'.' -f1).sam"
			case "$ALIGNER_TOOL" in
			"minimap2")
				minimap2 -a "$REF_NAME".mmi "$file" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"hisat2")
				hisat2 -x "$REF_GENOME" "$file" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bowtie2")
				bowtie2 -x "$REF_GENOME" -U "$file" -S "$SAMOUT" "$ALIGNER_OPTS"
				;;
			"bwa" | *)
				bwa mem "$REF_GENOME" "$file" -o "$SAMOUT" "$ALIGNER_OPTS"
				;;
			esac

		done
	fi
	# echo "Indexing somehting..."
	# echo "Aligning all files from $input_fastq_dir onto ref.fna"
}

function sam2bam() {
	echo -e "${c_cyan}BAM:${c_reset} Converting, sorting and Indexing SAM file into BAM file format"
	for file in aligned_out/sam/*; do
    BAM_OUT="aligned_output/bam/$(basename "$file" .sam).bam"
		samtools view -Sb -o "$BAM_OUT" "$file"
		samtools sort -o "$BAM_OUT" "$file"
		samtools index "$BAM_OUT"
	done
}

function ann_genes() {
    echo -e "${c_cyan}Quantification:${c_reset} Quantify the lncRNA from annotation file ($ANN_FILE)"
	if [ "$FASTQ_READ" == "paired" ]; then
		featureCounts -p --countReadPairs -t exon -g gene_id -a "$ANN_FILE" -o counts.txt aligned_output/bam/*
	else
		featureCounts -t exon -g gene_id -a "$ANN_FILE" -o counts.txt aligned_output/bam/*
	fi
}

function diffge() {
	echo -e "${c_cyan}DGE:${c_reset} Deriving top differential gene expression using R DESeq2"
	Rscript --vanilla deseq.R
}

function bam2bed() {
    for file in aligned_out/bam/*; do
        BED_OUT="aligned_output/bed/$(basename "$file" .bam).bed"
        bedtools bamtobed -bedpe -i "$file" | bedtools sort -i stdin > "$BED_OUT"
    done
}

function mergebed() {
    BEDFILE="results/final.bed" #TODO FIXME
    filearg=""
    for file in $(ls -v aligned_output/bed/*); do
        # to chain the input argument together for all files. bedtools does not take * as expansion
        filearg="$filearg -i $file"
    done
    bedtools merge "$filearg" > "$BEDFILE"
}

function codingpotential() {
    cpat -r "$REF_GENOME" -g "$BEDFILE" -d "$LOGIT" -x "$HEXAMER" -o "results/cpat-out" "$CPAT_OPTS"
}

while [ $# -gt 0 ]; do
	case "$1" in
	-i | --input)
		INPUT_FASTQ_DIR="$2"
		check_fastq
		;;
	-h | --help)
		help_info
		;;
	-t | --trim)
		TRIMMER_TOOL="$2"
		TRIMMER_OPTS="$3"
		;;
	-r | --refgenome)
		REF_GENOME="$2"
		;;
	-a | --aligner)
		ALIGNER_TOOL="$2"
		ALIGNER_OPTS="$3"
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
  -c | --cpat)
      CPAT_OPTS="$2"
    ;;
		# *)
		#     help_info
		#     ;;
	esac
	shift
done

# aligner
print_fastq_files
