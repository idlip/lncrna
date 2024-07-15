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
    ${c_green}%s${c_reset} [options] [query]
    ${c_green}%s${c_reset} [query] [options]
    ${c_green}%s${c_reset} [options] [query] [options]

    Options:
      ${c_yellow}-i, --input
        Input Fastq files directory

        " "${0##*/}" "${0##*/}" "${0##*/}" # to printf the bash file name
	  exit 0
}

# check distro used and determine the install command to be used
_check_distro() {
    _distro="$(grep -i "^ID=" /etc/os-release | cut -d'=' -f2)"
    case $_distro in
        "debian" | "ubuntu" | "linuxmint" | "pop" | "kali")
            install_pkg="sudo apt install"
            ;;
        "arch")
            install_pkg="sudo pacman -S"
            ;;
        "centos" | "fedora")
            install_pkg="sudo dnf...."
            ;;
        "nixos")
            install_pkg="echo really want to"
            ;;
    esac
}
_check_distro

# array of dependency tools required
_dependency=(
    "fastp"
    "bwa"
    "bowtie2"
    "samtools"
    "lolcat"
)

# check dependency package and install it
_check_pkgs() {
    for pkg in "${_dependency[@]}"; do
        if [ ! -x "$(command -v "$pkg")" ]; then
            echo "Package $pkg not installed, Installing it now."
            echo "Please provide the sudo password to run '$install_pkg $pkg' for installing the package:"
            $install_pkg "$pkg"
        fi
    done
    echo "All packages are ready!"
}
_check_pkgs

pre_process() {
    echo "fastp on here"
}





while [ $# -gt 0 ]; do
	  case "$1" in
	      -h | --help)
		        help_info
		        ;;
        -i | --input)
            input_dir=$1
            ;;
        -p)
            pre_process
            ;;
	  esac
done
