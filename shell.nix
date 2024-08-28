with import <nixpkgs> {};

let
  pypkgs = (python312.withPackages ((ps: with ps; [
    venvShellHook
    pip
    biopython
    (callPackage ./cpat.nix {})
    numpy
    pysam
    (callPackage ./forkhtseq.nix {})
    (callPackage ./gffutil.nix {})
  ])));


  rpkgs = with rPackages; [
    lintr styler
    # languageserver
    # data_table
    # dplyr
    ggplot2 ggpubr
    tidyverse # includes 9 packages that are down below
    GDCRNATools topGO ALL
    Biostrings
    BiocManager
    DESeq2 # gene exp analysis
    # rmarkdown
    # add more packages (https://search.nixos.org)
    # from CRAN or Bioconductor
  ];

  in

pkgs.mkShell {

  nativeBuildInputs = [ pkgs.bashInteractive ];

  # EnvVars = The thung
  # NIX_LD_LIBRARY_PATH = lib.makeLibraryPath [
  # pkgs
  # ];
  # NIX_LD = lib.fileContents "${stdenv.cc}/nix-support/dynamic-linker";

  buildInputs = with pkgs; [
    shellcheck bash-language-server shfmt
    fastp samtools bwa bowtie2
    subread star
    vscodium bedtools
    yq
    pypkgs

    (rWrapper.override {
      packages = [
        rpkgs
      ];
    })

  ];

  shellHook = ''
    echo "Environment to detect lncRNA is loaded!"
  '';

}
