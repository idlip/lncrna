with import <nixpkgs> {};

let
  rpkgs = with rPackages; [
    # lintr
    # languageserver
    # data_table
    # dplyr
    ggplot2
    # tidyverse # includes 9 packages that are down below
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
    subread

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
