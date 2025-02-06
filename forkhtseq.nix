{
  lib,
  buildPythonPackage,
  fetchFromGitHub,
  fetchpatch,
  swig,
  cython,
  matplotlib,
  numpy,
  pandas,
  pysam,
  setuptools,
  pytestCheckHook,
  nix-update-script,
}:
buildPythonPackage rec {
  pname = "fork-htseq";
  version = "1.0";
  pyproject = true;

  src = ./htseq;
  #   fetchurl {
  #   owner = "htseq";
  #   repo = "htseq";
  #   rev = "release_${version}";
  #   hash = "sha256-7ocrmuj9LOtPz9XbI5rKGcdE5JbFz/pZh00Nie65XxE=";
  # };

  nativeBuildInputs = [ swig ];

  build-system = [
    cython
    numpy
    pysam
    setuptools
  ];

  dependencies = [
    numpy
    pysam
  ];

  optional-dependencies = {
    htseq-qa = [ matplotlib ];
  };

  pythonImportsCheck = [ "HTSeq" ];

  nativeCheckInputs = [
    pandas
    pytestCheckHook
  ] ++ optional-dependencies.htseq-qa;

  preCheck = ''
    rm -r src HTSeq
    export PATH=$out/bin:$PATH
  '';

  passthru.updateScript = nix-update-script {
    extraArgs = [
      "--version-regex"
      "release_(.+)"
    ];
  };

  doCheck = false;

  meta = with lib; {
    homepage = "https://htseq.readthedocs.io/";
    description = "Framework to work with high-throughput sequencing data";
    maintainers = with maintainers; [ unode ];
  };
}
