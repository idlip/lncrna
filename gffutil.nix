{
lib,
buildPythonPackage,
fetchFromGitHub,
setuptools,
pyfaidx,
    argh,
    argcomplete,
    simplejson
}:

buildPythonPackage rec {
  pname = "gffutils";
  version = "0.13";
  format = "setuptools";

  src = fetchFromGitHub {
    owner = "daler";
    repo = pname;
    rev = "refs/tags/v${version}";
    hash = "sha256-aoTThkYyX80TV8Ahnl+jtPFU0AfUjitT+JmrVeY6EZE=";
  };

  propagatedBuildInputs = [
    setuptools
    pyfaidx
    argh
    argcomplete
    simplejson
  ];

  # The test needs MuJoCo that is not free library.
  doCheck = false;

  # pythonImportsCheck = [ "gym" ];

  meta = with lib; {
    description = "Toolkit for developing and comparing your reinforcement learning agents";
    homepage = "https://gym.openai.com/";
    license = licenses.mit;
    maintainers = with maintainers; [ hyphon81 ];
  };
}
