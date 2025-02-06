{
lib,
buildPythonPackage,
fetchFromGitHub,
setuptools,
numpy,
pysam,
}:

buildPythonPackage rec {
  pname = "cpat";
  version = "3.0.5";
  format = "pyproject";

  src = fetchFromGitHub {
    owner = "liguowang";
    repo = pname;
    rev = "refs/tags/v${version}";
    hash = "sha256-PH0u9ePn6nA6/w+dnZpJW+0vXdd2i06fO24R21Vfgo4=";
  };

  propagatedBuildInputs = [
    setuptools
    numpy
    pysam
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
