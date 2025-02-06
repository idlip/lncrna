{ lib
, fetchFromGitHub
, rustPlatform
}:

rustPlatform.buildRustPackage rec {
  pname = "bed2gtf";
  version = "1.9.2";

  src = fetchFromGitHub {
    owner = "alejandrogzi";
    repo = "bed2gtf";
    rev = "refs/tags/v.${version}";
    hash = "sha256-lspGle5729sE6H7W5EF73RwrXKpJP/+U8alWrG6ASww=";
  };

  cargoHash = "sha256-pGJf2pKDuUVqmMGKqloKAymfxjIuBMu2F+aWllrk2+4=";

  # buildInputs = rpathLibs;

  # postInstall = ''
  #   install -D -m 0555 leftwm/doc/leftwm.1 $out/share/man/man1/leftwm.1
  # '';

  # dontPatchELF = true;

  meta = {
    description = "Tiling window manager for the adventurer";
    homepage = "https://github.com/leftwm/leftwm";
    license = lib.licenses.mit;
    platforms = lib.platforms.linux;
    maintainers = with lib.maintainers; [ yanganto ];
    changelog = "https://github.com/leftwm/leftwm/blob/${version}/CHANGELOG.md";
    mainProgram = "leftwm";
  };
}
