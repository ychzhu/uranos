{ stdenv
, lib
, qtbase
, root
, tbb
, wrapQtAppsHook
}:

stdenv.mkDerivation rec {
  pname = "uranos";
  version = "v1.0.9";

  src = ./.;

  buildInputs = [ qtbase root tbb ];
  nativeBuildInputs = [ wrapQtAppsHook ];

  buildPhase = ''
    runHook preBuild
    cd src
    qmake
    make
    runHook postBuild
  '';
}
