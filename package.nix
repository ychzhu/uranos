{ stdenv
, lib
, qtbase
, root
, tbb
, wrapQtAppsHook
}:

stdenv.mkDerivation rec {
  pname = "uranos";
  version = "1.0.9";

  src = ./.;

  buildInputs = [ qtbase root tbb ];
  nativeBuildInputs = [ wrapQtAppsHook ];

  configurePhase = ''
    runHook preConfigure
    qmake PREFIX=$out
    runHook postConfigure
  '';
}
