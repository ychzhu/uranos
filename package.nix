{ stdenv
, lib
, qtbase
, root
, tbb
, unzip
, wrapQtAppsHook
}:

stdenv.mkDerivation rec {
  pname = "uranos";
  version = "1.0.9";

  src = ./.;

  buildInputs = [ qtbase root tbb ];
  nativeBuildInputs = [ wrapQtAppsHook unzip ];

  qtWrapperArgs = [ "--unset LD_LIBRARY_PATH" ];

  configurePhase = ''
    runHook preConfigure
    qmake PREFIX=$out
    runHook postConfigure
  '';

  postInstall = ''
    cd $out/share/uranos
    unzip ENDFdata.zip
    rm ENDFdata.zip
  '';
}
