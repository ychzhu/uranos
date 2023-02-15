{ stdenv
, lib
, qtbase
, wrapQtAppsHook
}:

stdenv.mkDerivation rec {
  pname = "uranos";
  version = "v1.0.9";

  src = ./.;

  buildInputs = [ qtbase ];
  nativeBuildInputs = [ wrapQtAppsHook ];

}
