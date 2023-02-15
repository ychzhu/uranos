let
  pkgs = import <nixpkgs> {};
in libsForQt5.callPackage ./package.nix
