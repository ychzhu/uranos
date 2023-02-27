let
  pkgs = import <nixpkgs> {};
in pkgs.libsForQt5.callPackage ./package.nix { }
