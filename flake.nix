{
  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-unstable;

  outputs = { self, nixpkgs, ... }:
  let
    pkgs = import nixpkgs { system = "x86_64-linux"; };
  in {
    packages.x86_64-linux.uranos = pkgs.libsForQt5.callPackage ./package.nix { };
    packages.x86_64-linux.default = self.packages.x86_64-linux.uranos;
  };
}
