{
  description = "A devShell with rust tools";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    fenix.url = "github:nix-community/fenix";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { nixpkgs, fenix, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ fenix.overlays.default ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        fx = fenix.packages.${system};
        rust-toolchain = fx.combine [
          fx.latest.cargo
          fx.latest.rustc
          fx.latest.rust-analyzer
          fx.latest.clippy
          fx.latest.rustfmt
          fx.latest.rust-src
          fx.latest.miri
        ];
      in
      with pkgs;
      {
        devShells.default = mkShell {
          buildInputs = [
            cargo-flamegraph
            cargo-fuzz
            cargo-nextest
            cargo-outdated
            cargo-show-asm
            taplo
          ];

          nativeBuildInputs = [
            rust-toolchain
          ];
        };
      }
    );
}
