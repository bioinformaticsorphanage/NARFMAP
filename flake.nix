{
  description = "Fork of DRAGMAP, the first of it's name, merger of PRs, cutter of releases, takeoverer of worlds";

  inputs.nixpkgs.url = "https://flakehub.com/f/NixOS/nixpkgs/0.1.*.tar.gz";
  inputs.flake-parts.url = "github:hercules-ci/flake-parts";
  inputs.treefmt-nix.url = "github:numtide/treefmt-nix";

  outputs = inputs @ {flake-parts, ...}:
    flake-parts.lib.mkFlake {inherit inputs;} {
      imports = [
        inputs.treefmt-nix.flakeModule
      ];
      flake = {
        # TODO Put your original flake attributes here.
      };
      systems = [
        # systems for which you want to build the `perSystem` attributes
        "x86_64-linux"
        # TODO "aarch64-linux" "x86_64-darwin" "aarch64-darwin"
      ];
      perSystem = {
        config,
        pkgs,
        ...
      }: {
        #
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            boost
            gnumake
            gtest
            zlib
          ];

          env = {
            BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
            BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";
            GTEST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.gtest}/include";
            GTEST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.gtest}/lib";
            GTEST_ROOT = "${pkgs.gtest}";
            LD_LIBRARY_PATH = "${pkgs.lib.getLib pkgs.gtest}/lib";
            LIBCLANG_PATH = "${pkgs.lib.getLib pkgs.llvmPackages.libclang.lib}/lib";
          };
        };

        packages.default = config.packages.dragmap;

        packages.dragmap = pkgs.stdenv.mkDerivation {
          pname = "NARFMAP";
          version = "1.4.3";

          src = ./.;

          buildInputs = [
            pkgs.boost
            pkgs.gnumake
            pkgs.gtest
            pkgs.zlib
          ];

          buildPhase = ''
            make
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp build/release/compare $out/bin
            cp build/release/dragen-os $out/bin
          '';

          env = {
            HAS_GTEST = 0;
            # FIXME nix and static compiling doesn't seem to be working
            # STATIC = 1;
            BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
            BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";
            GTEST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.gtest}/include";
            GTEST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.gtest}/lib";
            GTEST_ROOT = "${pkgs.gtest}";
            LD_LIBRARY_PATH = "${pkgs.lib.getLib pkgs.gtest}/lib";
            LIBCLANG_PATH = "${pkgs.lib.getLib pkgs.llvmPackages.libclang.lib}/lib";
          };

          meta = with pkgs.lib; {
            homepage = "https://github.com/bioinformaticsorphanage/NARFMAP";
            description = ''
              The Dragen mapper/aligner Open Source Software.
            '';
            licencse = licenses.gpl3;
            platforms = with platforms; linux;
            maintainers = [maintainers.edmundmiller];
          };
        };
        treefmt.config = {
          projectRootFile = "flake.nix";

          programs.alejandra.enable = true;
          programs.clang-format.enable = true;
          programs.deadnix.enable = true;
          programs.prettier.enable = true;
          programs.prettier.settings = {
            editorconfig = true;
            embeddedLanguageFormatting = "auto";
          };
        };
      };
    };
}
