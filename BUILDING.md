# Building URANOS

## Linux

First, clone the URANOS repository and `cd` into it:
```console
$ git clone https://github.com/mkoehli/uranos
$ cd uranos
```

### Building with Nix
Currently, the easiest way to build URANOS from scratch is via the provided Nix derivations. If you already have Nix installed, run:
```console
$ nix-build
```
Nix will automatically fetch all dependencies for you. The resulting program will be located at `result/bin/uranos-gui`.

To install URANOS on your machine, run:
```console
$ nix-env --install -f default.nix
```


Nix can be installed on your machine on almost any Linux distribution with the following command:
```console
$ sh <(curl -L https://nixos.org/nix/install) --daemon
```
For more information on Nix, see [https://nixos.org](https://nixos.org).

### Without Nix
If you do not want to use Nix, you will need to manually install the following dependencies:

* [Qt 5](https://qt.io)
* [ROOT](https://root.cern.ch/)
* libtbb (probably available through your distribution package manager)

To build URANOS, run the following commands inside the URANOS source tree:
```console
$ qmake
$ make
```

To install URANOS, run:
```console
$ sudo make install
```

## Windows

(work in progress)
