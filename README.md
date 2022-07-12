# URANOS - Ultra Rapid Neutron-Only Simulation
[![URANOS at Zenodo](https://img.shields.io/static/v1?label=Code&message=10.5281/zenodo.6578668&color=blue)](https://doi.org/10.5281/zenodo.6578668) [![URANOS at the Geoscientific Model Development](https://img.shields.io/static/v1?label=Manuscript&message=10.5194/gmd-2022-93&color=yellow)](https://doi.org/10.5194/gmd-2022-93)

![splashScreenUranos](https://user-images.githubusercontent.com/106179070/170102368-93e66f49-12ab-44a9-860a-2bd1977f715c.jpg)

URANOS is a Monte Carlo toolkit specifically tailored for environmental sciences. It can be used to simulate cosmogenic neutron radiation near the Earth's surface and to study its response to environmental factors, such as soil water content, snow, or biomass.


Key features:

- stand-alone executable,
- parallelization (can run performantly in multiple parallel instances), 
- graphical user interface (optional) for the command line.
- 3D voxel engine to allow the generation of complex environments

## Use, discuss, contribute

In the **uranos-users mailing list**, users can ask questions, discuss, and contribute to the development of URANOS user community. [Subscribe now!](https://www.ufz.de/index.php?en=41538)

## Cite as

- **The software code** is available [on GitHub](https://github.com/mkoehli/uranos) and has been released by [Zenodo, doi:10.5281/zenodo.6578668](https://doi.org/10.5281/zenodo.6578668).
- **A manuscript** on URANOS has been published as a preprint in the journal *Geoscientific Model Development*:
    > Köhli, M., Schrön, M., Zacharias, S., and Schmidt, U.: URANOS v1.0 – the Ultra Rapid Adaptable Neutron-Only Simulation for Environmental Research, Geosci. Model Dev. Discuss. [preprint], [doi:10.5194/gmd-2022-93](https://doi.org/10.5194/gmd-2022-93), *in review*, 2022. 

## Publications using URANOS

- [See the list of publications using URANOS](doc/PUBLICATIONS.md)

## Usage

### Prerequisites

- [Visual Studio 2019 community version](https://my.visualstudio.com/Downloads?q=visual%20studio%202019&wt.mc_id=o~msft~vscom~older-downloads) (prerequisite for ROOT6). Install with the development option: `C++ application`
- [ROOT 6.22.08](https://root.cern/download/root_v6.22.08.win32.vc16.exe) (prerequisite for running URANOS)
- Basic libraries (possibly missing on some systems):
    - Microsoft Visual C++ 2015-2022 redistributional package
    - Windows 10/11 SDK (With C++ development files)

### Preparation

1. Download the zipped binary package and the ENDF data zip file. 
2. Unzip the URANOS files into a folder of your choice.
3. Unzip the ENDF data files, which also contains the `Input Spectrum Calibration File`, into a folder of your choice. 
4. Run `URANOSGUI.exe`.
5. Insert the full path and filename of the `IncomingSpectrum.root` into `Input Spectrum Calibration File` under URANOS' *Setup* tab. 
6. Insert the full path to the ENDF files into `Cross Section Folder` under URANOS' *Setup* tab. 
7. Set a `work directory` where the configuration files of your scenario are to be stored. Mind a trailing slash.
8. Set an `output directory` where the simulation results will be stored.  Mind a trailing slash.

### First Run

In order to perform a simulation, a physics environment has to be configured. For most settings the basic default values are already set. Press `Load Minimal Config` to generate a standard layer setup. Then press `Simulate`. If the simulation has started, the GUI can be closed, the calculation will be executed in a separate command line window.

### Individual Scenarios

In order to configure your own simulation:
1. Define a layer structure of pre-defined materials according to the examples on the website or the minimal config. 
2. Press `Save` to store the geometry configuration in the `work folder`, and `Load` to load it to the GUI. 
3. Put a PNG image file into the `work folder` and name it by the number of the layer you want to define, for example `6.png` for the layer number 6.
 By pressing `Use Layer Maps` URANOS will search for these PNG files in the work folder.
 
 *Note:* The PNG files are a convenient way to define your input material composition at the horizontal scale. They have to be in grayscale (or similar RGB colors) and in an aspact ratio of 1:1 (quadratic). Predefined Materials are encoded by these grayscale values, see the file `Material Codes.txt` in the URANOS folder. The PNG will be stretched to the full domain size and each pixel will be extruded in the layer to a 3D pixel (voxel) of the given material.

 
