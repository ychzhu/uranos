# URNAOS
URANOS - the Ultra Rapid Neutron-Only Simulation is a Monte Carlo toolkit specifically tailored for Environmental Sciences

# Run Notes
URANOS v1.00

URANOS runs as a stand-alone executable. No installation is necessary. 
It can also run in several instances. The GUI itself is the User Interface which allows to prepare the simulation. If the simulation has started, the GUI can be closed, the calculation will be executed in a separate command line window.

# Software/Windows Prerequisites:

- Microsoft Visual C++ 2015-2022 redistributional package
- Windows SDK (With C++ development files)
- ROOT (MSI Installer provided on the website, approx. 100 MB) 6.22.08 (a virus scan warning might appear here)

# Usage
Run Preparation:
- Unzip the Uranos files into a folder of your choice
- Run URANOSGUI.exe

Preparing the physics environment:
- Copy the IncomingSpectrum.root into a folder of your choice and provide the full path, including the filename, in the "setup" ribbon in URANOS for the "Input Spectrum Calibration File"
- Unzip the ENDF files into a folder of your choice and provide the full path in the "setup" ribbon in URANOS for the "Cross Section Folder"

Additionally:
- Set a "work directory" where configurations are stored. 
- Set an "output directory" where the simulation results will be stored. 
- Don't forget to end all folders by a slash or backslash

First Run:
In order to perform a simulation, a physics environment has to be configured.
For most settings default values are set. Press "Load Minimal Config" to generate a standard layer setup. Then press "Simulate".

In order to configure your own simulation:
- define a layer structure of pre-defined materials according to the examples on the website or the minimal config. 
- press "Save" to store the geometry configuration in the "work folder", "Load" will load it to the GUI 
- put a png file into the "work folder" and name it by the number of the layer you want to define, for example "6.png" for the layer number 6.
 By pressing "Use Layer Maps" URANOS will search for these png files in the "work folder".
 The png files have to be grayscale and in an aspact ratio of 1:1 (quadratic). Predefined Materials are encoded by these grayscale values, see the file "Material Codes" in the URANOS folder.
 The png will be stretched to the full domain size and each pixel will be extruded in the layer to a 3D pixel (voxel) of the given material.



