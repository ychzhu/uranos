# Information about examples

The examples are taken from actual simulation scenarios and can be used to either simulate the specific topology or used as a base to create own simulations.

The folders within 'examples' contain each 
- the geometry configuration "UranosGeometryConfig.dat", which defines the layer structure and 
- the "URANOS.cfg", which contains the settings for the simulation run from that folder. The up to five folder paths in the top of the file need to be adapted to the reference computer system.
- the ASCII matrices or PNG files whith the material definitions used to extrude the voxel geometry. Please refer to the "MaterialCodes.txt" for the assignment of numbers to specific materials.

The exmpamples for the detector models are designed to be used with the batchrun feature in order to generate the energy-dependent detector response functions.
