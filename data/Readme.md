# Information about data

- The ENDF files contain cross sections required to run URANOS and evaluate the interactions with materials. The ENDF data folder is to be set as the 'Cross Section Folder' in URANOS in the Setup tab. 
- The IncomingSpectrum.root is to be declared as the Input Spectrum Calibration File in URANOS in the Setup tab.
- The parameter files, included in the ENDF zip file are required the evaluate the Sato 2016 incoming spectrum.
- The response functions can be used to emulate different detectors in URANOS when selecting the 'Physics Model' for scoring. By setting a path to the respective file URANOS uses the response function of the named model. Otherwise the standard detector geometry with 25 mm moderator thickness is used.
