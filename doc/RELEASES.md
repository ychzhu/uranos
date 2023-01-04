# URANOS Release Notes

## URANOS version 1.05 (January 4, 2023)

### Changelog
- bugfixes
- added material: 231 = Polyvinyl chloride (1.45 g/cm^3)


## URANOS version 1.04 (January 2, 2023)

### Changelog
- additional high resolution track output added
- URANOS CentOS 7 64bit version release with QT 5.9.7 and ROOT 6.22.08.


## URANOS version 1.03 (December 30, 2022)

### Changelog
- URANOS Linux version release. 
- URANOS has currently been tested on Ubuntu 20 64bit with QT 5.12.8 and ROOT 6.22.08.
- Compiled binaries are now available for Ubuntu 20. For running under Linux, please also update the ENDF data.


## URANOS version 1.02 (December 26, 2022)

### Changelog
- bugfixes for command line options


## URANOS version 1.01 (December 5, 2022)

### Changelog
- smaller bugfixes
- updated QCustomPlot to version 2.1.1


## URANOS version 1.00 beta (September 22, 2022)

### Changelog
- bugfixes for material settings


## URANOS version 1.00a (May 24, 2022)

### New installation guidelines

Changed the program architecture to more recent versions of compilers and dependencies, also to prepare the Linux release.
URANOS now requires:
- Visual Studio 2019 community. Install with the development option "C++
application" (approx. 7 GB). Unfortunately, this installation is required to
run ROOT 6. [Link](https://my.visualstudio.com/Downloads?q=visual%20studio%202019&wt.mc_id=o~msft~vscom~older-downloads)
- ROOT (Installer provided on the website, approx. 100 MB) 6.22.08 (a virus
scan warning might appear here). [Link](https://root.cern/download/root_v6.22.08.win32.vc16.exe)
(In case you do not need to run older versions of URANOS it is recommended
to uninstall ROOT before installing the new version)

### Backward-incompatible changes:
- Due to the change in the high energy model, the intensity might scale differently in absolute numbers compared to v0.99.

### GUI improvements:
- "Save Configuration" button added.
