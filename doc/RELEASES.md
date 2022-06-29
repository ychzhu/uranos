# URANOS Release Notes

## URANOS v1.00a (May 24, 2022)

### New installation guidelines

Changed the program architecture to more recent versions of compilers and dependencies, also to prepare the Linux release.
URANOS now requires:
- Visual Studio 2019 community. Install with the development option "C++
application" (approx. 7 GB). Unfortunately, this installation is required to
run ROOT 6. Link: <https://my.visualstudio.com/Downloads?q=visual%20studio%202019&wt.mc_id=o~ms
ft~vscom~older-downloads>
- ROOT (Installer provided on the website, approx. 100 MB) 6.22.08 (a virus
scan warning might appear here). Link: <https://root.cern/download/root_v6.22.08.win32.vc16.exe>
(In case you do not need to run older versions of URANOS it is recommended
to uninstall ROOT before installing the new version)

### Backward-incompatible changes:
- Due to the change in the high energy model, the intensity might scale differently in absolute numbers compared to v0.99.

### GUI improvements:
- "Save Configuration" button added.
