#!python
# URANOS Setting changer
#   Usage: just doubleclick
# 1. Walks recursively through folder
# 2. Replace parts in Uranos.cfg by given replacements

# CONFIG HERE
# Walk from here recursesively:
work_on_folder = '.'
# Put your replacements here:
Replacements = [
    [
        "C:/path/to/work/",
        "C:/path/to/other/work/"
    ],[
        "URANOS-v101",
        "URANOS-v107"
    ]
]

#######################
import os, fnmatch

# Walk through folder and treat Uranos.cfg
def findReplace(directory, filePattern):
    i = 0
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            i += 1
            filepath = os.path.join(path, filename)
            print("%6.0f: %s" % (i, filepath), end="\r", flush=True)

            # Read Uranos.cfg
            with open(filepath) as f:
                s = f.read()

            # Replace
            for R in Replacements:
                s = s.replace(R[0], R[1])

            # Save Uranos.cfg
            with open(filepath, "w") as f:
                f.write(s)
    print("\n")
# %%
findReplace(work_on_folder, "*ranos.cfg")
