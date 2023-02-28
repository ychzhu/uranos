# %%
"""
Generate batch files for multiple URANOS jobs.
1. Walks recursively through a folder,
2. Appends to a list of scenarios every folder that includes Uranos.cfg
3. Equally split the list into a number x of batch files for easy parallelizazion.
4. Double click each batch file to run URANOS
"""
import os, fnmatch
import numpy as np

# %%
# Set the folder that includes folders of scenarios
work_on_folder = "my_project"
# Set URANOS executable
uranos_path = 'C:/Users/username/URANOS-v107/UranosGUI.exe'
# Split scenarios into X batch files
split_into_parts = 10
scenarios = []

# %%
def generate_batch_file(directory, filePattern):
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            scenarios.append('%s noGUI %s' % (uranos_path, filepath))
            #print(filepath, flush=True)
    #for c in range(copies):
    i = 0
    for scenario in np.array_split(scenarios, split_into_parts):
        i += 1
        filepath = 'uranos-batchrun-%s-%02d.bat' % (work_on_folder, i)
        print(filepath, flush=True)
        with open(filepath, "w") as f:
            f.write("\n".join(scenario))
        
# %%
generate_batch_file(work_on_folder, "*ranos.cfg")

