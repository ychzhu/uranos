#!python
# URANOS Scenario list generator
#   Usage: python generate_scenario_list.py -s 10
# 1. Walks recursively through folder (-w name) 
# 2. Collect folders containing Uranos.cfg but not *.root (i.e., unfinished scenarios)
# 3. Write all those folder names evenly distributed into x files (-s number) 

# CONFIG HERE
# Root folder 
root_folder_with_slashes = 'C:/path/to/work/'
# Distribute all scenarios across x files, or: -s 90
split_scenarios_per_file_default = 10
# Display overall file name of this endeavor, or: -w my_study
display_name_default = 'my_scenarios'

#######################
import os, fnmatch
import numpy as np
from glob import glob

from argparse import ArgumentParser
parser = ArgumentParser(
    prog = 'generate_scenario_list',
    description = 'Generate a (splitted) list of untouched scenarios, store into txt files.')
parser.add_argument("-s", "--split",
    dest="split_into_parts", default=split_scenarios_per_file_default,
    help="Split scenario list into given parts/files", type=int)
parser.add_argument("-r", "--root",
    dest="scenario_root_with_slashes", default=root_folder_with_slashes,
    help="Scenario root path with slashes")
parser.add_argument("-w", "--work",
    dest="work_on_folder", default='', 
    help="Folder name top level")
args = parser.parse_args()

scenarios = []

if args.work_on_folder != '':
    display_name = args.work_on_folder
else:
    display_name = display_name_default

num_completed = 0
num_todo = 0

# Walks through folder 
def generate_scenario_list(directory, filePattern):
    global num_completed, num_todo
    print('Walking through %s ...' % os.path.abspath(directory))
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        # Only consider folders containing Uranos.cfg
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            
            # Only consider folder not containing *.root
            root_files = glob(os.path.join(path, '*.root'))
            if len(root_files)>0:
                num_completed += 1
            else:
                scenarios.append('%s' % path.replace('\\','/').replace(args.scenario_root_with_slashes, ''))
                num_todo += 1
            print('%d completed, %d todo' % (num_completed, num_todo), end='\r', flush=True)
    #for c in range(copies):
    print("\n")
    
    # Split evenly across x files
    i = 0
    for scenario in np.array_split(scenarios, args.split_into_parts):
        i += 1
        filepath = 'uslist-%s-%02d.txt' % (display_name, i)
        print("> %s" % filepath, flush=True)
        with open(filepath, "w") as f:
            f.write("\n".join(scenario))
        
# %%
generate_scenario_list(args.work_on_folder, "*ranos.cfg")
#print('%d completed, %d todo' % (num_completed, num_todo))
