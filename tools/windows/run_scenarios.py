#!python
# URANOS Scenario list runner
#   Usage: python run_scenarios.py -f uranos_scenario_list-01.txt
# 1. Reads the given file which contains the folder names for the scenarios
# 2. For each scenario, run URANOS
# 3. Presents only important information from URANOS output
# 4. Logs last line to a monitring folder

# CONFIG HERE
# Root folder 
root_folder_with_slashes = 'C:/path/to/work/'
# Folder in which to save monitoring information
monitor_folder = 'monitor'
# URANOS settings
URANOS_app = 'C:/Users/username/URANOS-107/UranosGUI.exe'
URANOS_par = 'noGUI'
URANOS_cfg = 'Uranos.cfg'

#######################
import subprocess
import re
from time import time 

# Parse arguments
from argparse import ArgumentParser
parser = ArgumentParser(
    prog = 'run_scenarios',
    description = 'Batchrun through all scenarios in a txt file')
parser.add_argument(
    "-f", "--file", dest="file_scenarios_list",
    help="Read scenario list from FILE", metavar="FILE",
    default='uranos-scenarios-distzone_9-01.txt')
parser.add_argument(
    "-r", "--root", dest="scenario_root_with_slashes",
    help="Scenario root path with slashes",
    default=root_folder_with_slashes)
args = parser.parse_args()

monitor_file = "%s/%s" % (monitor_folder, args.file_scenarios_list)


# Execute URANOS
def execute(cmd):
    """
    Thanks to https://stackoverflow.com/a/4417735/2575273
    """
    #print('Executing %s' % cmd, flush=True)
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

# Read scenario list
with open(args.file_scenarios_list) as file:
    scenarios = [line.strip() for line in file]

# Loop through scenarios
items = len(scenarios)
item = 0
for scenario in scenarios:
    # Generate command
    command = [URANOS_app, URANOS_par, args.scenario_root_with_slashes +'/'+ scenario +'/'+ URANOS_cfg]    
    out = '%4.0f/%d %s: ... ' % (item, items, scenario)
    print(out, end="\r", flush=True)
            
    # Start timer
    item += 1
    time_start = int(time())
    
    # Execute and grep output line by line
    for line in execute(command):
        
        # Match lines
        match_runtime  = re.match(r'Runtime: (\d+) s', line.rstrip())
        match_progress = re.match(r'(\d.+) % completed', line.rstrip())
        if match_runtime:
            # Scenario finished
            out = '%4.0f/%d %s %5.1f%% (%5.1fh) ' % (item, items, scenario.ljust(50), 100, int(match_runtime.group(1))/3600)
            print(out, flush=True)
            with open(monitor_file, 'w') as f:
                f.write(out)
                
        elif match_progress:
            # Scenario running
            progress_pc = float(match_progress.group(1))
            time_now = int(time())
            if progress_pc == 0.0:
                eta = 0
            else:
                eta = ( 100/progress_pc * (time_now-time_start) +time_start -time_now )/3600
            
            out = '%4.0f/%d %s %5.1f%% (-%4.1fh) ' % (item, items, scenario.ljust(50), int(progress_pc*10)/10, eta)
            print(out, end="\r", flush=True)
            
            # Write line to seperate file to faciliate monitoring later
            with open(monitor_file, 'w') as f:
                f.write(out)
        else:
            pass
            #print('%4.0f/%d %s: %s ' % (item, items, scenario, line.rstrip()), end="\r", flush=True)
        
print("\n")

