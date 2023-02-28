#!python
# URANOS Monitor
#    Usage: just doubleclick
# 1. Reads all *.txt files in the same folder
# 2. Prints their first lines
# 3. Presents summary statistics
# 4. Uploads the output to an FTP server
# 5. Repeats every few seconds

# CONFIG HERE
# Refresh monitor window every x seconds
refresh_every_seconds = 30
# FTP settings
FTP_server = ''  
FTP_user = ''      
FTP_pswd = ''
FTP_folder = ''

#######################
import time
from glob import glob
import re
import sys
import os
import numpy as np
import ftplib
from io import BytesIO

def clear_screen():
    # use this when opening with python console:
    cls = lambda: os.system('cls')
    cls()
    # Use this in Git Bash:
    #print("\033c")
    
all_item  = 0 
all_items = 0
list_time_per_secenario = []

# Infinite loop
while(1):
    # Clear screen
    clear_screen()
    
    files = sorted(glob('*.txt'))
    all_item  = 0 
    all_items = 0
    content = []

    # Loop all files
    for file in files:
        
        # Read
        with open(file, 'r') as f:
            line = f.read()
        print('.', end='')
        
        # Read number of items
        match_item = re.match(r'^\s*(\d+)\/(\d+)\s', line)
        print('.', end='')
        if match_item:
            print('!', end='\r')
            
            item = float(match_item.group(1))-1
            items = float(match_item.group(2))
            if item==items-1 and '100.0%' in line:
                item += 1
            all_item += item 
            all_items += items
            
            if '100.0%' in line:
                # Grep total runtime from finished scenarios
                match_time = re.search(r'\(\s*(\d+\.\d)h\)', line)
                if match_time:
                    list_time_per_secenario.append(float(match_time.group(1)))
                
            
            progress_pc = item/items *100
            #progress_bar = '>' * int(progress_pc/10)
            out = '%3.0f%% %s: %s' % (progress_pc, file, line)
            content.append(out)
            print(out)
            
    processes = len(files)
    if all_items == 0:
        # If nothing yet to monitor
        out = 'Monitoring %d processes ...' % processes
        content.append(out)
        print(out)
            
    else:
        # Progress bar
        progress_pc = int(all_item/all_items*100)
        progress_bar = '>' * int(progress_pc*1.02)
        out = '|%s| %3.0f%%' % (progress_bar.ljust(102), progress_pc)
        content.append(out)
        print(out)

        # Summary statistics
        avg_time = np.nanmean(list_time_per_secenario)                if len(list_time_per_secenario)>0 else np.nan
        total_eta = (all_items - all_item) * avg_time /24 /processes  if len(list_time_per_secenario)>0 else np.nan
        out = 'Monitoring %d processes: %4.0f/%d scenarios finished, average duration: %3.1fh, total ETA: %4.1f days.' % (processes, all_item, all_items, avg_time, total_eta)
        content.append(out)
        print(out)
        
    
    # save output to monitor.log   
    with open('monitor.log', 'w') as f:
        f.write("\n".join(content))

    if FTP_server:
        # Upload to FTP
        try:
            print(' Uploading to %s:' % FTP_server, end='')
            remote = ftplib.FTP(FTP_server)
            remote.login(FTP_user, FTP_pswd)
            
            local_file = open('monitor.log', "rb")
            remote_file = '%s/monitor.log' % FTP_folder
            print('%s' % remote_file, end='')
            remote.storbinary("STOR %s" % remote_file, local_file)
            print(' ...OK (%.0f kB)' % (remote.size(remote_file)/1024), end='')
            xcmd = remote.sendcmd('SITE CHMOD 644 %s' % remote_file)
            #print(' OK (%s).' % xcmd)

        except ftplib.all_errors as e:
            print(str(e))
            print('= Sorry.')
            sys.exit()

        remote.close()
    
    # Flush output and sleep
    sys.stdout.flush()
    time.sleep(refresh_every_seconds)
    
