#!python
# URANOS Check status
#   Usage: python check_uranos.py
# 1. Reads remote monitoring file from FTP
# 2. Prints its content to the console
# (3. Potential improvement: auto-refresh similar to monitor.py, but that woudl occupy CPU)

# CONFIG HERE
# FTP settings
FTP_server = ''
FTP_user = ''
FTP_pswd = ''
FTP_file = ''

######################
import sys
import ftplib
from io import BytesIO

try:
    print('| Connecting... ', end='')
    remote = ftplib.FTP(FTP_server)
    print('Server OK.. ', end='')

    remote.login(FTP_user, FTP_pswd)
    print('Login OK.. ', end='')

    #remote.prot_p() # Enables secure communication, prerequisit for some providers
    #print('Encryption OK.')

    print('File... ', end='')
    memory = BytesIO()
    remote.retrbinary("RETR %s" % FTP_file, memory.write)
    print('%.0f kB. ' % (memory.getbuffer().nbytes/1024), end='')
    
except ftplib.all_errors as e:
    print(str(e))
    print('= Sorry.')
    sys.exit()

remote.close()
print('Closed.')

memory.seek(0)
for l in memory:
    print(l.decode("utf-8").rstrip())

