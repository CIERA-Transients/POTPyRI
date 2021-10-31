#Download and sort Keck data for KOA
#Requires full path name to wget bash script downloaded from KOA
#and new folder name to put data (e.g. MOSFIRE_20191123)

import os
import sys
import subprocess
import tarfile
import shutil
import glob

full_download_script = sys.argv[1]
data_path = os.path.dirname(full_download_script)
download_script = os.path.basename(full_download_script)
os.chdir(data_path)

subprocess.call(['bash',download_script])

tar_files = glob.glob('*.tar')
for tar in tar_files:
    tar_f = tarfile.open(tar)
    tar_f.extractall(path=data_path)
    tar_f.close()
    os.remove(tar)

KOA_ID = download_script.split('_')[1].split('.')[0]
files = glob.glob('KOA_'+KOA_ID+'/*/*/*/*.fits.gz')
new_dirname = sys.argv[2]
if not os.path.exists(new_dirname):
    os.mkdir('./'+new_dirname)
for f in files:
    shutil.move(f,'./'+new_dirname)
subprocess.call(['rm','-r','KOA_'+KOA_ID])
