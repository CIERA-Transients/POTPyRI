#Download and sort Gemini data
#Requires target name, instrument, download link 
#and gemini archive cookie (required for priority data)

import os
import sys
import subprocess
import tarfile
import shutil
import glob

target = sys.argv[1]
tel = sys.argv[2]
download_link = sys.argv[3]
gemini_cookie = sys.argv[4]

subprocess.run(['curl','-o','/data/Gemini/'+tel+'/'+target+'_sci.tar',download_link,'-H',
                    'Cookie: gemini_archive_session='+gemini_cookie])
subprocess.run(['curl','-o','/data/Gemini/'+tel+'/'+target+'_cal.tar',download_link+'#','-H',
                    'Cookie: gemini_archive_session='+gemini_cookie])

tar_files = glob.glob('/data/Gemini/'+tel+'/'+target+'*.tar')
if not os.path.exists('/data/Gemini/'+tel+'/'+target):
    os.mkdir('/data/Gemini/'+tel+'/'+target)
for tar in tar_files:
    tar_f = tarfile.open(tar)
    tar_f.extractall(path='/data/Gemini/'+tel+'/'+target+'/')
    tar_f.close()
    os.remove(tar)