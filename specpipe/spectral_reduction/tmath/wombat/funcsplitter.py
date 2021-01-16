import sys
filename=input('Enter filename to split: ')
filename=filename.strip()

try:
    infile=open(filename,'r')
except IOError:
    print('Bad filename')
    sys.exit(1)

inifile='__init__.py'
ini=open(inifile,'w')
    
for line in infile:
    if (line[0:3] == 'def'):
        try:
            outfile.close()
        except NameError:
            pass
        parts=line.split()
        loc=parts[1].find('(')
        ofilename=parts[1][0:loc]+'.py'
        ini.write('from .'+parts[1][0:loc]+' import '+parts[1][0:loc]+'\n')
        outfile=open(ofilename,'w')
    outfile.write(line)

    

outfile.close()
infile.close()
ini.close()
