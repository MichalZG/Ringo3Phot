#! /usr/bin/env python
import glob
import os
import shutil
import fnmatch

ext = '*.fits'
outDirName = 'results'
extToRm = ['csv', 'pdf', 'fits']
workPath = os.getcwd()
script_path = os.path.dirname(os.path.realpath(__file__))

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        dirs[:] = [d for d in dirs if (d.isdigit() and d.__len__() == 8)]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(workPath, root, name))
    return result


try:
    os.mkdir(os.path.join(workPath, outDirName))
except OSError:
    # unlink all files in results dir
    files = glob.glob(os.path.join(workPath, outDirName, '*'))
    for f in files:
        if os.path.islink(f):
            os.unlink(f)
	elif f.split('.')[-1] in extToRm:
	    os.remove(f)

    shutil.copy(os.path.join(script_path, 'ephem.csv'), os.path.join(workPath, outDirName))

fileList = find(ext, workPath)
for fileToLink in fileList:
    fileName = fileToLink.split('/')[-1]
    os.symlink(fileToLink,
               os.path.join(workPath, outDirName, fileName))
print 'Done'
