#! /usr/bin/env python
import glob
import os
import fnmatch


ext = '*.fits'
outDirName = 'results'
scriptPath = os.path.curdir


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path, topdown=True):
        dirs[:] = [d for d in dirs if (d.isdigit() and d.__len__() == 8)]
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


try:
    os.mkdir(os.path.join(scriptPath, outDirName))
except OSError:
    filesToRemove = glob.glob(os.path.join(scriptPath, outDirName, ext))
    for f in filesToRemove:
        os.remove(f)

fileList = find(ext, scriptPath)
for fileToLink in fileList:
    fileName = fileToLink.split('/')[-1]
    os.symlink(fileToLink,
               os.path.join(scriptPath, outDirName, fileName))
