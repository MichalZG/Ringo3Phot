#! /usr/bin/env python
import glob
import os
import fnmatch


ext = '*.fits'
outDirName = 'results'
workPath = os.getcwd()


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
    filesToRemove = glob.glob(os.path.join(workPath, outDirName, ext))
    for f in filesToRemove:
        os.remove(f)

fileList = find(ext, workPath)
for fileToLink in fileList:
    fileName = fileToLink.split('/')[-1]
    print fileToLink
    os.symlink(fileToLink,
               os.path.join(workPath, outDirName, fileName))
