#!/usr/bin/env python3

from glob import glob
from os import path

if __name__ == '__main__':
    this_dir = path.dirname(path.abspath(__file__))
    remove = len(this_dir) + 1
    sources = [
        fname[remove:]
        for dirname in ['Core', 'IO', 'Background', 'PSF', 'FitPSF','SubPixPhot']
        for fname in glob(path.join(this_dir, 'src', dirname, '*.cpp'))
    ]
    for exclude in ['src/FitPSF/SDKSourceBase.cpp',
                    'src/FitPSF/SDKUtil.cpp']:
        sources.remove(exclude)

    print(sources)
