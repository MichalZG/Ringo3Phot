#!/usr/bin/env python

import os
import glob


temp = glob.glob('[def]_e_*_1.fits')
for fil in temp:
    os.remove(fil)
