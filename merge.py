#/usr/bin/env python

import glob
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table, join

new_tables = glob.glob("*.fits")
old_tables = glob.glob("/media/disk1/data/bexb/results/*.fits")

for new_tab in new_tables:
    for old_tab in old_tables:
        if new_tab == old_tab.split("/")[-1]:
            t1 = fits.open(new_tab)
            t2 = fits.open(old_tab)
            nrows1 = t1[1].data.shape[0]
            nrows2 = t2[1].data.shape[0]

            nrows = nrows1 + nrows2
            hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
            for colname in t1[1].columns.names:
                hdu.data[colname][nrows1:] = t2[1].data[colname]
            hdu.writeto(old_tab, clobber='True')
