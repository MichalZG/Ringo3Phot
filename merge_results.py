#!/usr/bin/env python

import astropy.io.fits as fits
import fnmatch
import os
import shutil

from astropy.table import Table, join

sources = ("V", "4U", "SAX", "RX", "KS", "EXO", "LSI", "GRO")

# scan all output directories for fits files

files = []
for root, dirnames, filenames in os.walk("."):
    current_dir = os.path.basename(root)
    if current_dir == "output":
        print "+++++++++++++++++++++++++"
        print root + ":"
        print "+++++++++++++++++++++++++"
        for source in sources:
            print "--------------"
            print source
            print "--------------"
            filt = source + "*.fits"
            for filename in fnmatch.filter(filenames, filt):
                print filename
                full_path = os.path.join(root, filename)
                colour = filename.split(".")[-2].split("_")[-1]
                merged_filename = source + "_" + colour + ".fits"
                results_file = os.path.join("results", merged_filename)
                if os.path.isfile(results_file):
                    #print "exists"
                    source_fits = fits.open(full_path)
                    results_fits = fits.open(results_file)
                    source_table = source_fits[1].data
                    results_table = results_fits[1].data
                    source_unique_time = list(set([str(v) for v in source_table.field("TIME")]))[0]
                    #print "src:", source_unique_time
                    results_unique_times = list(set([str(v) for v in results_table.field("TIME")]))
                    #print "res:", results_unique_times
                    if source_unique_time not in results_unique_times:
                        t1 = source_fits
                        t2 = results_fits
                        nrows1 = t1[1].data.shape[0]
                        nrows2 = t2[1].data.shape[0]

                        nrows = nrows1 + nrows2
                        hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
                        for colname in t1[1].columns.names:
                            hdu.data[colname][nrows1:] = t2[1].data[colname]
                        hdu.writeto(results_file, clobber='True')
                        print "merged..."
                    else:
                        print "already merged..."
                else:
                    print "creating..."
                    shutil.copyfile(full_path, results_file)
