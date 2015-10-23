#!/usr/bin/env python

import os
import glob
import astropy.io.fits as fits
import warnings
import numpy as np
from alipy import pysex
import ConfigParser
from astropy.wcs import WCS

warnings.filterwarnings('ignore')
script_path = os.path.dirname(os.path.realpath(__file__))
work_dir = os.path.curdir


class Config():
    def __init__(self):
        config = ConfigParser.RawConfigParser()
        config.read(os.path.join(script_path, 'config.cfg'))
        self.config = config

        # files
        self.extension = config.get('files', 'extension')
        self.files_to_rm = config.get('files', 'files_to_rm').split(',')
        self.sources_file = config.get('files', 'sources_file')
        self.sources_dict = config.get('files', 'sources_dict')
        # keywords
        self.time_key = config.get('keywords', 'time_key')
        self.gain_key = config.get('keywords', 'gain_key')
        self.exp_key = config.get('keywords', 'exp_key')

        # dirs
        self.output_dir_name = config.get('dirs', 'output_dir')

        # sextractor
        self.phot_type = config.get('sextractor', 'phot_type')
        self.verbose_type = config.get('sextractor', 'verbose_type')
        self.detect_tresh = config.getfloat('sextractor', 'detect_tresh')
        self.analysis_thresh = config.getfloat('sextractor', 'analysis_thresh')
        self.backphoto_type = config.get('sextractor', 'backphoto_type')
        self.backphoto_thick = config.get('sextractor', 'backphoto_thick')
        self.assoc_radius = config.getfloat('sextractor', 'assoc_radius')


class Image:
    def __init__(self, image_name, image_star, image_filter, image_time,
                 image_exp_time, image_hdu, image_header, image_data):
        self.image_name = image_name
        self.image_star = image_star
        self.image_filter = image_filter
        self.image_time = image_time
        self.image_exp_time = image_exp_time
        self.image_hdu = image_hdu
        self.image_header = image_header
        self.image_data = image_data
        self.image_base_name = str(self.image_name).split(".")[0]
        self.image_base_base_name = str(self.image_base_name).split('/')[-1]

    # convert WCS to pix
    def world_2_pix(self):

        w = WCS(self.image_header)
        px_source, py_source = w.wcs_world2pix(float(star.ra_source),
                                               float(star.dec_source), 1)

        return np.array([[px_source, py_source, 1]])

    def flux_measure(self):

        file_to_sex = self.image_name
        coo_file_name = 'coo.coo'
        pix_coo_array = self.world_2_pix()
        # save txt file with pix coordinates for sextractor
        np.savetxt(coo_file_name, pix_coo_array, fmt='%f')

        # run sextractor
        cat = pysex.run(
            file_to_sex, keepcat=False,
            params=['FLUX_BEST', 'FLUXERR_BEST', 'VECTOR_ASSOC(1)'],
            conf_args={'VERBOSE_TYPE': cfg.verbose_type,
                       'BACKPHOTO_TYPE': cfg.backphoto_type,
                       'BACKPHOTO_THICK': cfg.backphoto_thick,
                       'DETECT_THRESH': cfg.detect_tresh,
                       'ANALYSIS_THRESH': cfg.analysis_thresh,
                       'GAIN': float(self.image_header[cfg.gain_key]),
                       'ASSOC_RADIUS': cfg.assoc_radius,
                       'ASSOC_PARAMS': '1,2',
                       'ASSOC_DATA': '3',
                       'ASSOC_NAME': coo_file_name})

        if len(cat) > 0:
            print cat[0][0], cat[1][0]
            print self.image_name, "flux measured"
            print "---------------------------------------------- "
            self.create_image_tab(cat)

        else:
            print self.image_name, "FLUX NOT MEASURED!"
            print "---------------------------------------------- "
            not_measured.append([self.image_name, self.image_star])

    def create_image_tab(self, cat):

        time = self.image_time
        color_filter = str(self.image_filter)
        rotangle = str(self.image_header['ROTANGLE'])
        rotskypa = str(self.image_header['ROTSKYPA'])

        try:
            rotor = str(self.image_header['RROTPOS'])
        except KeyError:
            rotor = 0

        try:
            airmass = str(self.image_header['AIRMASS'])
        except KeyError:
            airmass = 0

        try:
            seeing = str(self.image_header['SEEING'])
            moon_frac = str(self.image_header['MOONFRAC'])
            moon_dist = str(self.image_header['MOONDIST'])
        except KeyError:
            seeing = 0
            moon_frac = 0
            moon_dist = 0

        phase = (time - star.t0) / star.per - int((time - star.t0) / star.per)
        image_output = np.array([[time, cat[0][0], cat[1][0],
                                  color_filter, rotor, phase, airmass,
                                  rotangle, rotskypa, exp_time, seeing,
                                  moon_frac, moon_dist]])
        star.add_point(image_output)


class Star:
    def __init__(self, star_name, ra_source, dec_source, t0, per):
        self.star_name = star_name
        self.ra_source = ra_source
        self.dec_source = dec_source
        self.t0 = t0
        self.per = per
        self.star_images_list = []
        self.star_points = np.zeros(shape=(1, 13))

    def add_point(self, image_output):
        self.star_points = np.append(self.star_points, image_output, axis=0)

    def save_flux_table(self):
        for filter_name in ['D', 'E', 'F']:
            tab = self.star_points[np.where(
                self.star_points[:, 3] == filter_name)]

            c1 = fits.Column(name='TIME', format='D', array=tab[:, 0])
            c2 = fits.Column(name='COUNTS', format='D', array=tab[:, 1])
            c3 = fits.Column(name='COUNTS_ERR', format='D', array=tab[:, 2])
            c4 = fits.Column(name='FILTER', format='1A', array=tab[:, 3])
            c5 = fits.Column(name='ROTOR', format='1A', array=tab[:, 4])
            c6 = fits.Column(name='PHASE', format='D', array=tab[:, 5])
            c7 = fits.Column(name='AIRMASS', format='D', array=tab[:, 6])
            c8 = fits.Column(name='ROTANGLE', format='D', array=tab[:, 7])
            c9 = fits.Column(name='ROTSKYPA', format='D', array=tab[:, 8])
            c10 = fits.Column(name='EXPTIME', format='D', array=tab[:, 9])
            c11 = fits.Column(name='SEEING', format='D', array=tab[:, 10])
            c12 = fits.Column(name='MOON_FRAC', format='D', array=tab[:, 11])
            c13 = fits.Column(name='MOON_DIST', format='D', array=tab[:, 12])

            cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7,
                                 c8, c9, c10, c11, c12, c13])
            tbhdu = fits.new_table(cols)
            prihdr = fits.Header()
            prihdr['OBJECT'] = self.star_name
            prihdr['TIME'] = cfg.time_key
            prihdu = fits.PrimaryHDU(header=prihdr)
            thdulist = fits.HDUList([prihdu, tbhdu])
            thdulist.writeto(output_dir + (self.star_name) + '_' +
                             filter_name + '.fits', clobber=True)

            txt_header = ('Time, source_counts, err, filter, rotor, ' +
                          'phase, airmass,' 'rotangle, rotskypa, exptime, ' +
                          'seeing, moon_frac, moon_dist')
            # write output to CSV table
            np.savetxt(output_dir + (self.star_name) + '_' + filter_name +
                       '.csv', tab, delimiter=',', fmt="%s", header=txt_header)


def cleaning():
    print 'cleaning.....'
    for i in cfg.files_to_rm:
        files = glob.glob(os.path.join(work_dir, i))
        for j in files:
            os.remove(j)


def object_check(images):
    object_list = []
    for im in images:
        hdu = fits.open(im, mode='update')
        hdr = hdu[0].header

        if len(hdu) > 1:
            del(hdu[1:len(hdu)])
            hdu.flush()
        obj = str(hdr['OBJECT'].upper())

        if obj not in object_list:
            object_list.append(obj)

    return object_list


def sources_list_check(sources_list, object_list, sourcesDict):
    miss_object = []
    for obj in object_list:
        print objectInDict(obj, sourcesDict)
        if objectInSource(obj, sources_list):
            pass
        elif objectInDict(obj, sourcesDict) is not None:
            pass
        else:
            miss_object.append(obj)

    if len(miss_object) > 0:
        for miss in miss_object:
            print miss, 'not found in source file or source dict'
        exit()


def createObject(object_name, sources_list):
    coord_tab = sources_list[object_name]
    star = Star(object_name, coord_tab[0],
                coord_tab[1], coord_tab[2], coord_tab[3])
    return star


def objectInSource(object_name, sources_list):
    for source in sources_list:
        if object_name == source:
            return True


def objectInDict(object_name, sourcesDict):
    for source in sourcesDict:
        if object_name in source[1].split(','):
            return source[0]


def objectExist(object_name, star_list):
    if object_name in star_list:
        return True


def loadSourcesList():
    sources_list = {}
    sources_file = np.loadtxt(
        os.path.join(script_path, cfg.sources_file),
        dtype={'names': ('name', 'ra_source', 'dec_source', 't0', 'per'),
               'formats': ('S20', 'f8', 'f8', 'f8', 'f8')})
    for source in sources_file:
        sources_list.update({source[0]: [source[1], source[2],
                                         source[3], source[4]]})
    return sources_list


def loadSourcesDict():
    # read database sources file
    sourcesDict = np.loadtxt(
        os.path.join(script_path, cfg.sources_dict),
        dtype='string', delimiter=':')

    return sourcesDict


cfg = Config()
output_dir = os.path.join(work_dir, cfg.output_dir_name)


try:
    os.mkdir(output_dir)
except OSError:
    pass


star_list = {}
not_measured = []

sources_list = loadSourcesList()
sourcesDict = loadSourcesDict()

# create image list
images = sorted(glob.glob(os.path.join(work_dir, '*'+cfg.extension)))
# create object list
object_list = object_check(images)
# check if each object is in source list
sources_list_check(sources_list, object_list, sourcesDict)


for image in images:
    hdu = fits.open(image, mode='update')
    hdr = hdu[0].header
    data = hdu[0].data
    filter_name = image.split('/')[-1][0].upper()
    object_name = hdr['OBJECT'].upper()

    obs_time = hdr[cfg.time_key]
    exp_time = hdr[cfg.exp_key]

    im = Image(image, object_name, filter_name,
               obs_time, exp_time, hdu, hdr, data)

    if objectExist(object_name, star_list):
        star = star_list[object_name]
        im.flux_measure()
    elif objectInSource(object_name, sources_list):
        star = createObject(object_name, sources_list)
        star_list.update({star.star_name: star})
        im.flux_measure()
    else:
        trueObjectName = objectInDict(object_name, sourcesDict)
        if trueObjectName is not None:
            try:
                star = star_list[trueObjectName]
                im.flux_measure()
            except KeyError:
                star = createObject(trueObjectName, sources_list)
                star_list.update({star.star_name: star})
                im.flux_measure()
        else:
            print 'something wrong :D'
            exit()

# save output for each object
for obj in star_list.iteritems():
    obj[1].save_flux_table()
# remove temp files
cleaning()
