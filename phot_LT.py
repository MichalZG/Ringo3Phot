import os
import shutil
import glob
import astropy.io.fits as fits
import astropy.wcs as wcs
import warnings
import numpy as np
import astropy.wcs
import subprocess as sub
from alipy import pysex
import time
from astropy.wcs import WCS
warnings.filterwarnings('ignore')


#extension
extension = '.fits'

#frame keywords
time_key = 'MJD'
gain_key = 'GAIN'
exp_key = 'EXPTIME'

#dirs
work_dir = '/home/pi/kolchoz/LT/exo2013/'
output_dir = work_dir+'output/'

#sextractor param
detect_tresh = 2.0
analysis_thresh = 2.0
backphoto_thick = 64 #thickness of annulus

#align
n = 500 # Number of brightest stars of each image to consider (default 500)
visu = False
makepng = False
skipsaturated = False

script_path = os.path.dirname(os.path.realpath(__file__))
files_to_rm = ['*.coo']

class Image:
  def __init__(self, image_name, image_star, image_filter, image_time, image_exp_time, image_hdu, image_header, image_data):
    self.image_name = image_name
    self.image_star = image_star
    self.image_filter = image_filter
    self.image_time = image_time
    self.image_exp_time = image_exp_time
    self.image_hdu = image_hdu
    self.image_header = image_header
    self.image_data = image_data

    self.image_base_name = str(self.image_name).split(".")[0] # with dir, without extension
    self.image_base_base_name = str(self.image_base_name).split('/')[-1] # only file name without dir and extension

  def world_2_pix(self): #convert WCS to pix

    w = WCS(self.image_header)
    px_source, py_source = w.wcs_world2pix(float(star.ra_source), float(star.dec_source), 1)

    return np.array([[px_source, py_source, 1]])

  def flux_measure(self):

    hdu_flux = fits.open(work_dir+self.image_base_base_name+extension)
    hdr_flux = hdu_flux[0].header
    pix_coo_array = self.world_2_pix()
    np.savetxt(work_dir+self.image_base_base_name+'.coo', pix_coo_array, fmt='%f') #save txt file with pix coordinates for sextractor

    cat = pysex.run(work_dir+self.image_base_base_name+extension, keepcat=False, #run sextractor
		     params=['FLUX_BEST', 'FLUXERR_BEST', 'VECTOR_ASSOC(1)'], #values list in sextractor output
		     conf_args={'VERBOSE_TYPE':'QUIET',
				'BACKPHOTO_TYPE':'LOCAL',
	                        'BACKPHOTO_THICK':backphoto_thick,
		                'DETECT_THRESH':detect_tresh,
	                        'ANALYSIS_THRESH':analysis_thresh,
	                        'GAIN':float(self.image_header[gain_key]),
	                        'ASSOC_RADIUS':4.0,
	                        'ASSOC_PARAMS':'1,2',
	                        'ASSOC_DATA':'3',
                                'ASSOC_NAME':work_dir+self.image_base_base_name+'.coo'})

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
    exptime = str(self.image_header['EXPTIME'])

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
    image_output = np.array([[time, cat[0][0], cat[1][0],\
			                 color_filter, rotor, phase, airmass, rotangle,\
		                     rotskypa, exp_time, seeing, moon_frac, moon_dist]])
    star.add_point(image_output)


class Star:
  def __init__(self, star_name, ra_source, dec_source, t0, per):
    self.star_name = star_name
    self.ra_source = ra_source
    self.dec_source = dec_source
    self.t0 = t0
    self.per = per
    self.star_images_list = []
    self.star_points = np.zeros(shape=(1,13))

  def add_point(self, image_output):
    self.star_points = np.append(self.star_points, image_output, axis=0)

  def save_flux_table(self):
    tab = np.delete(self.star_points, (0), axis=0)

    c1 = fits.Column(name='TIME', format='D', array=tab[:,0])
    c2 = fits.Column(name='COUNTS', format='D', array=tab[:,1])
    c3 = fits.Column(name='COUNTS_ERR', format='D', array=tab[:,2])
    c4 = fits.Column(name='FILTER', format='1A', array=tab[:,3])
    c5 = fits.Column(name='ROTOR', format='1A', array=tab[:,4])
    c6 = fits.Column(name='PHASE', format='D', array=tab[:,5])
    c7 = fits.Column(name='AIRMASS', format='D', array=tab[:,6])
    c8 = fits.Column(name='ROTANGLE', format='D', array=tab[:,7])
    c9 = fits.Column(name='ROTSKYPA', format='D', array=tab[:,8])
    c10 = fits.Column(name='EXPTIME', format='D', array=tab[:,9])
    c11 = fits.Column(name='SEEING', format='D', array=tab[:,10])
    c12 = fits.Column(name='MOON_FRAC', format='D', array=tab[:,11])
    c13 = fits.Column(name='MOON_DIST', format='D', array=tab[:,12])

    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13])
    tbhdu = fits.new_table(cols)
    prihdr = fits.Header()
    prihdr['OBJECT'] = self.star_name
    prihdr['TIME'] = time_key
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(output_dir+(self.star_name)+'.fits', clobber=True) #write output to fits table

    txt_header = 'Time, source_counts, err, filter, rotor, phase, airmass, rotangle, rotskypa, exptime, seeing, moon_frac, moon_dist'
    np.savetxt(output_dir+(self.star_name)+'.csv', tab, delimiter=',', fmt="%s", header=txt_header) # write output to CSV table

def cleaning():
  print 'cleaning.....'
  for i in files_to_rm:
    files = glob.glob(work_dir+i)
    for j in files:
      os.remove(j)

def object_check(images):
  object_list = []
  for im in images:
    filter_name = im.split('/')[-1][0].upper()
    hdu = fits.open(im, mode='update')
    hdr = hdu[0].header

    if len(hdu) > 1:
      del(hdu[1:len(hdu)])
      hdu.flush()
    obj = str(hdr['OBJECT'].upper()+'_'+filter_name)

    if obj not in object_list:
      object_list.append(obj)

  return object_list

def sources_list_check(sources_file, object_list):
  temp = []
  miss_object = []
  for source in sources_file:
    temp.append(source[0])
    print source[0]
  for obj in object_list:
    if obj in temp:
      pass
    else:
      if obj not in miss_object:
	miss_object.append(obj)

  if len(miss_object) > 0:
    for miss in miss_object:
      print miss, 'not found in source file'
    exit()

def save_not_measured_list(not_measured):
  pass
  #np.savetxt(output_dir+'not_measured.dat', not_measured, delimiter=' ', fmt="%s") # write name of not measured images to txt


try:
  os.mkdir(output_dir)
except OSError:
  pass

star_list = {}
sources_list = {}
not_measured = []

#read database sources file
sources_file = np.loadtxt(script_path+'sources.dat',\
			 dtype={'names': ('name', 'ra_source', 'dec_source', 't0', 'per'),\
	                 'formats': ('S20', 'f8', 'f8', 'f8', 'f8')})


images = sorted(glob.glob(work_dir+"*"+extension)) #create image list
object_list = object_check(images) #create object list
sources_list_check(sources_file, object_list) #check if each object is in source list

for source in sources_file:
  sources_list.update({source[0]:[source[1], source[2], source[3], source[4]]})

for image in images:
  hdu = fits.open(image, mode='update')
  hdr = hdu[0].header
  data = hdu[0].data
  filter_name = image.split('/')[-1][0].upper()
  object_name = hdr['OBJECT'].upper()+'_'+filter_name

  obs_time = hdr[time_key]
  exp_time = hdr[exp_key]

  im = Image(image, object_name, filter_name, obs_time, exp_time, hdu, hdr, data)

  if object_name in star_list:
    star = star_list[object_name]
  else:
    coord_tab = sources_list[object_name]
    star = Star(object_name, coord_tab[0], coord_tab[1], coord_tab[2], coord_tab[3])
    star_list.update({object_name:star})

  im.flux_measure()


for obj in star_list.iteritems(): #save output for each object
  obj[1].save_flux_table()

cleaning() #remove temp files
save_not_measured_list(not_measured)
