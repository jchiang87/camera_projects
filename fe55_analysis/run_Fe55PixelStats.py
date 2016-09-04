from __future__ import absolute_import, print_function
import os
import glob
import logging
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest

plt.ion()

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

infiles = sorted(glob.glob('000-01_fe55_fe55_*_scti.fits'))
pickle_file = '000-01_fe55_stats.pkl'

if not os.path.isfile(pickle_file):
    fe55_pixel_stats = sensorTest.Fe55PixelStats(infiles, logger=logger)
    fe55_pixel_stats.to_pickle(pickle_file)
else:
    fe55_pixel_stats = sensorTest.Fe55PixelStats.read_pickle(pickle_file)
#    fe55_pixel_stats.set_selection_function('kalpha')
    fe55_pixel_stats.set_selection_function('amp')
    fe55_pixel_stats.sensor_id = '000-01'

fe55_pixel_stats.apflux_profile()

fe55_pixel_stats.dn_hists()

fe55_pixel_stats.pixel_hists('p3', 'p5')

figure, sresults = fe55_pixel_stats.pixel_diff_profile('x', 'p3', 'p5')
print(sresults)
