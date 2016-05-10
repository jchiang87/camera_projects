import glob
import matplotlib.pyplot as plt
import logging
import lsst.eotest.sensor as sensorTest

plt.ion()

logger = logging.getLogger()

infiles = sorted(glob.glob('000-01_fe55_fe55_*_scti.fits'))

fe55_pixel_stats = sensorTest.Fe55PixelStats(infiles, logger=logger)

fe55_pixel_stats.pixel_hists('p3', 'p5')
figure, sresults = fe55_pixel_stats.pixel_diff_profile('x', 'p3', 'p5')
print(sresults)

fe55_pixel_stats.pixel_hists('p1', 'p7')
figure, presults = fe55_pixel_stats.pixel_diff_profile('y', 'p1', 'p7')
print(presults)
