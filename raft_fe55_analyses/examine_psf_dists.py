import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
from fe55_files import fits_files
plt.ion()

chiprob_min = 0.1

slot = 'S10'
psf_catalog = os.path.basename(os.path.basename(fits_files[slot])[:len('ITL-3800C-111-Dev')]
                               + '_psf_results_nsig4_00.fits')
fit_range = (0.05, 0.95)
nsig = 5
for amp in range(1, 17):
    try:
        data = fits.open(psf_catalog)[amp].data
        index = np.where((data['CHIPROB'] > chiprob_min) &
                         (data['DN'] < 5000))
        dn = sorted(data['DN'][index])
        npts = len(dn)
        median = np.median(dn)
        stdev = np.std(dn[int(npts*fit_range[0]):int(npts*fit_range[1])])
        dn_range = median - nsig*stdev, median + nsig*stdev
        plt.figure()
        plt.hist(dn, bins=50, histtype='step', range=dn_range)
        plt.title('%s, %s' % (slot, imutils.channelIds[amp]))
        print imutils.channelIds[amp], min(dn), max(dn), \
            min(data['DN']), max(data['DN'])
    except StandardError as eobj:
        print amp, eobj
