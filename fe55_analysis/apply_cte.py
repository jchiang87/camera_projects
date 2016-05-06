"""
Apply simulated serial cte to simulated Fe55 data.
"""
from __future__ import absolute_import, print_function
import os
import glob
import numpy as np
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.eotest.sensor as sensorTest
import lsst.eotest.image_utils as imutils
import lsst.eotest.utilLib as eotest_utils
from lsst.eotest.fitsTools import fitsWriteto
from lsst.eotest.sensor.ctesim import fitsFile
from cti_amp_dict import pcti_, scti_

def apply_cte(infile, pcti=None, scti=None, verbose=False,
              amps=tuple(range(1, 17))):
    """
    Function to apply distinct levels of parallel cte and/or serial cte
    on each amplifier in an input sensor image.
    """
    if pcti is None:
        pcti = dict([(amp, 0) for amp in amps])
    if scti is None:
        scti = dict([(amp, 0) for amp in amps])
    segments = {}
    for amp in amps:
        if verbose:
            print("apply_cte: working on amp", amp)
        image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
        geom = sensorTest.makeAmplifierGeometry(infile)
        outimage = eotest_utils.ImageTools.applyCTI(image, geom.serial_overscan,
                                                    pcti[amp], scti[amp],
                                                    verbose)
        segments[amp] = outimage
    return fitsFile(segments, fits.open(infile))

infiles = sorted(glob.glob('sim_data/000-01_fe55_fe55*.fits'))
nfmin = 0
nfmax = min(len(infiles), 25)
for infile_ in infiles[nfmin:nfmax]:
    foo = apply_cte(infile_, pcti=pcti_, scti=scti_, verbose=True)
    tokens = os.path.basename(infile_).split('_')
    outfile = '_'.join(tuple(tokens[:-1] + ['scti.fits']))
    print(outfile)
    fitsWriteto(foo, outfile, clobber=True)
