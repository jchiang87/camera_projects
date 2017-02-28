from __future__ import print_function
from collections import OrderedDict
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
from lsst.eotest.sensor.MaskedCCD import MaskedCCDBiasImageException
from fe55_files import fits_files

def display_segment(ccd, amp, nx=10, ny=10):
    try:
        image = ccd.bias_subtracted_image(amp)
    except MaskedCCDBiasImageException:
        image = ccd[amp]
    bg_ctrl = afwMath.BackgroundControl(nx, ny, ccd.stat_ctrl)
    bg = afwMath.makeBackground(ccd[amp], bg_ctrl)
    image -= bg.getImageF()
    ds9.mtv(image)

# "Good" example
slot = 'S00'
amp = 1

# zero clusters
slot = 'S00'
amp = 4

# zero gain, zero error (high end DN tail ~1e9)
slot = 'S00'
amp = 16

# tiny gain, inf error, lots of clusters
slot = 'S01'
amp = 4

ccd = dict()
ccd['S00'] = sensorTest.MaskedCCD(fits_files['S00'])
ccd['S01'] = sensorTest.MaskedCCD(fits_files['S01'])

#display_segment(ccd['S00'], 1)
#display_segment(ccd['S00'], 4)
#display_segment(ccd['S00'], 16)
display_segment(ccd['S01'], 4)

#for amp in (1, 4, 16):
#    display_segment(ccd['S00'], amp)
#    ds9.incrDefaultFrame()
#
#display_segment(ccd['S01'], 4)
