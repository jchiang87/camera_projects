import lsst.afw.detection as afw_detect
import lsst.afw.image as afw_image
import lsst.afw.math as afw_math

import lsst.eotest.sensor as sensorTest

def bg_image(ccd, amp, nx=10, ny=10):
    bg_ctrl = afw_math.BackgroundControl(nx, ny, ccd.stat_ctrl)
    bg = afw_math.makeBackground(ccd[amp], bg_ctrl)
    return bg.getImageF()

dark = 'ITL-3800C-244_dark_dark_001_20171004104341.fits'
ccd = sensorTest.MaskedCCD(dark)


nsig = 4
for amp in ccd.keys()[:1]:
    image = ccd[amp]
    image -= bg_image(ccd, amp)
    stats = afw_math.makeStatistics(image, afw_math.MEDIAN | afw_math.STDEVCLIP,
                                    ccd.stat_ctrl)
    threshold = afw_detect.Threshold(stats.getValue(afw_math.MEDIAN) +
                                     nsig*stats.getValue(afw_math.STDEVCLIP))
    fp_set = afw_detect.FootprintSet(image, threshold)
    fps = [fp for fp in fp_set.getFootprints() if fp.getSpans()[0].getY() > 50]


outfile = 'my_footprints.fits'
with open(outfile, 'w') as output:
    fps[0].writeFits(outfile)
    for fp in fps[1:]:
        fp.writeFits(outfile, 'a')
