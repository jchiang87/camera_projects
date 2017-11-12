import glob
import numpy as np
import numpy.random as random
import astropy.io.fits as fits
import lsst.afw.detection as afw_detect
import lsst.afw.image as afw_image
import lsst.afw.math as afw_math
import lsst.eotest.sensor as sensorTest

def bg_image(ccd, amp, nx=10, ny=10):
    bg_ctrl = afw_math.BackgroundControl(nx, ny, ccd.stat_ctrl)
    bg = afw_math.makeBackground(ccd[amp], bg_ctrl)
    return bg.getImageF()

class IsMasked(object):
    def __init__(self, mask_file):
        self.ccd = sensorTest.MaskedCCD(mask_file)
    def __call__(self, amp, footprint):
        mask_imarr = self.ccd[amp].getImage().getArray()
        for span in footprint.getSpans():
            iy = span.getY()
            for ix in range(span.getX0(), span.getX1()+1):
                if mask_imarr[iy][ix] > 0:
                    return True
        return False

def make_mask(med_file, gains, outfile, ethresh=0.1, colthresh=20,
              mask_plane='BAD'):
    ccd = sensorTest.MaskedCCD(med_file)
    exptime = ccd.md.get('EXPTIME')
    pixels, columns = {}, {}
    for amp in ccd:
        bright_pixels \
            = sensorTest.BrightPixels(ccd, amp, exptime, gains[amp],
                                      ethresh=ethresh, colthresh=colthresh)
        pixels[amp], columns[amp] = bright_pixels.find()
    sensorTest.generate_mask(med_file, outfile, mask_plane, pixels=pixels,
                             columns=columns)
    return IsMasked(outfile)

sensor_id = 'ITL-3800C-244'

fe55_results = sensorTest.EOTestResults('ITL-3800C-244_eotest_results.fits')
gains = {amp: gain for amp, gain in
         zip(fe55_results['AMP'], fe55_results['GAIN'])}

mask_file = '%(sensor_id)s_bp_mask_file.fits' % locals()
med_file = 'ITL-3800C-244_median_dark_bp.fits'
#mask = make_mask(med_file, gains, mask_file)
is_masked = IsMasked(mask_file)

darks = glob.glob('%(sensor_id)s_dark_dark_*.fits' % locals())

nsig = 7

fp_id, x0, y0, pixel_values = [], [], [], []
index = -1
for dark in darks:
    print "processing", dark
    ccd = sensorTest.MaskedCCD(dark, mask_files=(mask_file,))
    for amp in ccd.keys():
        image = ccd[amp]
        image -= bg_image(ccd, amp)
        stats = afw_math.makeStatistics(image,
                                        afw_math.MEDIAN | afw_math.STDEVCLIP,
                                        ccd.stat_ctrl)
        threshold = \
            afw_detect.Threshold(stats.getValue(afw_math.MEDIAN) +
                                 nsig*stats.getValue(afw_math.STDEVCLIP))
        fp_set = afw_detect.FootprintSet(image, threshold)
        fps = [fp for fp in fp_set.getFootprints()]

        for fp in fps:
            if is_masked(amp, fp):
                continue
            index += 1
            for span in fp.getSpans():
                fp_id.append(index)
                iy = span.getY()
                ix0, ix1 = span.getX0(), span.getX1()
                x0.append(ix0)
                y0.append(iy)
                pixel_values.append(
                    np.array(image.getImage().getArray()[iy, ix0:ix1+1],
                             dtype=np.int)
                    )

hdu_list = fits.HDUList([fits.PrimaryHDU()])
columns = [fits.Column(name='fp_id', format='I', array=fp_id),
           fits.Column(name='x0', format='I', array=x0),
           fits.Column(name='y0', format='I', array=y0),
           fits.Column(name='pixel_values', format='PJ()',
                       array=np.array(pixel_values, dtype=np.object))]
hdu_list.append(fits.BinTableHDU.from_columns(columns))
hdu_list[-1].name = 'COSMIC_RAYS'
hdu_list.writeto('cosmic_ray_catalog.fits', overwrite=True)
