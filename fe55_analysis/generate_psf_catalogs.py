"""
Use the eotest code to generate the psf catalogs for the Fe55 data.
"""
from __future__ import absolute_import
import os
import glob
import logging
import lsst.eotest.sensor as sensorTest

logging.basicConfig()
log = logging.getLogger()

def generate_psf_catalogs(infiles, nfiles=None):
    "Generate the Fe55 catalogs, if necessary."
    for infile in infiles[:nfiles]:
        ext = '_'.join(os.path.basename(infile).split('_')[3:])
        psf_results_file = 'psf_results_%s' % ext
        if not os.path.isfile(psf_results_file):
            log.info("processing %s", os.path.basename(infile))
            fitter = sensorTest.PsfGaussFit(nsig=4, fit_xy=False)
            ccd = sensorTest.MaskedCCD(infile)
            for amp in ccd:
                log.info('Working on %i', amp)
                fitter.process_image(ccd, amp)
            fitter.write_results(psf_results_file)


if __name__ == '__main__':
#    infiles_ = sorted(glob.glob('data/E2V-CCD250-088_fe55_fe55_*.fits'))
#    infiles_ = sorted(glob.glob('000-01_fe55_fe55_??_scti.fits'))
    infiles_ = sorted(glob.glob('/nfs/farm/g/lsst/u1/jobHarness/jh_archive/ITL-CCD/ITL-3800C-089/vendorIngest/v0/16310/fe55/20160412145143/ITL-3800C-089_fe55_fe55_001_20160412145143.fits'))
    nfiles_ = len(infiles_)
    dn_range = (300, 370)
    generate_psf_catalogs(infiles_, nfiles=nfiles_)
