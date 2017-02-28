from __future__ import print_function
import glob
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raftTest
from lsst.eotest.sensor.MaskedCCD import MaskedCCDBiasImageException
import eTraveler.clientAPI.connection
import camera_components

def compute_fe55_thresholds(fe55_file, nsig=4, nx=10, ny=10):
    ccd = sensorTest.MaskedCCD(fe55_file)
    thresholds = dict()
    for amp in ccd:
        try:
            image = ccd.bias_subtracted_image(amp)
        except MaskedCCDBiasImageException:
            image = ccd[amp]
        bg_ctrl = afwMath.BackgroundControl(nx, ny, ccd.stat_ctrl)
        bg = afwMath.makeBackground(ccd[amp], bg_ctrl)
        image -= bg.getImageF()
        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)
        thresholds[amp] = median + nsig*stdev
    return thresholds

class Fe55_segment_processor(object):
    def __init__(self, fe55_files, results_files, psf_catalogs,
                 min_chiprob=0.1):
        self.results = OrderedDict()
        self.gains = OrderedDict()
        self.gain_errors = OrderedDict()
        self.thresholds = OrderedDict()
        for slot, results_file in results_files.items():
            self.results[slot] = sensorTest.EOTestResults(results_file)
            self.gains[slot] = dict(pair for pair in
                                    zip(self.results[slot]['AMP'],
                                        self.results[slot]['GAIN']))
            self.gain_errors[slot] = dict(pair for pair in
                                          zip(self.results[slot]['AMP'],
                                              self.results[slot]['GAIN_ERROR']))
            self.thresholds[slot] = compute_fe55_thresholds(fe55_files[slot])
        self._read_psf_catalogs(psf_catalogs)
        self.min_chiprob = min_chiprob
        self.use_threshold = False

    def _read_psf_catalogs(self, psf_catalogs):
        self.chiprobs = OrderedDict()
        for slot, catalog_file in psf_catalogs.items():
            chiprob_entries = OrderedDict()
            catalog = fits.open(catalog_file)
            for i, hdu in enumerate(catalog[1:17]):
                chiprob_entries[i+1] = hdu.data['CHIPROB']
            self.chiprobs[slot] = chiprob_entries

    def __call__(self, slot, ccd, amp, xy_bounds=None):
        if self.use_threshold:
            return self.thresholds[slot][amp]
        return len(np.where(self.chiprobs[slot][amp] > self.min_chiprob)[0])

    def write_gain_stats(self):
        for slot in self.gains:
            print(slot, self.results[slot].infile[:len('ITL-3800C-012-Dev')])
            print("segment   gain    gain_error    #det_clusters   #selected    fp_threshold")
            for amp in self.gains[slot]:
                print("%-2s   %9.3f   %11.2e   %10i   %10i   %13.2f" %
                      (imutils.channelIds[amp],
                       self.gains[slot][amp],
                       self.gain_errors[slot][amp],
                       len(self.chiprobs[slot][amp]),
                       self(slot, None, amp),
                       self.thresholds[slot][amp]))
            print("")

user = 'jchiang'
db_name = 'Dev'
prod_server = True
conn = eTraveler.clientAPI.connection.Connection(user, db_name, prod_server)
raft_id = 'LCA-11021_RTM-004_ETU2-Dev'
UNIT_TYPE = 'LCA-11021_RTM'
raft = camera_components.Raft.create_from_connection(conn, raft_id, UNIT_TYPE,
                                                     no_batched='false')

fits_files = OrderedDict()
with open('LCA-11021_RTM-004_ETU2-Dev_fe55_file_list.txt') as input:
    for line in input:
        slot, fits_file = line.strip().split()
        fits_files[slot] = fits_file

results_files = OrderedDict()
psf_catalogs = OrderedDict()
for slot, sensor_id in raft.items():
    results_files[slot] = glob.glob('%s_eotest_results*.fits' % sensor_id)[0]
    psf_catalogs[slot] = glob.glob('%s_psf_results_nsig4*.fits' % sensor_id)[0]

seg_processor = Fe55_segment_processor(fits_files, results_files, psf_catalogs)
seg_processor.write_gain_stats()

#fe55_mosaic = raftTest.RaftMosaic(fits_files, gains=seg_processor.gains,
#                                  segment_processor=seg_processor)
#fe55_mosaic.plot(title='%s, Fe55 stats' % raft_id)
#plt.savefig('%s_fe55_stats.png' % raft_id)
#
#seg_processor.use_threshold = True
#fe55_threshold_mosaic = raftTest.RaftMosaic(fits_files,
#                                            gains=seg_processor.gains,
#                                            segment_processor=seg_processor)
#fe55_threshold_mosaic.plot(title='%s, Fe55 cluster threshold' % raft_id)
#plt.savefig('%s_fe55_threshold.png' % raft_id)
