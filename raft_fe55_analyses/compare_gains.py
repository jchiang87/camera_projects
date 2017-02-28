import os
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
from fe55_files import sensors

ts8_results_path = '/nfs/farm/g/lsst/u1/mirror/BNL-test/test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4642D/collect_raft_results/v0/26554'

fe55_gains = dict()
ptc_gains = dict()
fig = plt.figure()
for slot, sensor_id in sensors.items():
    results = sensorTest.EOTestResults('%s_eotest_results_00.fits' % sensor_id)
    fe55_gains[slot] = dict([pair for pair in
                             zip(results['AMP'], results['GAIN'])])
    ts8_resfile = os.path.join(ts8_results_path,
                               sensor_id + '_eotest_results.fits')
    if not os.path.isfile(ts8_resfile):
        raise RuntimeError('file not found:', ts8_resfile)
    ts8_results = sensorTest.EOTestResults(ts8_resfile)
    ptc_gains[slot] = dict([pair for pair in zip(ts8_results['AMP'],
                                                 ts8_results['PTC_GAIN'])])

    plt.errorbar(results['GAIN'], ts8_results['PTC_GAIN'],
                 xerr=results['GAIN_ERROR'], yerr=ts8_results['PTC_GAIN_ERROR'],
                 fmt='.', label=slot)
xmin, xmax, ymin, ymax = 0, 2, 0, 2
plt.xlabel('Fe55 Gain (e-/ADU)')
plt.ylabel('PTC Gain (e-/ADU)')
plt.plot([xmin, xmax], [ymin, ymax], 'k:')
plt.axis((xmin, xmax, ymin, ymax))
plt.legend()
