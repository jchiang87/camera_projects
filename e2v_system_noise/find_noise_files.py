import os
from collections import namedtuple, OrderedDict
import subprocess
import numpy as np
import siteUtils
import lsst.eotest.sensor as sensorTest

def make_vendor_file_list(outfile='file_list.txt'):
    command = '(find /nfs/farm/g/lsst/u1/vendorData/E2V -name \*_Noise_Multiple_Samples_Summary\*.csv) >& %(outfile)s' % locals()
    subprocess.call(command, shell=True)

AmpNoise = namedtuple('AmpNoise', ['read_noise_adu', 'read_noise_e',
                                   'total_noise_adu', 'total_noise_e'])

class SensorNoise(dict):
    def __init__(self, infile):
        super(SensorNoise, self).__init__()
        with open(infile) as input_:
            for line in input_:
                if line.startswith('Amp'):
                    continue
                tokens = line.strip().split(',')
                self[int(tokens[0])] = AmpNoise(*(float(x) for x in tokens[1:]))

def find_noise_data(list_of_files):
    e2v_data = OrderedDict()
    with open(list_of_files) as file_list:
        for line in file_list:
            if not line.startswith('/nfs'):
                continue
            sensor_id = line.split('/')[8]
            job_id = line.split('/')[10]
            e2v_data[sensor_id] = SensorNoise(line.strip())

    return e2v_data

def get_eotest_results(sensor_id, data_product='EOTEST_RESULTS'):
    query = ' && '.join(('LSST_NUM=="%(sensor_id)s"',
                         'DATA_PRODUCT=="%(data_product)s"')) % locals()
    datasets = siteUtils.datacatalog_query(query)
    return sorted([item for item in datasets.full_paths()])


if __name__ == '__main__':
    vendor_list = 'file_list.txt'
    #make_vendor_file_list(vendor_list)

    os.environ['LCATR_DATACATALOG_FOLDER'] = \
        '/LSST/mirror/SLAC-prod/prod/e2v-CCD'
    e2v_data = find_noise_data('file_list.txt')

    for item in e2v_data.keys():
        sensor_id = item
        eotest_results_file = get_eotest_results(sensor_id)[-1]
#        print sensor_id, eotest_results_file

        eotest_results = sensorTest.EOTestResults(eotest_results_file)

        print sensor_id
        lsst_total_noise = eotest_results['TOTAL_NOISE']/eotest_results['GAIN']
        e2v_system_noise = [np.sqrt(e2v_data[sensor_id][i].total_noise_adu**2
                                    - e2v_data[sensor_id][i].read_noise_adu**2)
                            for i in range(1, 17)]
        print "amp  total_noise     system_noise"
        print "     (ADU, eotest)   (ADU, e2v)"
        for i in range(16):
            amp = i + 1
            print '%2i       %.2f           %.2f' % (amp, lsst_total_noise[i],
                                                     e2v_system_noise[i])
        print
