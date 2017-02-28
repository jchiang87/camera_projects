import os
from DataCatalog import DataCatalog
import lsst.eotest.sensor as sensorTest
import eTraveler.clientAPI.connection
import camera_components

def rename_file(infile, i):
    outfile = infile.replace('.fits', '_%02i.fits' % i)
    os.rename(infile, outfile)

user = 'jchiang'
db_name = 'Dev'
prod_server = True
conn = eTraveler.clientAPI.connection.Connection(user, db_name, prod_server)

UNIT_TYPE = 'LCA-11021_RTM'
raft_id = 'LCA-11021_RTM-004_ETU2-Dev'
run_num = '4642D'

folder = os.path.join('/LSST/mirror/BNL-test/test', UNIT_TYPE, raft_id, run_num)
datacat = DataCatalog(folder=folder, site='slac.lca.archive')
def fe55_files(sensor_id):
    query = ' && '.join(('LSST_NUM="%s"' % sensor_id,
                         'TESTTYPE="FE55"', 'IMGTYPE="FE55"'))
    datasets = datacat.find_datasets(query)
    return sorted(datasets.full_paths())

def mask_files(sensor_id):
    return ('/nfs/farm/g/lsst/u1/jobHarness/jh_archive-test/LCA-11021_RTM/LCA-11021_RTM-004_ETU2-Dev/4656D/fe55_raft_analysis/v0/26513/ITL-3800C-145-Dev_rolloff_defects_mask.fits',)

task = sensorTest.Fe55Task()
task.config_fit_xy = False

raft = camera_components.Raft.create_from_connection(conn, raft_id, UNIT_TYPE,
                                                     no_batched='false')

with open('%s_fe55_file_list.txt' % raft_id, 'w') as output:
    for slot, sensor_id in raft.items():
        for i, fe55_file in enumerate(fe55_files(sensor_id)):
            print slot
            output.write("%s  %s\n" % (slot, fe55_file))
#            task.run(sensor_id, (fe55_file,), mask_files(sensor_id))
#            rename_file('%s_eotest_results.fits' % sensor_id, i)
#            rename_file('%s_psf_results_nsig4.fits' % sensor_id, i)
