import os
import glob
from lsst.eotest.fitsTools import fitsWriteto
from lsst.eotest.sensor.ctesim import ctesim_cpp as ctesim

infiles = sorted(glob.glob('sim_data/000-01_fe55_fe55*.fits'))

nfiles = 5

pcti = 0
scti = 1e-3
for infile in infiles[:nfiles]:
    foo = ctesim(infile, pcti=pcti, scti=scti, verbose=True)
    tokens = os.path.basename(infile).split('_')
    outfile = '_'.join(tuple(tokens[:-1] + ['s%(scti)4.0e.fits' % locals()]))
    print outfile
    fitsWriteto(foo, outfile, clobber=True)
