import os
import glob
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import lsst.eotest.sensor as sensorTest
import lsst.eotest.image_utils as imutils
from profile_plot import profile_plot

plt.ion()

logging.basicConfig()
log = logging.getLogger()

def xy_data(infile, amp, dn_range):
    psf_results = fits.open(infile)
    data = psf_results[amp].data
    selection = (dn_range[0] < data['DN']) & (data['DN'] < dn_range[1])
    return data['XPOS'][selection], data['YPOS'][selection], \
        data['CHIPROB'][selection]

def find_peak_pixel(ix, iy, imarr):
    try:
        max_val = max(imarr[iy-1:iy+2, ix-1:ix+2].flatten())
    except ValueError:
        return ix, iy
    index = np.where(imarr[iy-1:iy+2, ix-1:ix+2] == max_val)
    return index[1][0]+ix-1, index[0][0]+iy-1

def generate_psf_catalogs(infiles, nfiles=None):
    "Generate the Fe55 catalogs, if necessary."
    for i, infile in enumerate(infiles[:nfiles]):
        ext = '_'.join(os.path.basename(infile).split('_')[3:])
        psf_results_file = 'psf_results_%s' % ext
        if not os.path.isfile(psf_results_file):
            log.info("processing %s" % os.path.basename(infile))
            fitter = sensorTest.PsfGaussFit(nsig=4, fit_xy=False)
            ccd = sensorTest.MaskedCCD(infile)
            for amp in ccd:
                log.info('Working on %i' % amp)
                fitter.process_image(ccd, amp)
            fitter.write_results(psf_results_file)

def create_data_frame(infiles, nfiles=None, dn_range=None):
    "Create the pixel value data frame."
    log.info(str(dn_range))
    rows = []
    for i, infile in enumerate(infiles[:nfiles]):
        log.info("processing file: %s" % infile)
        ccd = sensorTest.MaskedCCD(infile)
        ext = '_'.join(os.path.basename(infile).split('_')[3:])
        psf_results_file = 'psf_results_%s' % ext
        for amp in range(1, 17):
            xpos, ypos, chiprob = xy_data(psf_results_file, amp, dn_range)
            mi = ccd.bias_subtracted_image(amp)
            imarr = mi.getImage().getArray()
            log.info("for amp %i, number of clusters: %i" % (amp, len(xpos)))
            for x, y, cp in zip(xpos, ypos, chiprob):
                ix, iy = int(np.round(x)), int(np.round(y))
                # Pixel layout.  p4 contains the Fe55 peak. p0 is closest
                # to the output node.
                #   6 7 8
                #   3 4 5
                #<  0 1 2
                try:
                    log.debug(imarr[iy-1:iy+2, ix-1:ix+2])
                    my_row = (amp, ix, iy, x-ix, y-iy, imarr[iy-1, ix],
                              imarr[iy, ix-1], imarr[iy, ix],
                              imarr[iy, ix+1], imarr[iy+1, ix], cp)
                    rows.append(my_row)
                except IndexError:
                    pass

    columns = 'amp ix iy xloc yloc p1 p3 p4 p5 p7 chiprob'.split()
    data_frame = pd.DataFrame(data=rows, columns=columns)
    log.info('number of data frame rows: %i' % len(data_frame))
    return data_frame

if __name__ == '__main__':
    import sys

    log.setLevel(logging.INFO)

#    infiles = sorted(glob.glob('data/E2V-CCD250-088_fe55_fe55_*.fits'))
#    nfiles = 10
#    dn_range = (760, 840)
#    df_file = 'cluster_pixel_data.pkl'

#    infiles = sorted(glob.glob('sim_data/000-01_fe55_fe55_*.fits'))
#    nfiles = 5
#    dn_range = (300, 370)
#    df_file = 'cluster_pixel_simdata.pkl'

    infiles = sorted(glob.glob('000-01_fe55_fe55_??_s1e-04.fits'))
    nfiles = 5
    dn_range = (300, 370)
    df_file = 'cluster_pixel_sim_s1e-04.pkl'

    generate_psf_catalogs(infiles, nfiles=nfiles)

    if not os.path.isfile(df_file):
        data_frame = create_data_frame(infiles, nfiles=nfiles,
                                       dn_range=dn_range)
        data_frame.to_pickle(df_file)
    else:
        data_frame = pd.read_pickle(df_file)

    try:
        xlocmin, xlocmax = float(sys.argv[1]), float(sys.argv[2])
        ylocmin, ylocmax = float(sys.argv[3]), float(sys.argv[4])
    except IndexError:
        xlocmin, xlocmax = -0.5, 0.5
        ylocmin, ylocmax = -0.5, 0.5

    chiprobmin = data_frame['chiprob'] > 0.01

    # plot the peak locations
    peak_selection = (chiprobmin
                      & (data_frame['xloc'] > xlocmin)
                      & (data_frame['xloc'] < xlocmax)
                      & (data_frame['yloc'] > ylocmin)
                      & (data_frame['yloc'] < ylocmax))
    peaks = data_frame[peak_selection]
    plt.rcParams['figure.figsize'] = 3, 3
    fig = plt.figure()
    plt.plot(peaks['xloc'], peaks['yloc'], 'k.')
    plt.axis((-0.6, 0.6, -0.6, 0.6))

    # Plot the distributions
    amps = range(1, 17)
    pixels = ('p1', 'p3', 'p5', 'p7')
    colors = ('blue', 'red', 'red', 'blue')
    linestyles = ('solid', 'solid', 'dashed', 'dashed')
    plt.rcParams['figure.figsize'] = 16, 12
    fig = plt.figure()
    range_ = (-10, 30)
    bins = 20
    for amp in amps:
        extname = 'AMP%02i' % amp
        panel = 2*(amp-1) + 1
        subplot = (4, 8, panel)
        axes = fig.add_subplot(*subplot)
        selection = (data_frame['amp'] == amp) & peak_selection
        df = data_frame[selection]
        #log.info("df size after selection: %i" % len(df))
        hists = []
        for pixel, color, linestyle in zip(pixels, colors, linestyles):
            if min(df[pixel].values) < range_[1]:
                hists.append(plt.hist(df[pixel].values, color=color,
                                      linestyle=linestyle, range=range_,
                                      bins=bins, histtype='step'))
        axes.set_title(extname)

        # Plot the residiual trailing distributions:
        # p5-p3 (serial, red), p7-p1 (parallel, blue)
        subplot = (4, 8, panel+1)
        axes = fig.add_subplot(*subplot)
        p1, p3, p5, p7 = tuple(hists)
        xvals = (p1[1][1:] + p1[1][:-1])/2.
        plt.errorbar(xvals, p5[0] - p3[0], yerr=np.sqrt(p5[0] + p3[0]),
                     marker='.', color='red')
        plt.errorbar(xvals, p7[0] - p1[0], yerr=np.sqrt(p7[0] + p1[0]),
                     marker='.', color='blue')
        axes.set_title(extname)

    amp = 1
    selection = (
        peak_selection
        & (data_frame['amp'] == amp)
#        & (data_frame['p5'] > 5)
#        & (data_frame['p5'] < 25)
        )
    df = data_frame[selection]
    plt.rcParams['figure.figsize'] = 5, 5
    fig = plt.figure()
    ax = fig.add_subplot(111)
    profile_plot(ax, df['ix'].values, df['p5'], color='red')
    profile_plot(ax, df['ix'].values, df['p3'], color='blue')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.loglog(df['p3'], df['p5'], 'k.')

    amp = 5
    selection = (
        peak_selection
        & (data_frame['amp'] == amp)
#        & (data_frame['p5'] > 5)
#        & (data_frame['p5'] < 25)
        )
    df = data_frame[selection]
    plt.rcParams['figure.figsize'] = 5, 5
    fig = plt.figure()
    ax = fig.add_subplot(111)
    profile_plot(ax, df['ix'].values, df['p5'], color='red')
    profile_plot(ax, df['ix'].values, df['p3'], color='blue')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.loglog(df['p3'], df['p5'], 'k.')
