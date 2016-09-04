import numpy as np
import scipy.optimize
import matplotlib
#matplotlib.use('cairo')
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
from MultiPanelPlot import MultiPanelPlot
plt.ion()

class TrailedCharge(object):
    def __init__(self, ccd, amp, lastskip=4, gain=4, bias_per_pixel=None):
        self.nrows = ccd.amp_geom.imaging.getMaxY() + 1
        mi = ccd[amp]
        imarr = mi.Factory(mi, ccd_high.amp_geom.imaging).getImage().getArray()
        oscan = mi.Factory(mi, ccd_high.amp_geom.serial_overscan).getImage().getArray()[:self.nrows, :-lastskip]

        self.q_lastcol = np.sum(imarr[:, -1])
        self.ncols = imarr.shape[1]
        self.oscan_values = np.array([np.sum(column) for column
                                      in oscan.transpose()])
        self.oscan_errors = np.sqrt(self.oscan_values*gain)/gain
        self.oscan_pixels = np.arange(1, len(self.oscan_values) + 1)
        self.bias_per_pixel = bias_per_pixel

    def cti(self, oscan_cols=2):
        """
        CTI estimate based on section 9.4, LCA-10103.
        """
        bias_est = np.mean(self.oscan_values[oscan_cols:])
        So = sum(self.oscan_values[:oscan_cols] - bias_est)
        Si = self.q_lastcol - bias_est
        return So/Si/self.ncols

    def model(self, oscan_pix, q0, tau, cti, bias_per_pix):
        """
        Model of overscan pixel values:
        Exponential + CTI + bias level.
        """
        bias = bias_per_pix*self.nrows
        return (q0*(self.q_lastcol - bias)*self.nrows*np.exp(-oscan_pix/tau)
                + (self.q_lastcol - bias)*self.ncols*cti**oscan_pix*(1. - cti)
                + bias)

    def resids(self, pars, pixels, values, errors):
        """
        (Data - model)/errors for chi-square calculation.
        """
        try:
            q0, tau, cti, bias_level = pars
        except ValueError:
            q0, tau, cti = pars
            bias_level = self.bias_per_pixel
        return (values - self.model(pixels, q0, tau, cti, bias_level))/errors

    def __call__(self, pars):
        """
        Return the chi-square as the objective function to minimize.
        """
        my_resids = self.resids(pars, self.oscan_pixels, self.oscan_values,
                                self.oscan_errors)
        return np.sum(my_resids**2)

    def plot_model(self, pars, color='blue', marker=':', label=None,
                   linewidth=1):
        """
        Plot the model fit.
        """
        try:
            q0, tau, cti, bias_level = pars
        except ValueError:
            q0, tau, cti = pars
            bias_level = self.bias_per_pixel
        return plt.plot(self.oscan_pixels,
                        self.model(self.oscan_pixels, q0, tau,
                                   cti, bias_level)/self.nrows,
                        marker, color=color, label=label, linewidth=linewidth)

    def plot_fit(self, pars, color='blue', label=None, linewidth=1):
        """
        Plot the overscan column sums (ADU / pixel) vs overscan pixel.
        """
        handle = plt.errorbar(self.oscan_pixels, self.oscan_values/self.nrows,
                              yerr=self.oscan_errors/self.nrows,
                              fmt='.', color=color)
        self.plot_model(pars, color=color, label=label, linewidth=linewidth)
        return handle

class MultiObjectiveFunctions(object):
    """
    Class to combine objective functions.
    """
    def __init__(self, funcs):
        self.funcs = funcs
    def __call__(self, pars):
        """
        Sum over the objective functions evaluated at the passed parameters.
        """
        return np.sum(func(pars) for func in self.funcs)

if __name__ == '__main__':
#    sensor_id = 'ITL-3800C-013'
#    sensor_id = 'ITL-3800C-022'
#    sensor_id = 'ITL-3800C-098'
    sensor_id = 'ITL-3800C-058'
    ccd_high = sensorTest.MaskedCCD('%s_superflat_high.fits' % sensor_id)
    ccd_low = sensorTest.MaskedCCD('%s_superflat_low.fits' % sensor_id)

    my_plots = MultiPanelPlot(4, 4, figsize=(12, 12))
    for amp in ccd_high:
        tc_low = TrailedCharge(ccd_low, amp, lastskip=4)
        tc_high = TrailedCharge(ccd_high, amp, lastskip=4)

        # Fix bias levels to last column value
        tc_low.bias_per_pixel = tc_low.oscan_values[-1]/tc_low.nrows
        tc_high.bias_per_pixel = tc_high.oscan_values[-1]/tc_high.nrows

        # Initial parameter estimates
        cti_0 = tc_high.cti()
        tau = 1./np.log((tc_high.oscan_values[1] - tc_high.oscan_values[-1])
                        /(tc_high.oscan_values[2] - tc_high.oscan_values[-1]))
        q0 = ((tc_high.oscan_values[1]/tc_high.nrows - tc_high.bias_per_pixel)
              /np.exp(-2/tau)/tc_high.q_lastcol)
        print amp, q0, tau

        p0 = q0, tau, cti_0
        bounds = ((0, None), (0, None), (0, None))
        result_high = scipy.optimize.minimize(tc_high, p0, method='L-BFGS-B',
                                              bounds=bounds)
        result_low = scipy.optimize.minimize(tc_low, p0, method='L-BFGS-B',
                                             bounds=bounds)

        tc_combined = MultiObjectiveFunctions([tc_low, tc_high])
        result = scipy.optimize.minimize(tc_combined, p0,
                                         method='L-BFGS-B', bounds=bounds)
        print amp, result.fun, result.x, '%.4e' % tc_high.cti()
        print amp, result_low.fun, result_low.x, '%.4e' % tc_low.cti()
        print amp, result_high.fun, result_high.x, '%.4e' % tc_high.cti()
        print
        my_plots.add_subplot()
        tc_low.plot_fit(result_low.x, label='low flux-only', linewidth=1)
        tc_high.plot_fit(result_high.x, color='red', label='high flux-only',
                         linewidth=1)
        tc_low.plot_model(result.x, color='green', marker='--',
                          label='joint')
        tc_high.plot_model(result.x, color='green', marker='--')
        if tc_high.oscan_values[0]/tc_low.oscan_values[0] > 5:
            plt.yscale('log')
        axis_range = list(plt.axis())
        axis_range[:2] = 0.5, axis_range[1] + 0.5
        plt.axis(axis_range)
        plt.legend(loc=0, fontsize='x-small')
    my_plots.set_title(sensor_id)
    my_plots.set_xlabel('overscan pixel')
    my_plots.set_ylabel('ADU / pixel')
    plt.savefig('%s_overscan_fits.png' % sensor_id)
