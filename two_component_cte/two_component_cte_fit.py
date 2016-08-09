import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
plt.ion()

class TrailedCharge(object):
    def __init__(self, ccd, amp, lastskip=4, gain=4):
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

    def cti(self, oscan_cols=2):
        """
        CTI estimate based on section 9.4, LCA-10103.
        """
        bias_est = np.mean(self.oscan_values[oscan_cols:])
        So = sum(self.oscan_values[:oscan_cols] - bias_est)
        Si = self.q_lastcol - bias_est
        return So/Si/self.ncols

    def model(self, oscan_pix, q0, tau, cti, bias_level):
        """
        Model of overscan pixel values:
        Exponential + CTI + bias level.
        """
        return (q0*self.nrows*np.exp(-oscan_pix/tau)
                + self.q_lastcol*self.ncols*cti**oscan_pix*(1. - cti)
                + bias_level*self.nrows)

    def resids(self, pars, pixels, values, errors):
        """
        (Data - model)/errors for chi-square calculation.
        """
        q0, tau, cti, bias_level = pars
        return (values - self.model(pixels, q0, tau, cti, bias_level))/errors

    def __call__(self, pars):
        """
        Return the chi-square as the objective function to minimize.
        """
        my_resids = self.resids(pars, self.oscan_pixels, self.oscan_values,
                                self.oscan_errors)
        return np.sum(my_resids**2)

    def plot_model(self, pars, color='blue', marker=':', fmt='none'):
        """
        Plot the model fit.
        """
        q0, tau, cti, bias_level = pars
        plt.plot(self.oscan_pixels,
                 self.model(self.oscan_pixels, q0, tau, cti, bias_level)/self.nrows,
                 marker, color=color)
        return plt.errorbar(self.oscan_pixels,
                            self.model(self.oscan_pixels, q0, tau, cti,
                                       bias_level)/self.nrows,
                            color=color, fmt=fmt)

    def plot_fit(self, pars, color='blue'):
        """
        Plot the overscan column sums (ADU / pixel) vs overscan pixel.
        """
        q0, tau, cti, bias_level = pars
#        fig = plt.figure()
        handle = plt.errorbar(self.oscan_pixels, self.oscan_values/self.nrows,
                              yerr=self.oscan_errors/self.nrows,
                              fmt='.', color=color)
        self.plot_model(pars, color=color)
        plt.xlabel('overscan pixel')
        plt.ylabel('ADU / pixel')
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
    sensor_id = 'ITL-3800C-013'
    ccd_high = sensorTest.MaskedCCD('ITL-3800C-013_superflat_high.fits')
    ccd_low = sensorTest.MaskedCCD('ITL-3800C-013_superflat_low.fits')

#    fig = plt.figure()
#    for amp in ccd_high:
    for amp in (3, 8):
        tc_low = TrailedCharge(ccd_low, amp, lastskip=4)
        tc_high = TrailedCharge(ccd_high, amp, lastskip=4)

        # Initial parameter estimates
        bias_level = 1000.
        cti_0 = tc_high.cti()
        tau = 1./np.log((tc_low.oscan_values[1] - bias_level*tc_low.nrows)
                        /(tc_low.oscan_values[2]- bias_level*tc_low.nrows))
        q0 = (tc_low.oscan_values[1]/tc_low.nrows - bias_level)/np.exp(-2/tau)

        p0 = q0, tau, cti_0, bias_level
        bounds = ((0, None), (0, None), (0, None), (0, None))

        tc_combined = MultiObjectiveFunctions([tc_low, tc_high])
        result = scipy.optimize.minimize(tc_combined, p0, method='L-BFGS-B',
                                         bounds=bounds)

        result_high = scipy.optimize.minimize(tc_high, p0, method='L-BFGS-B',
                                              bounds=bounds)
        result_low = scipy.optimize.minimize(tc_low, p0, method='L-BFGS-B',
                                             bounds=bounds)
        print amp, result_low.fun, result_low.x, '%.4e' % tc_low.cti()
        print amp, result_high.fun, result_high.x, '%.4e' % tc_high.cti()
        print amp, result.fun, result.x, '%.4e' % tc_high.cti()
        print
#        ax = fig.add_subplot(4, 4, amp)
        fig = plt.figure()
        handles = [tc_low.plot_fit(result_low.x)]
        handles.append(tc_high.plot_fit(result_high.x, color='red'))
        handles.append(tc_high.plot_model(result.x, color='green', marker='--',
                                          fmt='.'))
        if tc_high.oscan_values[0]/tc_low.oscan_values[0] > 5:
            plt.yscale('log')
        plt.title('%s, amp %i' % (sensor_id, amp))
        plt.legend(handles, ['low flux', 'high flux', 'joint fit'], loc=0)
        plt.savefig('%s_amp_%02i_overscan_fits.png' % (sensor_id, amp))
#    plt.savefig('%s_overscan_fits.png' % sensor_id)
