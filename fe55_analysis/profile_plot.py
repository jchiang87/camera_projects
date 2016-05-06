import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

__all__ = ['profile_plot']

def profile_plot(ax, x, y, bins=20, xbounds=None, marker='.',
                 color='black'):
    if xbounds is None:
        xmin, xmax = min(x), max(x)
    else:
        xmin, xmax = xbounds
    binedges = np.linspace(xmin, xmax, bins+1)

    df = pd.DataFrame(dict(x=x, y=y))
    df['bin'] = np.digitize(df['x'], binedges)
    bincenters = (binedges[1:] + binedges[:-1])/2.
    profile = pd.DataFrame({'bincenters' : bincenters,
                            'N' : df['bin'].value_counts(sort=False)},
                           index=range(1, bins+1))
    for bin in profile.index:
        profile.ix[bin, 'ymean'] = df.ix[df['bin']==bin, 'y'].mean()
        profile.ix[bin, 'ymeanerr'] = df.ix[df['bin']==bin, 'y'].std()

    ax.errorbar(profile['bincenters'], profile['ymean'],
                yerr=profile['ymeanerr'], xerr=(xmax-xmin)/bins/2,
                fmt=None, marker=marker, ecolor=color)
    return profile

if __name__ == '__main__':
    plt.ion()

    x = np.linspace(0, 2*np.pi, 1000)
    amp = 0.5
    offset = 2.
    ymodel = amp*np.sin(x) + offset
    ysig = 0.3
    y = np.random.normal(ymodel, ysig)

    ax = plt.figure().add_subplot(111)
    profile_plot(ax, x, y, bins=30)

