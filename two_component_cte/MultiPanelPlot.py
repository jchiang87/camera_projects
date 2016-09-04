import numpy as np
import matplotlib.pyplot as plt

class MultiPanelPlot(object):
    def __init__(self, nx, ny, figsize=(10, 10),
                 subplot_offsets=(0.025, 0.025)):
        plt.rcParams['figure.figsize'] = figsize
        self.fig = plt.figure()
        self.frame_axes = self.fig.add_subplot(111, frameon=False)
        self.frame_axes.get_xaxis().set_ticks([])
        self.frame_axes.get_yaxis().set_ticks([])
        self.nx, self.ny = nx, ny
        self.counter = 0
        self.subplot_offsets = subplot_offsets
        self.subplot_axes = []
    def add_subplot(self):
        self.counter += 1
        if self.counter > self.nx*self.ny:
            raise RuntimeError("No more subplot slots are available.")
        axes = self.fig.add_subplot(self.nx, self.ny, self.counter)
        bbox = axes.get_position()
        points = bbox.get_points()
        points[0] += self.subplot_offsets[0]
        points[1] += self.subplot_offsets[1]
        bbox.set_points(points)
        axes.set_position(bbox)
        self.subplot_axes.append(axes)
    def set_xlabel(self, *args, **kwds):
        self.frame_axes.set_xlabel(*args, **kwds)

    def set_ylabel(self, *args, **kwds):
        self.frame_axes.set_ylabel(*args, **kwds)

    def set_title(self, title):
        self.frame_axes.text(0.5, 1. + 3*self.subplot_offsets[1],
                             title,
                             horizontalalignment='center',
                             verticalalignment='top',
                             transform=self.frame_axes.transAxes,
                             size='large')

if __name__ == '__main__':
    plt.ion()
    my_plots = MultiPanelPlot(3, 3)
    x = np.linspace(1, 100, 20)
    for i in range(my_plots.nx*my_plots.ny):
        my_plots.add_subplot()
        plt.errorbar(x, np.random.normal(1, size=len(x)))
    my_plots.set_xlabel('x values')
    my_plots.set_ylabel('y values')
    my_plots.set_title('my plots')


