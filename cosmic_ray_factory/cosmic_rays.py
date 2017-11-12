from collections import namedtuple, defaultdict
import numpy as np
import numpy.random as random
import astropy.io.fits as fits

CR_Span = namedtuple('CR_Span', 'x0 y0 pixel_values'.split())

class CosmicRays(list):
    def __init__(self, *args, **kwds):
        super(CosmicRays, self).__init__(*args, **kwds)

    def paint(self, image_array, num_crs=30):
        for i in range(num_crs):
            image_array = self.paint_cr(image_array)
        return image_array

    def paint_cr(self, image_array, index=None, pixel=None):
        if index is None:
            cr = random.choice(self)
        else:
            cr = self[index]
        if pixel is None:
            pixel = (random.randint(image_array.shape[1]),
                     random.randint(image_array.shape[0]))
        for span in cr:
            for dx, value in enumerate(span.pixel_values):
                try:
                    image_array[pixel[1] + span.y0 - cr[0].y0,
                                pixel[0] + span.x0 - cr[0].x0 + dx] += value
                except IndexError:
                    pass
        return image_array

    def read_catalog(self, catalog_file, extname='COSMIC_RAYS'):
        with fits.open(catalog_file) as catalog:
            cr_cat = catalog[extname]
            crs = defaultdict(list)
            for i, span in enumerate(cr_cat.data):
                crs[span[0]].append(CR_Span(*tuple(span)[1:]))
        self.extend(crs.values())

if __name__ == '__main__':
    crs = CosmicRays()
    crs.read_catalog('cosmic_ray_catalog.fits')

    imarr = np.zeros((2000, 512), dtype=np.int)

    imarr = crs.paint(imarr, num_crs=100)

    foo = fits.HDUList([fits.PrimaryHDU()])
    foo[0].data = imarr
    foo.writeto('foo.fits', overwrite=True)
