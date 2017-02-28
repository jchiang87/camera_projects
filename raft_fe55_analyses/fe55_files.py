import os
from collections import OrderedDict

fits_files = OrderedDict()
sensors = OrderedDict()
with open('LCA-11021_RTM-004_ETU2-Dev_fe55_file_list.txt') as input:
    for line in input:
        slot, fits_file = line.strip().split()
        fits_files[slot] = fits_file
        sensors[slot] = os.path.basename(fits_file)[:len('ITL-3800C-133-Dev')]
