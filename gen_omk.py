import time

import numpy as np

import cambfast
from cambfast import get_spectrum_camb

lmax = 50
npts = 100

cf = cambfast.CAMBfast()
cf.add_parameter('omk', 0.05, -0.3, 0.3, npts, fixed=False)
cf.generate_interp(lmax, CMB_unit='muK')
cf.write_funcs(f'omk_lmax{lmax}.npz')

