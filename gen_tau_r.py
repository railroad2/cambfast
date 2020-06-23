import time

import numpy as np

import cambfast
from cambfast import get_spectrum_camb

lmax = 40
npts = 30

cf = cambfast.CAMBfast()
cf.ini_file = '../ini/planck_2018.ini'
cf.add_parameter('As', 2.0e-9, 1.7e-9, 2.3e-9, npts, fixed=False)
cf.add_parameter('tau', 0.05, 0.02, 0.08, npts, fixed=False)
cf.add_parameter('r', 0.05, -0.3, 0.3, npts, fixed=True)
cf.generate_interp(lmax, CMB_unit='muK')
cf.write_funcs(f'As_tau_r_lmax{lmax}_npts{npts}.npz')

