import os.path 
import unittest
import time

import numpy as np
import pylab as plt

import cambfast
from gbpipe.spectrum import get_spectrum_camb

class test_cambfast_H0(unittest.TestCase): 

    def __init__(self):
        self.lmax = 2000
        self.cf = cambfast.cambfast('H0', 50, 80, nsample=30, lmax=self.lmax)

    def setUp(self):
        self.t0 = time.time()

    def test1_inter(self):
        par_in = np.linspace(53, 79, 10)
        for i in par_in:
            dls_inter = self.cf.get_spectrum(i)
            plt.loglog(dls_inter.T, 'r-')

    def test2_camb(self):
        par_in = np.linspace(53, 79, 10)
        for i in par_in:
            dls_camb = get_spectrum_camb(lmax=self.lmax, H0=i)
            plt.loglog(dls_camb.T, 'b-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


class test_cambfast_tau(unittest.TestCase): 

    def __init__(self):
        self.lmax = 2000
        self.pname = 'tau'
        self.pmin = 0.030
        self.pmax = 0.080
        self.cf = cambfast.cambfast(self.pname, self.pmin, self.pmax, nsample=30, lmax=self.lmax)

    def setUp(self):
        self.t0 = time.time()

    def test1_inter(self):
        par_in = np.linspace(0.04, 0.07, 10)
        for i in par_in:
            dls_inter = self.cf.get_spectrum(i)
            plt.loglog(dls_inter.T, 'r-')

    def test2_camb(self):
        par_in = np.linspace(0.04, 0.07, 10)
        kwargs = {}
        for i in par_in:
            kwargs[self.pname] = i
            dls_camb = get_spectrum_camb(lmax=self.lmax, **kwargs)
            plt.loglog(dls_camb.T, 'b-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


class test_cambfast_tau_withfile(unittest.TestCase): 

    lmax = 100
    pname = 'tau'
    pmin = 0.030
    pmax = 0.080
    fname = './tau_test.npz'

    def setUp(self):
        self.t0 = time.time()

    def test1_writefile(self):
        if os.path.isfile(self.fname):
            pass
        else:
            cf = cambfast.cambfast(self.pname, self.pmin, self.pmax, nsample=30, lmax=self.lmax)
            cf.write_funcs('tau_test.npz')
        
    def test2_inter(self):
        cf = cambfast.cambfast(filename=self.fname)
        par_in = np.linspace(0.04, 0.07, 100)
        for i in par_in:
            dls_inter = cf.get_spectrum(i, lmax=47)
            plt.loglog(dls_inter.T, 'r-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(test_cambfast_tau_withfile))
    return suite


if __name__=='__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())

