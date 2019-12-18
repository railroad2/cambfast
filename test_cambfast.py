import os.path 
import unittest
import time

import numpy as np
import pylab as plt

import cambfast
from cambfast import get_spectrum_camb


class test_cambfast(unittest.TestCase): 
    @classmethod
    def setUpClass(cls):
        cls.cf = cambfast.CAMBfast() 
        cls.cf.add_parameter('tau', 0.0522, 0.02, 0.1, 3, fixed=False)
        cls.cf.add_parameter('r', 0.01, 0.0, 0.1, 3, fixed=False)
        cls.cf.add_parameter('H0', 55., 50., 70., 3, fixed=False)
        cls.cf.generate_interp(1000, CMB_unit='muK') 
        cls.lmax = 2000
        cls.cf.write_funcs('test.npz')

    def setUp(self):
        self.t0 = time.time()

    def test1_cambfast_tau(self):
        par_in = np.linspace(0.03, 0.04, 10)
        for par in par_in:
            cls = self.cf.get_spectrum(tau=par, isDl=True)
            plt.loglog(cls.T, 'r-')

    def test2_cambfast_H0(self):
        par_in = np.linspace(55, 67, 10)

        for par in par_in:
            cls = self.cf.get_spectrum(H0=par, isDl=True)
            plt.loglog(cls.T, 'g-')

    def test3_camb_single(self):
        par_in = np.linspace(0.03, 0.04, 10)
        kwargs = {}
        for i in par_in:
            kwargs['tau'] = i
            dls_camb = get_spectrum_camb(lmax=self.lmax, r=0.01, **kwargs, isDl=True, CMB_unit='muK')
            plt.loglog(dls_camb.T, 'b-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


class test_cambfast_precomputed(unittest.TestCase): 
    @classmethod
    def setUpClass(cls):
        cls.cf = cambfast.CAMBfast() 
        cls.cf.load_funcs('test.npz')
        cls.lmax = 2000

    def setUp(self):
        self.t0 = time.time()

    def test1_cambfast_tau(self):
        par_in = np.linspace(0.03, 0.04, 10)
        for par in par_in:
            cls = self.cf.get_spectrum(tau=par, isDl=True)
            plt.loglog(cls.T, 'r-')

    def test2_cambfast_H0(self):
        par_in = np.linspace(55, 67, 10)

        for par in par_in:
            cls = self.cf.get_spectrum(H0=par, isDl=True)
            plt.loglog(cls.T, 'g-')

    def test3_camb_single(self):
        par_in = np.linspace(0.03, 0.04, 10)
        kwargs = {}
        for i in par_in:
            kwargs['tau'] = i
            dls_camb = get_spectrum_camb(lmax=self.lmax, r=0.01, **kwargs, isDl=True, CMB_unit='muK')
            plt.loglog(dls_camb.T, 'b-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(test_cambfast_precomputed))
    #suite.addTest(unittest.makeSuite(test_cambfast))
    return suite


if __name__=='__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())

