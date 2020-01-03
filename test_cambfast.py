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
        cls.cf.add_parameter('tau', 0.0522, 0.02, 0.1, 10, fixed=False)
        cls.cf.add_parameter('r', 0.01, 0.0, 0.1, 10, fixed=False)
        cls.cf.add_parameter('As', 2e-9, 1e-9, 3e-9, 10, fixed=False)
        cls.lmax = 60
        cls.cf.generate_interp(cls.lmax, CMB_unit='muK') 
        cls.cf.write_funcs('tau_r_As_lmax60.npz')

    def setUp(self):
        self.t0 = time.time()

    def test1_cambfast_tau(self):
        par_in = np.linspace(0.03, 0.04, 50)
        for par in par_in:
            cls = self.cf.get_spectrum(tau=par, isDl=True)
            plt.loglog(cls.T, 'r-')

    def test2_cambfast_As(self):
        par_in = np.linspace(1e-9, 3e-9, 50)

        for par in par_in:
            cls = self.cf.get_spectrum(H0=par, isDl=True)
            plt.loglog(cls.T, 'g-')

    """
    def test3_camb_single(self):
        par_in = np.linspace(0.03, 0.04, 50)
        kwargs = {}
        for i in par_in:
            kwargs['tau'] = i
            dls_camb = get_spectrum_camb(lmax=self.lmax, r=0.01, **kwargs, isDl=True, CMB_unit='muK')
            plt.loglog(dls_camb.T, 'b-')
    """

    def testz_show(self):
        #plt.show()
        pass

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


class test_cambfast_precomputed(unittest.TestCase): 
    @classmethod
    def setUpClass(cls):
        cls.cf = cambfast.CAMBfast() 
        cls.cf.load_funcs('tau_r_As_lmax60.npz')
        cls.lmax = 60

    def setUp(self):
        self.t0 = time.time()

    def test1_cambfast_tau(self):
        par_in = np.linspace(0.03, 0.04, 10)
        for par in par_in:
            cls = self.cf.get_spectrum(tau=par, isDl=True)
            #plt.loglog(cls.T, 'r-')

    def test2_cambfast_As(self):
        par_in = np.linspace(1e-9, 3e-9, 10)

        for par in par_in:
            cls = self.cf.get_spectrum(As=par, isDl=True)
            #plt.loglog(cls.T, 'g-')

    def test2_1_cambfast_tau(self):
        par_in = np.linspace(0.03, 0.04, 10)
        for par in par_in:
            cls = self.cf.get_spectrum(tau=par, isDl=True)
            #plt.loglog(cls.T, 'r-')

    def test3_camb_single(self):
        par_in = np.linspace(0.03, 0.04, 10)
        kwargs = {}
        for i in par_in:
            kwargs['tau'] = i
            dls_camb = get_spectrum_camb(lmax=self.lmax, r=0.01, **kwargs, isDl=True, CMB_unit='muK')
            #plt.loglog(dls_camb.T, 'b-')

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

