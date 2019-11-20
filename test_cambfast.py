import os.path 
import unittest
import time

import numpy as np
import pylab as plt

import cambfast
from cambfast import get_spectrum_camb

class test_cambfast_H0(unittest.TestCase): 
    @classmethod
    def setUpClass(cls):
        cls.lmax = 2000
        cls.cf = cambfast.CAMBfast('H0', 50, 80, nsample=30, lmax=cls.lmax)

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
    @classmethod
    def setUpClass(cls): 
        cls.lmax = 2000
        cls.pname = 'tau'
        cls.pmin = 0.030
        cls.pmax = 0.080
        cls.cf = cambfast.CAMBfast(cls.pname, cls.pmin, cls.pmax, nsample=30, lmax=cls.lmax, r=0.1)

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
            dls_camb = get_spectrum_camb(lmax=self.lmax, r=0.01, **kwargs)
            plt.loglog(dls_camb.T, 'b-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


class test_cambfast_tau_withfile(unittest.TestCase): 
    @classmethod
    def setUpClass(cls):
        cls.lmax = 2000
        cls.pname = 'tau'
        cls.pmin = 0.020
        cls.pmax = 0.080
        cls.fname = './tau_test.npz'
        cls.nsample = 30

    def setUp(self):
        self.t0 = time.time()

    def test1_writefile(self):
        if os.path.isfile(self.fname):
            pass
        else:
            cf = cambfast.CAMBfast(self.pname, self.pmin, self.pmax, nsample=self.nsample, lmax=self.lmax, r=0.01, CMB_unit='muK')
            cf.write_funcs('tau_test.npz')
        
    def test2_inter(self):
        cf = cambfast.CAMBfast(filename=self.fname)
        par_in = np.linspace(0.04, 0.07, 20)
        for i in par_in:
            dls_inter = cf.get_spectrum(i, lmax=1000)
            plt.loglog(dls_inter.T, 'r-')

    def testz_show(self):
        plt.show()

    def tearDown(self):
        print ("Elapsed time =", time.time() - self.t0)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(test_cambfast_tau_withfile))
    #suite.addTest(unittest.makeSuite(test_cambfast_tau))
    return suite


if __name__=='__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())

