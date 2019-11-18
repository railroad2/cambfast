import time
import numpy as np
import camb
from scipy.interpolate import interp2d

from gbpipe.utils import dl2cl

class cambfast():

    def __init__(self, pname=None, pmin=None, pmax=None, lmax=None, nsample=10, CMB_unit=None, filename=None):
        self.funcTT = []
        self.funcEE = []
        self.funcBB = []
        self.funcTE = []
        self.respar = []

        if filename is not None:
            self.load_funcs(filename)

        else:
            if (np.array([pname, pmin, pmax, lmax]) == None).any():
                raise ValueError("All the arguments: pname, pmin, pmax, lmax should be given.")

            self.pname = pname
            self.pmin = pmin
            self.pmax = pmax
            self.nsample = nsample
            self.lmax = lmax
            self.CMB_unit = CMB_unit

            self.__generate_interp()

    def __generate_interp(self):
        pars = camb.CAMBparams()
        kwargs_cosmology = {}
        if self.pname != 'H0':
            kwargs_cosmology['H0'] = 67.5

        self.ell = np.arange(self.lmax+1)
        parr = np.linspace(self.pmin, self.pmax, self.nsample)

        dls_TT = []
        dls_EE = []
        dls_BB = []
        dls_TE = []
        for par in parr:
            kwargs_cosmology[self.pname] = par
            pars.set_cosmology(**kwargs_cosmology)
            pars.InitPower.set_params()
            pars.WantTensors=True
            results = camb.get_results(pars)
            dls = results.get_total_cls(lmax=self.lmax, CMB_unit=self.CMB_unit).T

            dls_TT.append(dls[0])
            dls_EE.append(dls[1])
            dls_BB.append(dls[2])
            dls_TE.append(dls[3])

        dls_TT = np.array(dls_TT)
        dls_EE = np.array(dls_EE)
        dls_BB = np.array(dls_BB)
        dls_TE = np.array(dls_TE)

        self.funcTT = interp2d(self.ell, parr, dls_TT) 
        self.funcEE = interp2d(self.ell, parr, dls_EE) 
        self.funcBB = interp2d(self.ell, parr, dls_BB) 
        self.funcTE = interp2d(self.ell, parr, dls_TE) 

    def get_spectrum(self, par, lmax=None, isDl=True):
        dls = []
        if lmax is None:
            lmax = self.lmax
        elif lmax > self.lmax:
            lmax = self.lmax

        ell = np.arange(lmax+1)
        dls.append(self.funcTT(ell, par))
        dls.append(self.funcEE(ell, par))
        dls.append(self.funcBB(ell, par))
        dls.append(self.funcTE(ell, par))

        dls = np.array(dls)

        if isDl:
            return dls
        else:
            cls = dl2cl(dls)
            return cls

    def write_funcs(self, fname):
        np.savez(fname, 
                 pname = self.pname, 
                 pmin = self.pmin,
                 pmax = self.pmax,
                 nsample = self.nsample,
                 lmax = self.lmax,
                 CMB_unit = self.CMB_unit,
                 funcTT = self.funcTT, 
                 funcEE = self.funcEE, 
                 funcBB = self.funcBB, 
                 funcTE = self.funcTE, 
                 allow_pickle = True)
        print ("The functions have been saved in {}".format(fname))

    def load_funcs(self, fname):
        dat = np.load(fname, allow_pickle=True)
        self.pname = dat['pname']
        self.pmin = dat['pmin']
        self.pmax = dat['pmax']
        self.nsample = dat['nsample']
        self.lmax = dat['lmax']
        self.CMB_unit = dat['CMB_unit']
        self.funcTT = dat['funcTT'].item()
        self.funcEE = dat['funcEE'].item()
        self.funcBB = dat['funcBB'].item()
        self.funcTE = dat['funcTE'].item()

