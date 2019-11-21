import os
import time

import camb
import numpy as np
from scipy.interpolate import interp2d


args_cosmology = ['H0', 'cosmomc_theta', 'ombh2', 'omch2', 'omk', 
                  'neutrino_hierarchy', 'num_massive_nutrinos',
                  'mnu', 'nnu', 'YHe', 'meffsterile', 'standard_neutrino_neff', 
                  'TCMB', 'tau', 'deltazrei', 'bbnpredictor', 'theta_H0_range'] 

args_InitPower = ['As', 'ns', 'nrun', 'nrunrun', 'r', 'nt', 'ntrun', 'pivot_scalar', 
                  'pivot_tensor', 'parameterization']


class CAMBfast():

    def __init__(self, pname=None, pmin=None, pmax=None, lmax=None, nsample=10, CMB_unit=None, filename=None, **kwargs):
        self.funcTT = []
        self.funcEE = []
        self.funcBB = []
        self.funcTE = []
        self.CMB_unit = CMB_unit
        self.pothers = kwargs

        self.attrs = ['pname', 'pmin', 'pmax', 'nsample', 'lmax', 'CMB_unit', 
                      'funcTT', 'funcEE', 'funcBB', 'funcTE', 'pothers']

        self.fn_pre = './precomputed/{0}_{1}_{2}_{3}_{4}.npz'.format(pname, pmin, pmax, nsample, CMB_unit)

        if filename is not None:
            self.load_funcs(filename)

        elif os.path.isfile(self.fn_pre):
            self.load_funcs(self.fn_pre)

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
            try:
                os.mkdir('precomputed')
            except:
                pass

            self.write_funcs(self.fn_pre)


    def __generate_interp(self):
        self.ell = np.arange(self.lmax+1)
        parr = np.linspace(self.pmin, self.pmax, self.nsample)

        dls_TT = []
        dls_EE = []
        dls_BB = []
        dls_TE = []
        kwargs = self.pothers.copy()
        for par in parr:
            kwargs[self.pname] = par
            dls = get_spectrum_camb(lmax=self.lmax, CMB_unit=self.CMB_unit, **kwargs)

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

        kwargs = {}
        for aname in self.attrs:
            kwargs[aname] = getattr(self, aname)

        np.savez(fname, **kwargs)

        print ("The functions have been saved in {}".format(fname))

    def load_funcs(self, fname):
        dat = np.load(fname, allow_pickle=True)
        for name in dat.files:
            setattr(self, name, dat[name].item())


def dl2cl(dls):
    cls = dls.copy()

    if (len(cls) > 10):
        cls = cls.T

    ell = np.arange(len(cls[0]))

    for i in range(len(cls)):
        cls[i][1:] = cls[i][1:] * (2. * np.pi) / (ell[1:] * (ell[1:] + 1))

    if (len(dls) > 10):
        cls = cls.T

    return cls


def get_spectrum_camb(lmax, 
                      isDl=True, cambres=False, TTonly=False, unlensed=False, CMB_unit=None, 
                      **kwargs):
   
    ## arguments to dictionaries
    kwargs_cosmology={}
    kwargs_InitPower={}
    wantTensor = False

    for key, value in kwargs.items():  
        if key in args_cosmology: 
            kwargs_cosmology[key]=value
            if key == 'r':
                wantTensor = True
        elif key in args_InitPower:
            kwargs_InitPower[key]=value
        else:
            print_warning('Wrong keyword: ' + key)

    ## for camb > 1.0
    if not ('H0' in kwargs_cosmology.keys()):
        kwargs_cosmology['H0'] = 67.5

    ## call camb
    pars = camb.CAMBparams()
    pars.set_cosmology(**kwargs_cosmology)
    pars.InitPower.set_params(**kwargs_InitPower)
    pars.WantTensors = True
    results = camb.get_results(pars)

    if (TTonly):
        if unlensed:
            dls = results.get_unlensed_total_cls(lmax=lmax, CMB_unit=CMB_unit).T[0]
        else:
            dls = results.get_total_cls(lmax=lmax, CMB_unit=CMB_unit).T[0]
    else: 
        if unlensed:
            dls = results.get_unlensed_total_cls(lmax=lmax, CMB_unit=CMB_unit).T
        else:
            dls = results.get_total_cls(lmax=lmax, CMB_unit=CMB_unit).T

    if (isDl):
        res = dls
    else:
        cls = dl2cl(dls)
        res = cls

    if (cambres):
        return res, results
    else:
        return res

