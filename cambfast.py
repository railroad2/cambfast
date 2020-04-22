import os
import time

import camb
import numpy as np
from scipy.interpolate import interp2d, RegularGridInterpolator


args_cosmology = ['H0', 'cosmomc_theta', 'ombh2', 'omch2', 'omk', 
                  'neutrino_hierarchy', 'num_massive_nutrinos',
                  'mnu', 'nnu', 'YHe', 'meffsterile', 'standard_neutrino_neff', 
                  'TCMB', 'tau', 'deltazrei', 'bbnpredictor', 'theta_H0_range'] 

args_InitPower = ['As', 'ns', 'nrun', 'nrunrun', 'r', 'nt', 'ntrun', 'pivot_scalar', 
                  'pivot_tensor', 'parameterization']


class CAMBfast():

    class Parameter(): 

        def __init__(self, pname, pinit, pmin, pmax, nstep, fixed=False):
            assert (pmax >= pmin),  f'pmax (={pmax}) must be larger than pmin (={pmin}).'
            self.name = pname
            self.init = pinit
            self.min = pmin
            self.max = pmax
            self.nstep = nstep
            self.fixed = fixed

    def __init__(self, filename=None):
        self.pars = []
        self.pfixed = []
        self.lmax = 2000 
        self.CMB_unit = 'muK'
        self.ell = []
        self.funcTT = []
        self.funcEE = []
        self.funcBB = []
        self.funcTE = []

        self.attrs = ['pars', 'pfixed', 'lmax', 'CMB_unit', 'ell',
                      'funcTT', 'funcEE', 'funcBB', 'funcTE']

        if filename is not None:
            self.load_funcs(filename)

    def add_parameter(self, pname, pinit, pmin, pmax, nstep, fixed=False):
        par = self.Parameter(pname, pinit, pmin, pmax, nstep, fixed)
        if par.fixed:
            self.pfixed.append(par)
        else:
            self.pars.append(par)

    def generate_interp(self, lmax, CMB_unit):
        print ('Generating interpolators')
        self.lmax = lmax
        self.CMB_unit = CMB_unit
        self.ell = np.arange(self.lmax+1)

        parr = []
        for par in self.pars:
            parr.append(np.linspace(par.min, par.max, par.nstep))

        kw_pfixed = {}
        for par in self.pfixed:
            kw_pfixed[par.name] = par.init

        grid = np.meshgrid(*parr, indexing='ij')
        shape = list(grid[0].shape)

        grids = np.reshape(grid, (len(grid), np.prod(shape))) 

        def __spectrum(*p):
            print (*p)
            kw_par = {} 
            for par, pval in zip(self.pars, p):
                kw_par[par.name] = pval

            return get_spectrum_camb(self.lmax, isDl=True, CMB_unit=self.CMB_unit, 
                                     **kw_par, **kw_pfixed) 

        dat = np.array(list(map(__spectrum, *grids)))
        newshape = np.append(shape, np.shape(dat)[-1])

        datTT = np.reshape(dat[:, 0, :], newshape)
        datEE = np.reshape(dat[:, 1, :], newshape)
        datBB = np.reshape(dat[:, 2, :], newshape)
        datTE = np.reshape(dat[:, 3, :], newshape)

        self.funcTT = RegularGridInterpolator((*parr, self.ell), datTT)
        self.funcEE = RegularGridInterpolator((*parr, self.ell), datEE)
        self.funcBB = RegularGridInterpolator((*parr, self.ell), datBB)
        self.funcTE = RegularGridInterpolator((*parr, self.ell), datTE)

    def get_spectrum(self, lmax=None, isDl=True, **kwpars):
        dls = []
        if lmax is None:
            lmax = self.lmax
        elif lmax > self.lmax:
            lmax = self.lmax

        pars = []
        for p in self.pars:
            try:
                pars.append(kwpars[p.name])
            except KeyError:
                pars.append(p.init)

        ell = np.arange(lmax+1)
        try:
            dls.append(self.funcTT((*pars, ell))[:lmax+1])
            dls.append(self.funcEE((*pars, ell))[:lmax+1])
            dls.append(self.funcBB((*pars, ell))[:lmax+1])
            dls.append(self.funcTE((*pars, ell))[:lmax+1])

            dls = np.array(dls)
        except ValueError:
            dls = get_spectrum_camb(lmax=lmax, isDl=True, **kwpars)

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
            if name[:4]=='func':
                setattr(self, name, dat[name].item())
            else:
                setattr(self, name, dat[name])


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
            print('Wrong keyword: ' + key)

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


