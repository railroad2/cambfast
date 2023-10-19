# cambfast

The cambfast package computes cmb angular power spectra for the given variables
by linear interpolations between precomputed spectra. 
The camb package (https://camb.info/) is used to compute the precomputed spectra.

## usage

* Generating precomputed spectra.

```
import cambfast

cf = cambfast.CAMBfast()
cf.add_parameter('tau', 0.05, 0.02, 0.1, 10, fixed=False)
lmax = 30
cf.generate_interp(lmax, CMB_unit='muK')
cf.write_funcs(<precomputed file name>)
```

* Compute the spectra using the precomputed spectra.

```
cf = cambfast.CAMBfast(filename=<precomputed file name>)
print(cf.pname)
dls = cf.get_spectrum(tau=0.05)
```

* See `test_cambfast.py` for more.
    
The generated precomputed files are written in the 'precomputed' directory. 
Please note that the parameters other than the added parameters are 
taken from the Plank 2018 result (Planck collaboration, A&A 641, A6 (2020)).

