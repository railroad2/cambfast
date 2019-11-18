cambfast

A package to compute cmb angular power spectra for a single variable 
by linear interpolations between CAMB spectra.

usage
    >>> import cambfast

    >>> cf = cambfast.cambfast('tau', 0.03, 0.08, nsample=10, lmax=2000, As=2.092e-9, r=0.01)
    >>> dls = cf.get_spectrum(0.05)

    or

    >>> cf = cambfast.cambfast(filename=<pre-computed file name>)
    >>> print(cf.pname)
    >>> dls = cf.get_spectrum(0.05)
    
    
     

