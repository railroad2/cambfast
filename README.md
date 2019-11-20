cambfast

A package to compute cmb angular power spectra for a single variable 
by linear interpolations between CAMB spectra.

usage

    import cambfast

    cf = cambfast.Cambfast('tau', 0.03, 0.08, nsample=10, lmax=2000, As=2.092e-9, r=0.01)
    dls = cf.get_spectrum(0.05)

or

    cf = cambfast.Cambfast(filename=<pre-computed file name>)
    print(cf.pname)
    dls = cf.get_spectrum(0.05)
    
It generated precomputed files in 'precomputed' directory. 
Please note that the parameters other than the main parameter are default values
unless they are not defined as extra arguments.
    
     

