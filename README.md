# FRB Redshift Estimator

Code as described in [Pol, Lam, McLaughlin, Lazio, and Cordes (2019)](http://adsabs.harvard.edu/abs/2019arXiv190307630P). 

Our code can be run in one of two ways: from the command line and from directly within python (compatible with versions 2/3). The code also needs the NE2001 code (found [here](http://astrosun2.astro.cornell.edu/~cordes/NE2001/)) to be installed before the program can be executed.

## Command Line

To calculate the redshift of an FRB, 

	usage: frbz.py [-h] [--NE2001 NEDIR] [--unweighted] [--hostdm HOSTDM]
		       [--mwdm MWDM]
		       dm dmerr ...

	FRB Redshift Estimator

	positional arguments:
	  dm               Observed DM [pc cm^-3]
	  dmerr            Error on observed DM [pc cm^-3]
	  galcoord         If --mwdm is not provided, two values separated by a space:
		           Galactic latitude and Galactic longitude [deg]

	optional arguments:
	  -h, --help       show this help message and exit
	  --NE2001 NEDIR   Path pointing to the NE2001 bin.NE2001/ directory location
	  --unweighted     Use uniform weighted distribution (versus matter weighted
		           distribution
	  --hostdm HOSTDM  Host DM [pc cm^-3]
	  --mwdm MWDM      Milky Way DM [pc cm^-3]

## Python function

calcz(dm,dmerr,mwarg,hostdm=0.0,weighted=True,evaluate=True,NEDIR=\"NE2001/bin.NE2001/\"):

    Calculates the PDF for the redshift of the FRB
    dm        : DM value [pc cm^-3]
    dmerr     : DM error [pc cm^-3]
    mwarg     : Either:
              : tuple of (Galactic longitude, Galactic latitude) in [deg], or
              : Milky Way DM [pc cm^-3]

    weighted  : use matter weighted distribution if true
    evaluate  : if true, returns the DM value with minus and plus errors
    -----
    Returns
    * A PDF instance of the redshift if evaluate = False
    * The redshift value, the minus uncertainity, the plus uncertainty
    
