# DM_IGM

Documentation in progress...

Our code can be run in one of two ways: from the command line and from directly within python (compatible with versions 2/3).

## Command Line

To calculate the redshift of an FRB, 

usage: frbz.py [-h] [--NE2001 NEDIR] [--unweighted] [--hostdm HOSTDM]
               gl gb dm dmerr

FRB Redshift Estimator

positional arguments:
  gl               Galactic latitude [deg]
  gb               Galactic longitude [deg]
  dm               Observed DM [pc cm$^{-3}$ ]
  dmerr            Error on observed DM [pc cm$^{-3}$ ]

optional arguments:
  -h, --help       show this help message and exit
  --NE2001 NEDIR   Directory pointing to NE2001's bin.NE2001/
  --unweighted     Use uniform weighted distribution (versus matter weighted
                   distribution
  --hostdm HOSTDM  Host DM [pc cm$^{-3}$]

## Python function

calcz(gl,gb,dm,dmerr,hostdm=0.0,weighted=True,evaluate=True,NEDIR=\"NE2001/bin.NE2001/\"):

    Calculates the PDF for the redshift of the FRB
    dm        : DM value [pc cm^-3]
    dmerr     : DM error [pc cm^-3]
    gl        : Galactic longitude [deg]
    gb        : Galactic latitude [deg]

    weighted  : use matter weighted distribution if true
    evaluate  : if true, returns the DM value with minus and plus errors
    -----
    Returns
    * A PDF instance of the redshift if evaluate = False
    * The redshift value, the minus uncertainity, the plus uncertainty
    