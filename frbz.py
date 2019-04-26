#!/usr/bin/python
import numpy as np
import glob
import pdf
PDF = pdf.PDF
import scipy.interpolate as interpolate
import NE2001
import argparse
import sys

#Test addition of pyymw16
import pyymw16



# ==================================================
# From a pdf, get the mode and innter 95% confidence intervals
# ==================================================
def getbounds(p):
    p.resample(1000)# enforce linear bins for proper cdf-ing

    mode,b,c,d = p.likelihood_evaluator(median=False,pm=False)
    xm,xp = p.likelihood_evaluator(values=[0.025,0.975])

    return mode,mode-xm,xp-mode

# ==================================================
# Other utility functions
# ==================================================
def gaussian(x,b,c):
    a = 1.0/np.sqrt(2*np.pi*c**2) #normalized, though pdf will renormalize anyway
    return a*np.exp(-0.5*(x-b)**2 / c**2)

def normalize_area(array,x=None,full=False):
    if x is None:
        x=np.arange(len(array))
    area=np.trapz(array,x=x)
    if full:
        return array/area,area
    return array/area



# ==================================================
# Primary function
# ==================================================
def calcz(dm,dmerr,mwarg,hostdm=0.0,weighted=True,evaluate=True,NEDIR="NE2001/bin.NE2001/", ymw = False):
    """
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
    """
    
    if weighted:
        filenames = sorted(glob.glob("slices_weighted/*npz"))
    else:
        filenames = sorted(glob.glob("slices_unweighted/*npz"))

    '''
    --------------------------------------------------
    Read in integrated redshift slices
    --------------------------------------------------
    '''
    zs = np.zeros(len(filenames))
    pdfs = []

    for i,filename in enumerate(filenames):
        npz = np.load(filename)
        z = filename.split("/")[-1][2:-4]
        zs[i] = float(z)
        if weighted:
            bins = 10**npz['arr_0'] #lengths the same
        else:
            bins = npz['arr_0'] #length 201
            bins = (bins[1:] + bins[:-1])/2.0 
        counts = npz['arr_1']
        p = interpolate.interp1d(bins,counts,fill_value=0,kind='linear',bounds_error=False)
        pdfs.append(p)


    '''
    --------------------------------------------------
    Make 2D DM,z histogram and relevant internal functions
    --------------------------------------------------
    '''
    dms = np.arange(0,15000,1.0)
    histogram2d = np.zeros((len(zs),len(dms))) #f(DM|z)

    for i,z in enumerate(zs):
        p = pdfs[i]
        histogram2d[i] = p(dms)
    histogram2d /= np.max(histogram2d)
    spline = interpolate.RectBivariateSpline(zs,dms,histogram2d) # use RectBivariateSpline instead of interp2d for speed over grid

    def get_fDMz(z):
        retval = spline(z,dms)[0]
        return normalize_area(retval)

    def bayes_single(zvals,DM): #f(z|DM) = f(DM|z)f(z)/f(DM)
        retval = np.zeros_like(zvals)

        if isinstance(DM,PDF): #This is a PDF
            for i,z in enumerate(zvals):
                x,xm,xp = getbounds(DM)
                retval[i] = spline(z,x)[0] / DM.calc(x)
        else:
            for i,z in enumerate(zvals):
                retval[i] = spline(z,DM)[0]

        return normalize_area(retval)

    
    '''
    --------------------------------------------------
    Determine Milky Way contribution to the DM
    --------------------------------------------------
    '''
    if isinstance(mwarg,(tuple,list)):
        gl, gb = mwarg
        if ymw:
            mwdm = pyymw16.dist_to_dm(gl, gb, 50e3)[0].value
        else:
            mwdm = NE2001.NE2001(gl,gb,50.0,ndir=-1,DIR=NEDIR)['output']['DM']
    else:
        mwdm = mwarg
    xs = np.arange(0,200,0.01)
    mw_pdf = PDF(xs,gaussian(xs,mwdm,0.2*mwdm)) #20% error

    '''
    --------------------------------------------------
    Determine CGM contribution to the DM
    --------------------------------------------------
    '''
    xs = np.arange(0,200,0.01) #same as above
    cgm_pdf = PDF(xs,gaussian(xs,65.0,15.0)) #50-80 pc cm^-3 range

    '''
    --------------------------------------------------
    Determine IGM contribution to the DM
    --------------------------------------------------
    '''
    xs = np.arange(dm-6*dmerr,dm+6*dmerr,0.01)
    obs_pdf = PDF(xs,gaussian(xs,dm,dmerr))


    
    dm_pdf = obs_pdf - mw_pdf - cgm_pdf - hostdm
    zpdf = PDF(zs,bayes_single(zs,dm_pdf))

    if evaluate:
        z,zm,zp = getbounds(zpdf) 
        return z,zm,zp
    else:
        return zpdf



    



if __name__ == "__main__":
    # FRB 121102 test
    #print(calcz(557.0,2.0,(174.95,-0.225138),dmhost=300.0,weighted=True,NEDIR='/home/michael/Research/Programs/NE2001/bin.NE2001/'))

    parser = argparse.ArgumentParser(description="FRB Redshift Estimator")

    parser.add_argument('--NE2001',dest='NEDIR',default='NE2001/bin.NE2001/',help="Path pointing to the NE2001 bin.NE2001/ directory location")
    parser.add_argument('--unweighted',dest="unweighted",action="store_true",default=False,help="Use uniform weighted distribution (versus matter weighted distribution")
    #Additional flag for YMW16 model
    parser.add_argument('--ymw', dest='ymw',action='store_true',default=False,help="Use YMW model instead of NE2001")
    parser.add_argument('--hostdm',type=float,default=0.0,help="Host DM [pc cm^-3]")#,nargs=1,dest="hostdm",action="store_const",default=0)
    parser.add_argument('--mwdm',type=float,default=None,help="Milky Way DM [pc cm^-3]")
    parser.add_argument('dm', action="store",type=float,help="Observed DM [pc cm^-3]")
    parser.add_argument('dmerr', action="store",type=float,help="Error on observed DM [pc cm^-3]")
    parser.add_argument('galcoord',action="store",type=float,nargs=argparse.REMAINDER,help="If --mwdm is not provided, two values separated by a space: Galactic latitude and Galactic longitude [deg]")
    
    results = parser.parse_args()
    weighted = not results.unweighted
    ymw = results.ymw
    if results.mwdm is not None: # Do not use NE2001
        print("z=%0.3f-%0.3f+%0.3f"%(calcz(results.dm,results.dmerr,results.mwdm,hostdm=results.hostdm,weighted=weighted,NEDIR=results.NEDIR, ymw = ymw)))
    else:
        gb, gl = results.galcoord
        print("z=%0.3f-%0.3f+%0.3f"%(calcz(results.dm,results.dmerr,(gb,gl),hostdm=results.hostdm,weighted=weighted,NEDIR=results.NEDIR, ymw = ymw)))


