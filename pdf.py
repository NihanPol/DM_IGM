import numpy as np
import scipy.interpolate as interpolate
import scipy.special as special
import matplotlib.pyplot as plt
import copy


EPS = special.erf(1.0/np.sqrt(2))/2.0




def sigfig(x,xm,xp,eps=0.00000001): # eps fixes rounding errors
    argm = -int(np.floor(np.log10(np.abs(xm)))) #argument to np.around()
    argp = -int(np.floor(np.log10(np.abs(xp))))
    arg = np.max([argm,argp])

    rx = np.around(x+eps, arg)
    rxm = np.around(xm+eps, arg)
    rxp = np.around(xp+eps, arg)
    #else: #Does not properly handle large integers!
    #    rx = int(x+eps)
    #    rxm = int(xm+eps)
    #    rxp = int(xp+eps) 
    if rx == 0.0:
        rx = np.abs(rx) #force +0.0

    try:
        rx = format(rx,'0.%if'%arg)
        rxm = format(rxm,'0.%if'%arg)
        rxp = format(rxp,'0.%if'%arg)
    except ValueError:
        print "Error:",rx,arg

    return rx, rxm, rxp

def pm_formatter(x,xm,xp):
    #xe = np.max([xm,xp])
    rx, rxm, rxp = sigfig(x,xm,xp)
    rx, rxm, rxp = sigfig(float(rx),float(rxm),float(rxp)) #lazy way to fix +/-0.009 going to +/-0.010, thus will return +/-0.01
    #print rx,rxm,rxp
    if rxm == rxp:
        return "$%s \\pm %s$"%(rx,rxm)
    else:
        return "$%s_{-%s}^{+%s}$"%(rx,rxm,rxp)


class PDF:
    def __init__(self,x,y,kind='linear',truncate=True):
        self.x = np.array(x,dtype=np.float)
        self.y = np.array(y,dtype=np.float) #ensure histograms are appropriately converted
        self.kind = kind
        if truncate:
            self.truncate() #free run
        else:
            self.run()

    def run(self):
        # Sort because otherwise interpolate may fail
        inds = np.argsort(self.x)
        self.x = self.x[inds]
        self.y = self.y[inds]

        #print self.y,"before"
        # Normalize to unit area
        self.y /= np.trapz(self.y,x=self.x)
        #print self.y,"after"
        self.f = interpolate.interp1d(self.x,self.y,kind=self.kind,fill_value=0.0,bounds_error=False)
        self.c = self.cdf()
        

    def calc(self,x):
        #if self.x[0] <= x <= self.x[-1]:
        if isinstance(x,(int,long,float)):
            if self.f.x[0] <= x <= self.f.x[-1]:
                return self.f(x).item() #returns zero-dimension array
            return 0.0
        if isinstance(x,(list,np.ndarray)):
            y = np.zeros_like(x)
            #inds = np.where(x < self.x[0])[0]
            #y[inds] = 0.0#self.y[0]
            #inds = np.where(x > self.x[-1])[0]
            #y[inds] = 0.0#self.y[-1]
            inds = np.where(np.logical_and(self.x[0]<=x,x<=self.x[-1]))[0] # no longer need this due to the fill_value=0?
            #print "inds",inds,x,self.x[0],self.x[-1]
            y[inds] = self.f(x[inds])
            return y
        #return None
    __call__ = calc


    def get_ranges(self,other,func=None,newdx=None,nsample=10000):#None):
        # For now assume x ranges are the same
        #maxrange = np.max(np.abs(self.x))+np.max(np.abs(other.x))
        if func is None:#????
            #maxrange = np.ptp(self.x) + np.ptp(other.x)
            maxrange = np.max(np.abs(self.x)) + np.max(np.abs(other.x))
        else:
            #maxrange = func(np.ptp(self.x),np.ptp(other.x)) 
            maxrange = func(np.max(self.x),np.max(other.x))  #maybe too much?
        if newdx is not None:
            dx = newdx
        else:
            dx = self.x[1]-self.x[0]
            otherdx = other.x[1]-other.x[0]
            dx = min([dx,otherdx]) #get the minimum of the two
        if nsample:
            xnew = np.linspace(-2*maxrange,2*maxrange+dx,nsample)
        else:
        #print "maxrange",maxrange
            xnew = np.arange(-2*maxrange,2*maxrange+dx,dx) #how to handle this?
        #print len(xnew)
        #xnew = self.x
        ynew = np.zeros_like(xnew)
        return xnew,ynew


    def sum(self,other,truncate=True):
        if isinstance(other,PDF):
            xnew,ynew = self.get_ranges(other)
            pre1 = self.calc(self.x) #Pre-define variables
            #pre1 = self.calc(xnew) #Pre-define variables
            # other.calc(xnew[i]-self.x) recalculates a lot of the same values for just a shift by the same amount dx
            dx = xnew[1] - xnew[0]
            length = len(self.x)
            #length = len(xnew)
            diffx = np.arange(xnew[0],xnew[-1]+length*dx,dx) - self.x[-1]
            pre2 = other.calc(diffx)
            #plt.plot(self.x,self.y)
            #plt.plot(other.x,other.y)
            #plt.plot(diffx,pre2) 



            #print other.x,diffx
            #print pre1,pre2,diffx

            for i in range(len(xnew)):
                #if i % 10000 == 0 and i != 0:
                #    print i,len(xnew)
                ynew[i] = np.trapz(pre1*other.calc(xnew[i]-self.x),x=self.x)
                #ynew[i] = np.trapz(pre1*pre2[i:length+i][::-1],x=self.x) #This is now failing, why?
            #plt.plot(xnew,ynew)
            #plt.show()
            pdf = PDF(xnew,ynew,kind=self.kind)
            if truncate:
                pdf.truncate()
            return pdf
        elif isinstance(other,(int,long,float)):
            newself = copy.deepcopy(self)
            newself.x += other
            if truncate:
                newself.truncate()
            return newself
    __add__ = sum
    __radd__ = sum
    __iadd__ = sum

    def minus(self,other,truncate=True):
        return self.sum(-other,truncate=truncate)
    __sub__ = minus
    def __rsub__(self,other,truncate=True):
        return(-self).sum(other)
    __isub__ = minus



    def product(self,other,truncate=True):
        if isinstance(other,PDF):
            xnew,ynew = self.get_ranges(other,func=lambda x,y:x*y)
            pre1 = self.calc(self.x) #pre-defined
            pre2 = np.abs(self.x)
            for i in range(len(xnew)):
                ynew[i] = np.trapz(pre1*other.calc(xnew[i]/self.x) / pre2,x=self.x)
            pdf = PDF(xnew,ynew,kind=self.kind)
            if truncate:
                pdf.truncate()
            return pdf
        elif isinstance(other,(int,long,float)):
            newself = copy.deepcopy(self)
            newself.x *= other
            if truncate:
                newself.truncate()
            return newself
        elif isinstance(other,(np.ndarray,list)):
            if len(other) != len(self.x):
                raise ValueError("operands could not be broadcast together with shapes (%i,) (%i,)"%(len(self.x),len(other)))
            newself = copy.deepcopy(self)
            newself.y *= np.array(other) #vastly different operation from a single value!
            if truncate:
                newself.truncate()
            return newself
    __mul__ = product
    __rmul__ = product
    __imul__ = product


    def inverse(self,other,truncate=True):
        pass
    
    def __neg__(self,truncate=True):
        return self.product(-1,truncate=truncate)
    __pos__ = lambda self: self
    

    def square(self,truncate=True):
        dx = self.x[1]-self.x[0]
        if np.max(self.x) <= 2:# and dx >= 0.05:
            dx = 0.01*dx
        xnew,ynew = self.get_ranges(self,func=lambda x,y:x*y,newdx=dx)
        
        inds = np.where(xnew>0)[0]
        sqrt = np.sqrt(xnew[inds])
        #print "sq",xnew,self.x,sqrt
        ynew[inds] = 1.0/(2*sqrt) * (self.calc(sqrt) + self.calc(-1*sqrt))

        #raise SystemExit
        '''
        for i in range(len(xnew)):
            if xnew[i] <= 0:
                continue
            sqrt = np.sqrt(xnew[i])
            ynew[i] = 1.0/(2*sqrt) * (self.calc(sqrt) + self.calc(-1*sqrt))
        '''
        pdf = PDF(xnew,ynew,kind=self.kind)
        if truncate:
            pdf.truncate()
        return pdf

        #return self.product(self,truncate=truncate)

    def __pow__(self,value,truncate=True):
        if value == 2:
            return self.square()
        pass

    def sqrt(self,truncate=True):
        xnew,ynew = self.get_ranges(self,func=lambda x,y:np.sqrt(x*y))
        #print xnew,ynew
        inds = np.where(xnew>0)[0]
        ynew[inds] = 2*np.abs(xnew[inds])*self.calc(xnew[inds]**2)

        '''
        for i in range(len(xnew)):
            if xnew[i] <= 0:
                continue
            ynew[i] = 2*np.abs(xnew[i])*self.calc(xnew[i]**2)
        '''
        pdf = PDF(xnew,ynew,kind=self.kind)
        if truncate:
            pdf.truncate()
        return pdf
        

    def cdf(self):
        # From utilities, pdf_to_cdf
        dx = self.x[1]-self.x[0]
        return np.cumsum(self.y)*dx

    def mean(self):
        return np.trapz(self.x*self.y,x=self.x)
    def median(self):
        ind = np.argmin(np.abs(self.c-0.50))
        return self.x[ind]
    def mode(self):
        return self.x[np.argmax(self.y)]
    def likelihood_evaluator(self,median=False,pm=True,values=None):
        # From utilities, likelihood_evaluator
        """
        median: if True, use the median value, otherwise the peak of the pdf
        pm: xminus and xplus are the plus/minus range, not the actual values
        Future: give it values to grab off the CDF (e.g. 2 sigma, 99%, etc)
        values: use this array
        """
        x = self.x
        y = self.y
        ycdf = self.c

        if values is None:
            if median:
                yb = 0.50   #Now take the median!
            else:
                indb = np.argmax(y)
                yb = ycdf[indb]

            ya = yb - EPS
            yc = yb + EPS
            yd = 0.95

            inda = np.argmin(np.abs(ycdf - ya))
            if median:
                indb = np.argmin(np.abs(ycdf - yb)) 
            indc = np.argmin(np.abs(ycdf - yc))
            indd = np.argmin(np.abs(ycdf - yd))

            inds = np.arange(inda,indc+1) #including indc    
            #print indc-inda,np.trapz(L[inds],x=Vrs[inds])
            xval = x[indb]
            if pm:
                xminus = x[indb] - x[inda]
                xplus = x[indc] - x[indb]
            else:
                xminus = x[inda]
                xplus = x[indc]
            x95 = x[indd]

            return xval,xminus,xplus,x95
        else:
            retval = np.zeros_like(values)
            for i,v in enumerate(values):
                indv = np.argmin(np.abs(ycdf - v))
                retval[i] = x[indv]
            return retval

    def __len__(self):
        return len(self.x)


    def truncate(self,eps=1e-10):#,buffer=5): #1e-15
        if np.all(self.y==0):
            return
        i = 0
        while True:
            if i == 30:
                raise SystemExit
            i += 1
            #print "eps",eps,self.y
            if np.all(self.y<eps):
                eps /= 10.0
                continue
            inds = np.where(self.y>eps)[0]
            if len(inds) == 0:
                eps *= 10.0
                #print "eps",eps
            #elif len(inds) == len(self.y):
            else:
                break
        self.x = self.x[inds]
        self.y = self.y[inds]

        self.run()
        #def xtruncate(self,statement):
        #if np.all(self.y==0):
        #    return
        #inds = np.where[statement

    def resample(self,n):
        newx = np.linspace(self.x[0],self.x[-1],n)
        newy = self.calc(newx)
        self.x = newx
        self.y = newy
        self.run()




    def plot(self):
        plt.plot(self.x,self.y,'k')
        plt.show()

    def save(self,filename):
        np.savez(filename,x=self.x,y=self.y)
