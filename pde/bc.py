#TODO: These functions need docstrings
from numpy import linalg
class hom(object):
    """Collection of static (or class) methods which put particular initial
        values along with homogeneous boundary conditions on grids given as
        the argument
        """
    @staticmethod
    def np_safe_abs(x):
        import numpy
        try:
            xn=linalg.norm(x,2)
        except ValueError:
            xn=abs(x)
        return xn
        
    @staticmethod
    def _tpl(grid, iv,ds='test'):
        grid.l.info('bc.hom: Homogeneous initial data template _tpl called.')
        #TODO: Research numpy routines which may automate this
        tax=grid._get_timelike_axis()
        def tmp(x):
            """{}""".format(ds)
            if x[tax]>0:
                return 0.0
            if x[tax]==0:
                return iv(x)
        return tmp
    
    @classmethod 
    def xpowalpha(cls,grid,alpha,norm=False,sum=False,zerocont=False,cutoff=False,cutoffval=0):
        #TODO: Should be trivial to combine this and dd_xpowalpha. I should not be lazy
        grid.l.info('bc.hom: Setting initial data to some power of x.')
        grid.l.debug('bc.hom: Parameters to xpowalpha: alpha={},norm={},sum={},zerocont={},cutoff={},cutoffval={}'.format(alpha,norm,sum,zerocont,cutoff,cutoffval))
        alpha=float(alpha)
        try:
            cv=map(float,cutoffval)
        except TypeError:
            cv=[float(cutoffval) for i in range(0,grid.u.ndim)]
        
        if norm:
            if cutoff:
                def tmp(x):
                    return cls.np_safe_abs([x[i]*float(cutoff and (x[i]>=cv[i])) for i in range(0,grid.u.ndim)])**alpha
            else:
                def tmp(x):
                    return cls.np_safe_abs(x)**alpha
            return cls._tpl(grid, tmp)
        if sum:
            if zerocont:
                def tmp(x):
                    return max(sum(float(cutoff and (x[i]>=cv[i]))*x[i]**alpha for i in range(0,grid.u.ndim)),0.0)
            else:
                def tmp(x):
                    return sum(float(cutoff and (x[i]>=cv[i]))*x[i]**alpha for i in range(0,grid.u.ndim))
        return cls._tpl(grid, tmp)
    
    @classmethod
    def dd_xpowalpha(cls,grid,alpha,cutoff=False,poscut=True,cv=0,negate=False):
        grid.l.info('bc.hom: Setting initial data to some power of x.')
        grid.l.debug('bc.hom: Parameters to dd_xpowalpha: alpha={},cutoff={},cv={}'.format(alpha,cutoff,cv))
        c=cutoff
        if alpha is 0:
            def tmp(x): return float(x[1]<=cv)
            return cls._tpl(grid, tmp) 
            
        if alpha<1 and (1/alpha % 2==0):
            grid.l.warning('bc.hom: Forcing cutoff value at least zero since alpha<1 is an even fractional power.')
            cv=min(cv,0)
            c=True
            if negate:
                poscut=True
            else:
                poscut=False

        if c:
            if negate:
                if poscut:
                    def tmp(x):
                        return sum(-1*pow(float(x[i]>=cv)*x[i],alpha) for i in range(0,grid.u.ndim))
                else:
                    def tmp(x):
                        return sum(-1*pow(-1*float(x[i]<cv)*x[i],alpha) for i in range(0,grid.u.ndim))
            else:
                if poscut:
                    def tmp(x):
                        return sum(pow(float(x[i]>=cv)*x[i],alpha) for i in range(0,grid.u.ndim))
                else:
                    def tmp(x):
                        return sum(pow(-1*float(x[i]<cv)*x[i],alpha) for i in range(0,grid.u.ndim))
        else:
            if negate:
                def tmp(x):
                    return sum(pow(-1*float(x[i]<0)*x[i],alpha)-pow(float(x[i]>=0)*x[i],alpha) for i in range(0,grid.u.ndim))
            else:
                def tmp(x):
                    return sum(pow(float(x[i]>=0)*x[i],alpha)-pow(-1*float(x[i]<0)*x[i],alpha) for i in range(0,grid.u.ndim))
        return cls._tpl(grid, tmp)
    
    @classmethod
    def test1(cls,grid):    
        """Homogeneous boundary conditions and prescribed initial conditions."""
        return cls.xpowalpha(grid,2,norm=True)
    
    @classmethod
    def barenblatt(cls,grid,t,sigma):
        grid.l.info('bc.hom: Setting initial data to a Barenblatt point source solution.')
        grid.l.debug('bc.hom: Parameters to barenblatt: sigma={},t={}'.format(sigma,t))
        def tmp(x):
            return cls._barenblatt(t, x, grid.u.ndim-1, sigma)
        return cls._tpl(grid, tmp)
        
    @classmethod
    def _barenblatt(cls,t,x,ndim,sigma,n0=1):
        """Point-source solution to the PDE
                ut-lapl(u^m)=0
            or
                ut-div(u^sigma*grad(u))=0
            where sigma=m-1.  LaTeX statement of solution:
             t^{\frac{-N}{N\sigma+2}}
             \left[ \frac{\sigma}{N\sigma+2}
              \left( \eta_{0}^{2}-\frac{|x|^{2}}{t^{\frac{2}{N\sigma+2}}}
              \right)_{+}
              \right]^{\frac{1}{\sigma}}
        """
        #TODO: Vectorize this code
        ndim,sigma=float(ndim),float(sigma)
        xn=cls.np_safe_abs(x)
        return pow(t,-ndim/(ndim*sigma+2.0))*pow(
                   (
                    sigma/(ndim*sigma+2.0)
                    )*max(
                                             (pow(n0,2.0)-(pow(t,-2.0/(ndim*sigma+2.0))*pow(xn,2.0)))
                                             ,0.0)
                   ,1/sigma)
                   
                
        
class inhom(object):
    pass