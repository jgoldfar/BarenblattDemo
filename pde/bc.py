#TODO: These functions need docstrings
from numpy import linalg

class hom(object):
    """
    Collection of static (or class) methods which put particular initial
    values along with homogeneous boundary conditions on grids given as
    the argument"""
    
    @staticmethod
    def np_safe_abs(x):
        """
        Safe abs function that attempts to calculate the 2-norm of the argument
        but falls back to the regular old abs function.
        """
        try:
            return linalg.norm(x,2)
        except ValueError:
            return abs(x)
        
    @staticmethod
    def _tpl(grid, iv, ds='test'):
        """
        Assign the initial value function `iv` to the grid `grid`.
        """
        #TODO: Research numpy routines which may automate this
        tax=grid._get_timelike_axis()
        def tmp(x):
            if x[tax]>0:
                return 0.0
            if x[tax]==0:
                return iv(x)
        return tmp

    @classmethod
    def dd_xpowalpha(cls,grid,alpha,cutoff=False):
        """
        Set initial data to (-x)^alpha. If `cutoff` is true, continue the function for
        positive values by zero before raising to the power `alpha`.
        """
        grid.l.info('bc.hom: Setting initial data to (-x)^alpha.')
        grid.l.debug('bc.hom: Parameters to dd_xpowalpha: alpha={},cutoff={}'.format(alpha,cutoff))
        if alpha is 0:
            def tmp(x): return float(x[1]<=0)
            return cls._tpl(grid, tmp) 

        if cutoff:
            def tmp(x):
                return sum(pow(-1*float(x[i]<0)*x[i],alpha) for i in range(0,grid.u.ndim))
        else:
            def tmp(x):
                return sum(pow(float(x[i]>=0)*x[i],alpha)-pow(-1*float(x[i]<0)*x[i],alpha) for i in range(0,grid.u.ndim))
        return cls._tpl(grid, tmp)
    
    @classmethod
    def barenblatt(cls,grid,t,sigma):
        """
        Set initial data to the Barentblatt point-source solution with prescribed `sigma`
        and at the moment `t`. Subsequent evolution should follow the Barenblatt profile.
        """
        grid.l.info('bc.hom: Setting initial data to a Barenblatt point source solution.')
        grid.l.debug('bc.hom: Parameters to barenblatt: sigma={},t={}'.format(sigma,t))
        def tmp(x):
            return cls._barenblatt(t, x, grid.u.ndim-1, sigma)
        return cls._tpl(grid, tmp)
        
    @classmethod
    def _barenblatt(cls,t,x,ndim,sigma,n0=1):
        """
        Point-source solution to the PDE
            ut-lapl(u^m)=0
        or
            ut-div(u^sigma*grad(u))=0
        where sigma=m-1. Calculated according to
        \[
         t^{\frac{-N}{N\sigma+2}}
         \left[ \frac{\sigma}{N\sigma+2}
          \left( \eta_{0}^{2}-\frac{|x|^{2}}{t^{\frac{2}{N\sigma+2}}}
          \right)_{+}
          \right]^{\frac{1}{\sigma}}
        \]
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
