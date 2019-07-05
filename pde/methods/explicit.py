"""Explicit methods to calculate PDE solutions on grids"""
from .base import method, differencer

class diffusion(method):
    """
    Calculate diffusion operator corresponding to the PDE
    u_t - (ktxu(t, x, u) u_x)_x - b*(u^beta)_x
    If `ktxu` is not given, it is set according to u^(alpha).
    """
    def __init__(self,grid,b=0,beta=1,alpha=1,ktxu=None):
        super(diffusion,self).__init__(grid)
        if ktxu is None:
            def tmp(t,x,u):
                if u <= 0:
                    return 0
                return pow(u, alpha)
            self.ktxu = tmp
        else:
            self.ktxu=ktxu
        self.b=b
        self.beta=beta
        
    def calculate(self):
        return self._calculate_python_1d(self.g, self.b, self.beta, self.ktxu)
    
    @staticmethod
    def _calculate_python_1d(grid, b, beta, ktxu):
        tax = grid._get_timelike_axis()
        oax=[l for l in range(grid.u.ndim) if l != tax]
        oax=oax[0] #FIXME: This supports only one space dimension.
        #FIXME: This implementation seems to have been subtly broken due to excessive (?)
        # but not careful encapsulation of details of grid. For now, we'll completely 
        # break this abstraction to see how this works...
        tvals = grid.x[tax]
        dt = tvals[1] - tvals[0]
        
        # This is "safe" since we didn't handle irregular space grids anyways
        xvals = grid.x[oax]
        dx = xvals[1] - xvals[0]
        
        cfl = dt/pow(dx, 2)
        for (i, ti) in enumerate(tvals):
            if i == 0:
                continue
            for (j, xj) in enumerate(xvals):
                if j == 0 or j == (len(xvals)-1):
                    continue
                yi=grid.u[i - 1, j]
                yip=grid.u[i - 1, j + 1]
                yim=grid.u[i - 1, j - 1]
                
                kyi = ktxu(ti,xj,yi)
                
                kyim = ktxu(ti,xj-dx,yim)
                
                kyip = ktxu(ti,xj+dx,yip)
                
                # Scheme 1 from Samarskii
                aiv = (kyim + kyi)/2
                aipv = (kyip + kyi)/2
                
                if yi<=0: # Outside support of solution
                    powyi = 0
                else: 
                    powyi = pow(yi, beta - 1)
                
                if yim <= 0:
                    powyim = 0
                else:
                    powyim = pow(yim, beta - 1)
                
                grid.u[i, j]=yi + (
                    cfl*(
                    aipv*(yip - yi) - aiv*(yi - yim)
                    ) + ((dt/dx) * b * (powyi-powyim))
                    )
        return grid

class ftcs_1d(method):
    #TODO: Implement FTCS for objective functions of form F(t,x,u,ux1,...,uxn,ux1x1,...,uxnxn)
    # We should be able to write a general method for equations of the form
    #   u_t=F(t,x,u,u_x,u_xx)
    # FTCS approximations are unstable for convection-type equations unless
    # you're _very_ careful with the scale of the problem.
    def __init__(self, grid, rhs):
        super(ftcs_1d,self).__init__(grid)
        import inspect
        if not inspect.isfunction(rhs):
            raise ValueError('rhs is not a function.')
        a=inspect.getargspec(rhs)
        self.l.debug('Argument specification retrieved using equivalent of ffn.__code__.co_varnames[:ffn.__code__.co_argcount]')
        self.l.debug('Argument specification of rhs: {}'.format(a))
        ap=set(['x','t','u','ux','uxx'])
        self.ap=ap
        self.l.debug('Checking if rhs takes any arguments from {}'.format(ap))
        if len(set(a[0]).intersection(ap))==0:
            raise ValueError('rhs must take at least one of x,t,u,ux,uxx as an argument.')
        self.passargs=set(a[0]).intersection(ap)
        self.l.info('Given rhs function takes arguments {}'.format(self.passargs))
        self.l.debug('Docstring for rhs: {}'.format(rhs.__doc__))
        self.ffn=rhs

    def _calc_pt(self,ind,tax=None,oax=None):
        if tax is None:
            tax=self.g._get_timelike_axis()

        if oax is None:
            oax=[l for l in range(self.g.N) if l != tax]
            oax=oax[0]
        pta=self.g._ax_step(ind,tax,-1)
        dt=self.g._dist(self.g.ind_to_x(ind),self.g.ind_to_x(pta))
        tmp=dict()
        if 'x' in self.passargs:
            pt=self.g.ind_to_x(pta)
            tmp['x']=pt[oax]

        if 't' in self.passargs:
            try:
                tmp['t']=pt[tax]
            except IndexError:
                pt=self.g.ind_to_x(pta)
                tmp['t']=pt[tax]
            
        if 'u' in self.passargs:
            tmp['u']=self.g.u[tuple(pta)]

        if 'ux' in self.passargs:
            tmp['ux']=differencer.dx(self.g, oax, pta, type=1)
        
        if 'uxx' in self.passargs:
            tmp['uxx']=differencer.uxx_centered(self.g, oax, pta)
        
        self.g.u[tuple(ind)]=dt*self.ffn(**tmp)+self.g.u[tuple(pta)]
        
    def calculate(self):
        tax=self.g._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]
        for i in self.g.get_ind_range(ax=tax,aset='active'):
            curr=self.g.get_plane_ind(val=i,ax=tax,aset='active')
            for ind in curr:
                self._calc_pt(ind,tax,oax)