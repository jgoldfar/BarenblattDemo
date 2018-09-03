"""Implements implicit PDE solution methods."""
from base import method
from numpy import eye, diagflat, linalg, amax, fabs,concatenate, dot,nan_to_num
#TODO: Add logging to these methods

class diffusion(method):
    """Attempts to solve problems of the form
        u_t=(k(u,x,t) u_x)_x + f(ux,u,x,t)
        Of course, convergence of the numerical scheme is not guaranteed!
    """
    def __init__(self,grid=None,kfn=None,ffn=None):#,b=None,beta=None,sigma=None):
        super(diffusion,self).__init__(grid)
#        ns=False
#        if kfn is None:
#            if sigma is None:
#                sigma=1.0
#                ns=True
#            else:
#                sigma=float(sigma)
#            print "Creating appropriate diffusion coefficient function."
#            def kfn(u=0,x=0,t=0):
#                return sigma*pow(u,sigma-1)
#            print """
#            def kfn(u,x,t):
#                return {}*pow(u,{}-1)
#            """.format(sigma,sigma)
#        if ffn is None:
#            if beta is None:
#                beta=1.0
#                ns=True
#            else:
#                beta=float(beta)
#            if b is None:
#                b=0.0
#                ns=True
#            else:
#                b=float(b)
#            print "Creating appropriate convection coefficient function."
#            def ffn(ux=0,u=0,x=0,t=0):
#                return b*beta*pow(u,beta-1.0)*ux
#            print """
#            def ffn(ux,u,x,t):
#                return {}*{}*pow(u,{}-1.0)*ux
#            """.format(b,beta,beta)
        self.kf=kfn
        self.ff=ffn
        #self.g.l.info("defaults={}, b={}, beta={}, sigma={}.".format(ns,b,beta,sigma))

    def __call__(self,scheme=1,**kwargs):
        self.calculate(scheme=scheme,**kwargs)
        
    def calculate(self,scheme=1,**kwargs):
        if self.g.u.ndim==2:
            if scheme is 3:
                c=self._calculate_1d_scheme_matus
            elif scheme is 2:
                c=self._calculate_1d_scheme_b
            else:
                c=self._calculate_1d_scheme_a
            c(self.g,self.kf,self.ff,**kwargs)
        else:
            raise NotImplementedError('Grids with more than one spatial dimension not supported.')

    @staticmethod
    def _make_trilinear_safe(subdiag,diag,supdiag):
        return nan_to_num(diagflat(subdiag,-1)+diagflat(supdiag,1)+diagflat(diag))
    
    @staticmethod
    def _find_diff(yv,xv,i,diff):
        #print diff
        if diff is 2:
            # Forward
            return (yv[i+1]-yv[i])/(xv[i+1]-xv[i])
        elif diff is 1:
            # Centered
            return (yv[i+1]-yv[i-1])/(xv[i+1]-xv[i-1])
        else:
            # Backward
            return (yv[i]-yv[i-1])/(xv[i]-xv[i-1])

    @staticmethod
    def _solve_sys_cat(A,b,lb,rb,allow_singular=1,logger=None):
        """Solves the system Ax=b for x, resorting to least-squares in the bad case.  Also
            pads the solution with the left and right endpoint values lb and rb."""
        try:
            t=concatenate([[lb],linalg.solve(A,b),[rb]])
            #print linalg.solve(A,b)
        except linalg.LinAlgError:
            print "Caught Singular system"
            if allow_singular is 'log':
                try:
                    print "Singular system; returning least-squares solution."
                    t=concatenate([[lb],dot(linalg.pinv(A),b),[rb]])
                except linalg.LinAlgError:
                    print "Error finding least-squares solution.  Giving zero."
                    print A
                    print b
                    return b*0
            else:
                raise
        return t           
  
    @staticmethod
    def _ai(kfn,xv,tv,vv,i,j,ai_param=1):
        #print kfn
        if ai_param is 3:
            return 2.0*(kfn(u=vv[i-1])*kfn(u=vv[i],x=xv[i],t=tv[j]))/(kfn(u=vv[i-1])+kfn(u=vv[i]))
        elif ai_param is 2:
            vav=(vv[i-1]+vv[i])/2.0
            #xav=(xv[i-1]+xv[i])/2.0
            return kfn(u=vav)
        else:
            return 0.5*(kfn(u=vv[i-1])+kfn(u=vv[i]))
    
    @classmethod
    def _calculate_1d_scheme_a(cls,grid,kfn,ffn,**kwargs):
        """Implements a modified version of implicit scheme a from Samarskii (pg 520.)"""
        
        tax=grid._get_timelike_axis()
        if tax is 0:
            oax=1
        else:
            oax = 0
        t=grid.x[tax]
        x=grid.x[oax]
#        
        
        # naive can be 0: Changes to RHS only, 1: Changes to matrix only. Setting naive to 1 
        # means that the derivative term will depend on the current time, while setting to
        # 0 puts the dependence of the spatial derivative term in the previous time.
        naive=kwargs.get('naive',1)
        #print naive
        # diff is the parameterization for derivative. 0: Backward difference, 1: Centered, 2: Forward
        diff=kwargs.get('diff',1)
        #print diff
        # allow_singular is a flag that tells how you would like to handle a singular LHS matrix.
        #    since this seems to be a common issue, the default behavior is to write a message to the
        #    log and use the moore-penrose pseudoinverse to find a least-squares solution to the
        #    same problem.  To set this behavior, choose 'log' or do nothing; to re-raise the error,
        #    set it to 'raise'.
        allow_singular=kwargs.get('allow_singular','log')
        #print allow_singular
        # ai_param tells which parameterization is chosen for the conduction coefficient in the 
        #    numerical scheme.
        ai_param=kwargs.get('ai_param',1)
        #TODO: Allow parameterization of the value passed to convection coefficient
        # gfn_param allows you to tweak the values being passed to ffn, which represents
        #    the convection coefficient.  Currently not implemented, but it makes sense
        #    to implement the same strategies as for ai (and probably use these here as well)
        #    but it may be worth it to experiment with others?
        if 'gfn_param' in kwargs:
            grid.l.warning('gfn_param currently not an implemented option.')
        gfn_param=kwargs.get('gfn_param',1)
        
        # slow_diff reduces the diffusion coefficient (making diffusion slower!)
        #    Larger values correspond to slower diffusion
        slow_diff=kwargs.get('slow_diff',10)
        
        # slow_conv reduces the convection coefficient (making convection slower)
        #    Currently not implemented.  Larger values correspond to slower convection.
        if 'slow_conv' in kwargs:
            grid.l.warning('slow_conv currently not an implemented option.')
        slow_conv=kwargs.get('slow_conv',1)
        #print ai_param
        def cl(kfn,ffn,x,t,y,i,j):
            #return i
            tau=t[j]-t[j-1]
            #print t[j]
            xf=x[i+1]-x[i]
            xb=x[i]-x[i-1]
            xav=0.5*(xf+xb)
            if naive is 1:
                return (-tau*cls._ai(kfn,x,t, y, i, j, ai_param))/(xav*xb*slow_diff)
            else:
                if diff is 2:   #Forward
                    return (-tau*cls._ai(kfn,x,t, y, i, j, ai_param))/(xav*xb*slow_diff)
                elif diff is 1: #Centered
                    return (-tau*cls._ai(kfn,x,t, y, i, j, ai_param))/(xav*xb*slow_diff)((tau*ffn(ux=1,u=y[i],t=t[j],x=x[i]))/(2*xb*slow_conv))
                else:           #Backward
                    return (-tau*cls._ai(kfn,x,t, y, i, j, ai_param))/(xav*xb*slow_diff)+((tau*ffn(ux=1,u=y[i],t=t[j],x=x[i]))/(xb*slow_conv))
        
        def c(kfn,ffn,x,t,y,i,j):
#            print "i:",i
            #return i+10
            tau=t[j]-t[j-1]
            xf=x[i+1]-x[i]
            xb=x[i]-x[i-1]
            xav=0.5*(xf+xb)
            
            if naive is 1:
                return 1+(tau/(xav*xf*slow_diff))*cls._ai(kfn,x,t, y, i+1, j, ai_param)+(tau/(xb*xav*slow_diff))*cls._ai(kfn, x,t,y, i, j, ai_param)
            else:
                if diff is 2:   #Forward
                    return 1+(tau/(xav*xf*slow_diff))*cls._ai(kfn,x,t, y, i+1, j, ai_param)+(tau/(xb*xav*slow_diff))*cls._ai(kfn, x,t,y, i, j, ai_param)+((tau(ffn(ux=1,u=y[i],x=x[i],t=t[i]))))
                elif diff is 1: #Centered
                    return 1+(tau/(xav*xf*slow_diff))*cls._ai(kfn,x,t, y, i+1, j, ai_param)+(tau/(xb*xav*slow_diff))*cls._ai(kfn, x,t,y, i, j, ai_param)
                else:           #Backward
                    return 1+(tau/(xav*xf*slow_diff))*cls._ai(kfn,x,t, y, i+1, j, ai_param)+(tau/(xb*xav*slow_diff))*cls._ai(kfn, x,t,y, i, j, ai_param)-((tau(ffn(ux=1,u=y[i],x=x[i],t=t[i]))))
        
        def cu(kfn,ffn,x,t,y,i,j):
            #return i+20
            tau=t[j]-t[j-1]
            xf=x[i+1]-x[i]
            xb=x[i]-x[i-1]
            xav=0.5*(xf+xb)
            #print xf," ",xb," ",xav
            if naive is 1:
                return (-tau*cls._ai(kfn,x,t, y, i+1, j, ai_param))/(xav*xf*slow_diff)
            else:
                if diff is 2:   #Forward
                    return ((-tau*cls._ai(kfn,x,t, y, i+1, j, ai_param))/(xav*xf*slow_diff))-((tau*ffn(ux=1,u=y[i],t=t[j],x=x[i]))/(xf*slow_conv))
                elif diff is 1: #Centered
                    return ((-tau*cls._ai(kfn,x,t, y, i+1, j, ai_param))/(xav*xf*slow_diff))-((tau*ffn(ux=1,u=y[i],t=t[j],x=x[i]))/(2*xf*slow_conv))
                else:           #Backward
                    return (-tau*cls._ai(kfn,x,t, y, i+1, j, ai_param))/(xav*xf*slow_diff)
        
        def b(kfn,ffn,x,t,y,i,j):
            #return float(i)
            tau=t[j]-t[j-1]
            if naive is 1:
                #print "b knows naive"
                yt=cls._find_diff(y, x, i, diff)
                tmp=tau*ffn(ux=yt,u=y[i],x=x[i],t=t[j])+y[i]
                if i is 1:
                    tmp=tmp-cl(kfn, ffn, x,t,y, 1, j)*y[0]
                elif i is n-1:
                    t=tmp-cu(kfn, ffn, x,t,y, n-1, j)*y[n]
                return tmp
            else:
                #print "b is naive"
                return y[i]
        g=grid.u
        n=len(g[0])

        for j in grid.get_ind_range(ax=tax,aset='active'):
            gl=g.copy()
            yi=gl[j-1]
            subdiag=[cl(kfn,ffn,x,t,yi,i,j) for i in xrange(2,n-1)]
            diag=[c(kfn,ffn,x,t,yi,i,j) for i in xrange(1,n-1)]
            supdiag=[cu(kfn,ffn,x,t,yi,i,j) for i in xrange(1,n-2)]
            bvec=nan_to_num([b(kfn,ffn,x,t,yi,i,j) for i in xrange(1,n-1)])
            #print p
            #print bvec
            A=cls._make_trilinear_safe(subdiag, diag, supdiag)
            g[j]=cls._solve_sys_cat(A, bvec, yi[0], yi[-1], allow_singular)
            #print grid.u
            #print "Solution finished for i=",i 
            
    @classmethod
    def _calculate_1d_scheme_b(cls,grid,kfn,ffn,**kwargs):
        """Implements a modified version of implicit scheme b from Samarskii (pg 520.)"""
        #TODO: Implement this in human-readable fashion
        
        naive=kwargs.get('naive',1) # Same as scheme a
        
        diff=kwargs.get('diff',1) # Same as scheme a
        
        # steperr is the error goal in iteration for each time step
        steperr=kwargs.get('steperr',10**-7)
        
        # nmax is the maximum number of iterations; usually converges rapidly (or diverges!) 
        nmax=kwargs.get('nmax',4)
        
        allow_singular=kwargs.get('allow_singular','log') # Same as scheme a
        
        ai_param=kwargs.get('ai_param',1) # Same as scheme a
        
        gfn_param=kwargs.get('gfn_param',1)
            
    @classmethod
    def _calculate_1d_scheme_matus(cls,grid,kfn,ffn,**kwargs):
        #Matus and Lemeschevsky's algorithm is an explicit algorithm, and is implemented in
        # one dimension in the explicit package; barring small numerical issues, either
        # the ftcs_1d or diffusion packages can cover this
        #TODO: Reference the algorithms in the explicit package 
        pass

class diffusion_old(method):
    #TODO: Verify the calculations involved in this code
    """Deprecated version of the implicit diffusion code""" 
    def __init__(self,grid=None,b=0,beta=1,sigma=1,ktxu=None):
        super(diffusion_old,self).__init__(grid)
        self.sigma=sigma
        self.b=b
        self.beta=beta
        
    def calculate(self,scheme='rev2'):
        if self.g.u.ndim==2:
            self._calculate_python_1d_rev2()
        else:
            raise NotImplementedError('Grids with more than one spatial dimension not supported.')
        
    def _calculate_python_1d_rev3(self,steperr=10**-7,itermax=3):
        """Implements scheme b in Samarskii's text with some slight modifications to give a 
            convection term rather than a reaction term.  Implements the convection as a 
            backward difference, and the convection terms appear on the RHS, not in the 
            coefficient matrix.
            
            This scheme is not rearranged or simplified by a computer algebra system. 
            """ 
            #TODO: This scheme can be optimized in the same way as rev2
        def aiv(l,v):
            #return 0.5*(self.ktxu(v[i-1])+self.kxtu(v[i])))
            return 0.5*(pow(v[l-1],self.sigma)+pow(v[l],self.sigma))
        def ci(h,tau,b,beta,sigma,y,l):
            return 1+((tau*aiv(l+1,y))/pow(h,2))+((tau*aiv(l,y))/pow(h,2))
        def cil(h,tau,b,beta,sigma,y,l):
            return -(tau*aiv(l,y))/pow(h,2)
        def ciu(h,tau,b,beta,sigma,y,l):
            return -tau*aiv(l+1,y)/pow(h,2)
            
        tax=self.g._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]#FIXME: There should be a better way to do this
        #FIXME: Assumes constants dt and dx: uniform grids only!
        g=self.g.u
        dt=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((1,0)))
        dx=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((0,1)))
        n=len(g[0])-2

        for i in self.g.get_ind_range(ax=tax,aset='active'):
            gl=g.copy()
            p=gl[i-1]
            pn=p+1
            left=g[i,0]
            right=g[i,-1]
            iter=1
            
            def di(h,tau,b,beta,sigma,y,yp,l):
                if l is 1:
                    return yp[l]-cil(h,tau,b,beta,sigma,y,l)*left-(b*beta*pow(y[l],beta-1)/h)*(y[l]-y[l-1])
                if l is n:
                    return yp[l]-ciu(h,tau,b,beta,sigma,y,l)*right-(b*beta*pow(y[l],beta-1)/h)*(y[l]-y[l-1])
                return yp[l]-(b*beta*pow(y[l],beta-1)/h)*(y[l]-y[l-1]) #I believe there is an issue with a missing tau here
            
            while 1>0:
                t1=[cil(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(2,n+1)]
                t2=[ciu(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n)]
                t3=[ci(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n+1)]
                d=nan_to_num([di(dx,dt,self.b,self.beta,self.sigma,p,g[i-1],j) for j in xrange(1,n+1)])
                a=nan_to_num(diagflat(t1,-1)+diagflat(t2,1)+diagflat(t3))
                #print a
                #print linalg.pinv(a,10**-6)
                try:
                    pn=concatenate([[left],linalg.solve(a,d),[right]])
                except linalg.LinAlgError:
                    try:
                        pn=concatenate([[left],dot(linalg.pinv(a),d),[right]])
                    except linalg.LinAlgError:
                        print "p:",p
                        print a
                        print d
                        pn=p*0
                        break                           
                
                if amax(fabs(pn-p))<steperr:
                    break
                if iter > itermax:
                    break
                iter += 1
                p=pn
            g[i]=pn

    def _calculate_python_1d_rev2(self,steperr=10**-5,itermax=20):
        """Implements scheme b in Samarskii's text with some slight modifications to give a 
            convection term rather than a reaction term.  Implements the convection as a 
            backward difference, so the convection terms appear with y_i and y_i-1 (ci and cil)
            For the stability of the scheme, we restrict the convection coefficient to
            be non-negative, but this still allows spurious oscillations in some expanding cases.
            
            This scheme is not rearranged or simplified by a computer algebra system. 
            """ 
        g=self.g.u
        n=len(g[0])-2
        def aiv(l,v):
            #return 0.5*(self.ktxu(v[i-1])+self.kxtu(v[i])))
            return 0.5*self.sigma*(pow(v[l-1],float(self.sigma)-1.0)+pow(v[l],float(self.sigma)-1.0))
            #return self.sigma*(pow((v[l-1]+v[l])*0.5,self.sigma-1))
        def kiv(l,v):
            if v[l]>0:
                if v[l+1]>0:
                    return 0.5*(pow(v[l],self.beta-1)+pow(v[l+1],self.beta-1))
                else:
                    return pow(v[l],self.beta-1)
            else:
                if v[l+1]>0:
                    return pow(v[l+1],self.beta-1)
                else:
                    return 0
#        def kiv(l,v):
#            return pow(v[l],float(self.beta)-1.0)                
        def ci(h,tau,b,beta,sigma,y,l):
            #print kiv(l,y)
            return 1+(tau/pow(h,2.0))*(aiv(l+1,y)+aiv(l,y))-((b*beta*kiv(l,y))/h)
            #return 1+((tau*aiv(l+1,y))/pow(h,2))+((tau*aiv(l,y))/pow(h,2))-(b*beta*pow(y[l],beta-1)/h)
        def cil(h,tau,b,beta,sigma,y,l):
            #return (b*beta*(pow(y[l],beta-1))/h)-((tau*aiv(l,y))/pow(h,2))
            #print kiv(l-1,y)
            return ((b*beta*kiv(l-1,y))/h)-((tau*aiv(l,y))/pow(h,2.0)) 
        def ciu(h,tau,b,beta,sigma,y,l):
            return -tau*aiv(l+1,y)/pow(h,2.0)
            
        tax=self.g._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]#FIXME: There should be a better way to do this
        #FIXME: Assumes constants dt and dx: uniform grids only!
        #FIXME: When possible, reduce the number of arguments to 
        dt=abs(self.g.x[0][1]-self.g.x[0][0])
#        print 'dt:',dt
        
        dx=abs(self.g.x[1][1]-self.g.x[1][0])
#        print 'dx:',dx

        for i in self.g.get_ind_range(ax=tax,aset='active'):
            gl=g.copy()
            p=gl[i-1]
            pn=p+1
            left=g[i,0]
            right=g[i,-1]
            iter=1
            
            def di(h,tau,b,beta,sigma,y,yp,l):
                if l is 1:
                    return yp[l]-cil(h,tau,b,beta,sigma,y,l)*left
                if l is n:
                    return yp[l]-ciu(h,tau,b,beta,sigma,y,l)*right
                return yp[l]
            
            while 1>0:
                t1=[cil(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(2,n+1)]
                t2=[ciu(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n)]
                t3=[ci(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n+1)]
                d=nan_to_num([di(dx,dt,self.b,self.beta,self.sigma,p,g[i-1],j) for j in xrange(1,n+1)])
                a=nan_to_num(diagflat(t1,-1)+diagflat(t2,1)+diagflat(t3))
                #print a
                #print linalg.pinv(a,10**-6)
                try:
                    #pn=concatenate([[left],dot(linalg.pinv(a,10**-6),d),[right]])
                    pn=concatenate([[left],linalg.solve(a,d),[right]])
                except linalg.LinAlgError:
                    self.g.l.warning("System is numerically singular; finding least-squares solution using pseudo-inverse.")
                    self.g.l.info('Encoutered at iteration {}'.format(iter))
                    self.g.l.debug("Last guess: {}\n Coefficient Matrix: {}\n RHS: {}".format(p,a,d))
                    try:
                        pn=concatenate([[left],dot(linalg.pinv(a,10**-6),d),[right]])
                    except linalg.LinAlgError:
                        #TODO: Make these work with logging system
                        print "p:",p
                        print a
                        print d
                        pn=p*0
                        raise
                
                if amax(fabs(pn-p))<steperr:
                    self.g.l.info("Accuracy goal reached for i={} in {} iterations.".format(i,iter))
                    break
                if iter > itermax:
                    self.g.l.info("Iteration maximum reached for i={}.  Final absolute error: {}".format(i, amax(fabs(pn-p))))
                    break
                iter += 1
                p=pn
            g[i]=pn

    def _calculate_python_1d_rev1(self,steperr=10**-7,itermax=50):
        #TODO: This scheme can also benefit from the changes made to rev2
        """Implements scheme b in Samarskii's text with some slight modifications to give a 
            convection term rather than a reaction term.  Implements the convection as a 
            backward difference, so the convection terms appear with y_i and y_i-1 (ciu and cil)
            
            This scheme has been simplified by a computer algebra system, which seems to negatively
            impact the numerical stability of the scheme.  Given the issues appearing with this
            method and the experiments we have run, it seems prudent to discard this scheme completely.
            """ 
        def aiv(i,v):
            #return 0.5*(self.ktxu(v[i-1])+self.kxtu(v[i])))
            return 0.5*(pow(v[i-1],self.sigma)+pow(v[i],self.sigma))
            #return pow(0.5*(v[i-1]+v[i]),self.sigma)
        def denom(h,tau,b,beta,sigma,y,l):
            return pow(h,2)+tau*(
                                    aiv(l+1,y)+aiv( l,y)-b*beta*h*pow(y[l],beta-1)
                                    )
        def cil(h,tau,b,beta,sigma,y,l):
            #return l
            #print "cil(",l,"): ",aiv(l+1,y)
            return tau*(aiv(l+1,y)-h*b*beta*pow(y[l],beta-1))/denom(h,tau,b,beta,sigma,y,l)
        def ciu(h,tau,b,beta,sigma,y,l):
            #return l
            #print "ciu(",l,"): ",aiv(l,y)
            return tau*aiv(l,y)/denom(h,tau,b,beta,sigma,y,l)
            
        tax=self.g._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]#FIXME: There should be a better way to do this
        #FIXME: Assumes constants dt and dx: uniform grids only!
        g=self.g.u
        dt=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((1,0)))
        dx=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((0,1)))
        n=len(g[0])-2

        for i in self.g.get_ind_range(ax=tax,aset='active'):
            gl=g.copy()
            p=gl[i-1]
            pn=p+1
            left=g[i,0]
            right=g[i,-1]
            iter=1
            
            def di(h,tau,b,beta,sigma,y,yp,l):
                
                if l is 1:
                    return cil(h,tau,b,beta,sigma,y,l)*left+(pow(h,2)*yp[l]/denom(h,tau,b,beta,sigma,y,l))
                if l is n:
                    return ciu(h,tau,b,beta,sigma,y,l)*right+(pow(h,2)*yp[l]/denom(h,tau,b,beta,sigma,y,l))
                return (pow(h,2)*(0.05*yp[l-1]+0.9*yp[l]+0.05*yp[l+1])/denom(h,tau,b,beta,sigma,y,l))
            
            while 1>0:
                t1=[cil(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(2,n+1)]
                t2=[ciu(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n)]
                d=[di(dx,dt,self.b,self.beta,self.sigma,p,g[i-1],j) for j in xrange(1,n+1)]
                a=eye(n)-diagflat(t1,-1)-diagflat(t2,1)
                pn=concatenate([[left],linalg.solve(a,d),[right]])
                #i=0
                #for k in pn:
                #    print i,": ",k
                #    i+=1 
                #print pn
                if amax(fabs(pn-p))<steperr:
                    break
                if iter > itermax:
                    break
#                print "iteration: ",iter," for slice i=",i
#                print amax(fabs(pn-p))," ",amax(fabs(pn-p))/amax(fabs(pn))
                iter += 1
                p=pn
                #exit()
            #print iter    
            iter=1
            g[i]=pn

    def _calculate_python_1d_schemea(self):
        """Implements scheme a in Samarskii's text with some slight modifications to give a 
            convection term rather than a reaction term.  Implements the convection as a 
            backward difference.
            
            This scheme has been simplified by a computer algebra system, which may detrimentally
            effect the numerical stability. 
            """ 
            #TODO: Rewrite this scheme in a more readable form.
        tax=self.g._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]#FIXME: There should be a better way to do this
        #print numpy.eye(5)-(numpy.diagflat(2*numpy.ones((1,4)),k=1)+numpy.diagflat(numpy.ones((1,4)),k=-1))
        #FIXME: Assumes constants dt and dx: uniform grids only!
        g=self.g.u
        dt=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((1,0)))
        dx=self.g._dist(self.g.ind_to_x((0,0)),self.g.ind_to_x((0,1)))
        for i in self.g.get_ind_range(ax=tax,aset='active'):
            #print g
            left=g[i,0]
            right=g[i,-1]
            p=g[i-1]
            #print "l=",l," r=",r
            #print p
            n=len(p)-2
            def denom(h,tau,b,beta,sigma,y,l):
#                return pow(h,2)+tau*(
#                                    pow(0.5*(y[l]+y[l+1]),sigma)+pow(0.5*(y[l]+y[l-1]),sigma)-b*beta*h*pow(y[l],beta-1)
#                                    )
                return pow(h,2)+tau*(
                                    pow(y[l],sigma)+0.5*(pow(y[l+1],sigma)+pow(y[l-1],sigma))-b*beta*h*pow(y[l],beta-1)
                                    )
            def cil(h,tau,b,beta,sigma,y,l):
                #print "cil called for i=",l 
                #return tau*(pow(0.5*(y[l]+y[l+1]),sigma)-h*b*beta*pow(y[l],beta-1))/denom(h,tau,b,beta,sigma,y,l)
                return tau*(0.5*(pow(y[l],sigma)+pow(y[l+1],sigma))-h*b*beta*pow(y[l],beta-1))/denom(h,tau,b,beta,sigma,y,l)
                #return l+20
            def ciu(h,tau,b,beta,sigma,y,l):
                #print "ciu called for i=",l
                #return tau*(pow(0.5*(y[l-1]+y[l]),sigma))/denom(h,tau,b,beta,sigma,y,l)
                return tau*(0.5*(pow(y[l-1],sigma)+pow(y[l],sigma)))/denom(h,tau,b,beta,sigma,y,l)
                #return l+10
            def di(h,tau,b,beta,sigma,y,l):
                #print "di called for i=",l
                add=0
                if l is 1:
                    add=cil(h,tau,b,beta,sigma,y,l)*left
                if l is n:
                    add=ciu(h,tau,b,beta,sigma,y,l)*right
                return add+(pow(h,2)*y[l]/denom(h,tau,b,beta,sigma,y,l))
                #return l
            #print n
            t1=[cil(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(2,n+1)]
            t2=[ciu(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n)]
            #print t1
            #print t2
            d=[di(dx,dt,self.b,self.beta,self.sigma,p,j) for j in xrange(1,n+1)]
            a=eye(n)-(diagflat(t1,-1)+diagflat(t2,1))
            #print a
            #print d
            #print g[i,1:-1]=
            g[i,1:-1]=linalg.solve(a,d)