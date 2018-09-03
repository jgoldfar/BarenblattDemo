"""Explicit methods to calculate PDE solutions on grids"""
from base import method, differencer
#from numpy import eye, diagflat, linalg, amax, fabs,concatenate
#from scipy.weave import converters

class diffusion(method):
	def __init__(self,grid=None,b=0,beta=1,alpha=None,ktxu=None):
		super(diffusion,self).__init__(grid)
		if ktxu is None:
			if alpha is None:
				def tmp(t,x,u):
					return 1
			else:
				def tmp(t,x,u):
					return alpha 
			self.ktxu = tmp
		else:
			self.ktxu=ktxu
		self.b=b
		self.beta=beta
		self.alpha=alpha
		
	@staticmethod
	def _grid_allows_stencil(grid):
		#TODO: Implement grid checking
		super(diffusion,self)._grid_allows_stencil()

	def calculate(self):
		if self.g.u.ndim==2:
			self._calculate_python_1d()
		else:
			raise NotImplementedError('Grids with more than one spatial dimension not supported.')
				
	def _calculate_python_1d(self):
		tax=self.g._get_timelike_axis()
		oax=[l for l in range(self.g.u.ndim) if l != tax]
		oax=oax[0]
		#print self.g.u
		for i in self.g.get_ind_range(ax=tax,aset='active'):
			#TODO: This can be written as a matrix multiplication (much quicker!)
			#Note that this algorithm was simplified by the help of a computer algebra
			# system; this may actually decrease the numerical accuracy of the 
			# scheme.
			curr=self.g.get_plane_ind(val=i,ax=tax,aset='active')
			for ind in curr:
				ua=self.g._ax_step(ind,tax,-1); xi=self.g.ind_to_x(ua)
				ub=self.g._ax_step(ua,oax,-1); xim=self.g.ind_to_x(ub)
				uc=self.g._ax_step(ua,oax,1); xip=self.g.ind_to_x(uc)
				dt=self.g._dist(self.g.ind_to_x(ind),xi)
				dx1=self.g._dist(xim,xi)
				#dx2=self.g._dist(xip,xi)
				yi=self.g.u[ua];t=xi[tax]
				yip=self.g.u[uc]
				yim=self.g.u[ub];
				self.g.u[ind]=(dt*(yim-yi + yip)*self.ktxu(t,xi[oax],yi) + dt*yim*self.ktxu(t,xim[oax],yim) + dt*(yip-yi)*self.ktxu(t,xip[oax],yip) + 2*dx1*(dx1*yi + self.b*self.beta*dt*( yim-yi)*pow(yi,self.beta-1)))/(2*pow(dx1,2) + dt*(self.ktxu(t,xi[oax],yi) + self.ktxu(t,xim[oax],yim)))
	
	def _calculate_python_1d_lin(self):
		tax=self.g._get_timelike_axis()
		oax=[l for l in range(self.g.u.ndim) if l != tax]
		oax=oax[0]
		dt=self.g._dist(self.g.ind_to_x((0,0)),self.g._ax_step((0,0),tax,1))
		dx1=self.g._dist(self.g.ind_to_x((0,0)),self.g._ax_step((0,0),oax,1))
		for i in self.g.get_ind_range(ax=tax,aset='active'):
			curr=self.g.get_plane_ind(val=i,ax=tax,aset='active')
			for ind in curr:
				ua=self.g._ax_step(ind,tax,-1); xa=self.g.ind_to_x(ua)
				ub=self.g._ax_step(ua,oax,-1)
				uc=self.g._ax_step(ua,oax,1)
				self.g.u[ind]=self.g.u[ua]+((self.alpha*dt)/(pow(dx1,2)))*(self.g.u[ub]+self.g.u[uc]-2*self.g.u[ua])


class ftcs_1d(method):
	#TODO: Implement FTCS for objective functions of form F(t,x,u,ux1,...,uxn,ux1x1,...,uxnxn)
	# We should be able to write a general method for equations of the form
	#	u_t=F(t,x,u,u_x,u_xx)
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