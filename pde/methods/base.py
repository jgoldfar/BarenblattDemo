""" Base method """

from ..grids.base import grid as gr
from ..log import plog

class method(object):
	def __init__(self, grid=None):
		super(method,self).__init__()
		self.l=plog(str(self.__class__))
		self.l.info(self.l.newinst_info('method',self.__class__))
		if grid is not None:
			self.set_grid(grid)

	def set_grid(self,grid):
		if issubclass(type(grid),gr):
			self.g=grid
		else:
			raise TypeError(self.l.type_err('grid',type(gr),type(grid)))

	def __call__(self):
		self.calculate()

	def _grid_allows_stencil(self):
		raise NotImplementedError('Method does not check if the grid\'s active area is valid.')
	
	def calculate(self):
		raise NotImplementedError('Method has no calculation subroutine.')
	
	def _get_dep_cone(self,ind):
		raise NotImplementedError('Method cannot give the cone of dependence at a point.')
	
class differencer(object):
	"""Contains standard class methods to calculate differences on grids."""
	#TODO: Test whether calculating differences on a whole slice will save time
	#TODO: Test whether using weave here has any timing advantages
	#TODO: Build calculation of dx/dt into this code
	#TODO: Make calculation of dx/dt quicker for uniform grids (detect subclassing)
	@staticmethod
	def _indleft(grid,axis,ind):
		return tuple(grid._ax_step(ind,axis,-1))
	
	@staticmethod
	def _indright(grid,axis,ind):
		return tuple(grid._ax_step(ind,axis,1))
	
	@classmethod
	def _uleft(cls,grid,axis,ind):
		return grid.u[cls._indleft(grid,axis,ind)]
	
	@classmethod
	def _uright(cls,grid,axis,ind):
		return grid.u[cls._indright(grid,axis,ind)]

	@classmethod
	def _central_diff(cls,grid,axis,ind,dh=None):
		il=cls._indleft(grid,axis,ind)
		ir=cls._indright(grid,axis,ind)
		dx1=grid._dist(grid.ind_to_x(il),grid.ind_to_x(ind))
		dx2=grid._dist(grid.ind_to_x(ir),grid.ind_to_x(ind))
		return (grid.u[ir]-grid.u[il])/(dx1+dx2)
	
	@classmethod
	def _back_diff(cls,grid,axis,ind,dh=None):
		il=cls._indleft(grid,axis,ind)
		dx=grid._dist(grid.ind_to_x(il),grid.ind_to_x(ind))
		return (grid.u[tuple(ind)]-grid.u[il])/dx
	
	@classmethod
	def _forward_diff(cls,grid,axis,ind,dh=None):
		ir=cls._indright(grid,axis,ind)
		dx=grid._dist(grid.ind_to_x(ir),grid.ind_to_x(ind))
		return (grid.u[ir]-grid.u[tuple(ind)])/dx
	
	@classmethod
	def dt_backward(cls,grid,ind):
		return cls.dx(grid,grid._get_timelike_axis(),ind,0)
	
	@classmethod
	def dx(cls,grid,axis,ind,type=1):
		if type==2:
			return cls._forward_diff(grid,axis,ind)
		elif type==1:
			return cls._central_diff(grid,axis,ind)
		elif type==0:
			return cls._back_diff(grid,axis,ind)
	
	@classmethod
	def uxx_centered(cls,grid,axis,ind):
		il=cls._indleft(grid, axis, ind)
		ir=cls._indright(grid,axis,ind)
		dx1=grid._dist(grid.ind_to_x(ir),grid.ind_to_x(ind))
		dx2=grid._dist(grid.ind_to_x(il),grid.ind_to_x(ind))
		if abs(dx1-dx2)>1e-10:
			raise NotImplementedError('Nonuniform grids not currently supported for uxx-central difference.\n\Grid distance difference: {}'.format(dx1-dx2))
		return (1/(dx1**2))*(grid.u[il]+grid.u[ir]-2*grid.u[ind])