"""Provides plotting routines for grid objects"""
from numpy import zeros, tile, linspace, transpose
import matplotlib.pyplot as plt
#from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from grids.base import trivgrid
from .log import plog

class plotobj(object):
	def __init__(self,bk=None,fname=None):
		self.fname=fname
		super(plotobj,self).__init__()
		self.l=plog(self.__class__)
		
		
	def save(self, fname=None,dpi=None):
		if fname is None:
			if self.fname is None:
				raise ValueError('No filename given.')
		ext=fname[-3:len(fname)]
		if ext=='svg':
			self.l.info('Switching backend to SVG')
			plt.switch_backend('SVG')
			fmt='svg'
		elif ext=='.ps':
			self.l.info('Switching backend to PS')
			plt.switch_backend('PS')
			fmt='ps'
		elif ext=='eps':
			self.l.info('Switching backend to PS to save eps file')
			plt.switch_backend('PS')
			fmt='eps'
		elif ext=='pdf':
			self.l.info('Switching backend to PDF')
			plt.switch_backend('PDF')
			fmt='pdf'
		elif ext=='png':
			self.l.info('Default backend selected to save png file')
			fmt='png'
		
		self.l.info('Saving current plots to {}'.format(fname))
		self._save(fname,fmt)

	@staticmethod
	def _save(fname,fmt):
		plt.savefig(fname, format=fmt,bbox_inches='tight')
	
	@staticmethod
	def _start():
		plt.close('all')
		return plt.figure()
	
	def start(self):
		self.l.warning('Closing all plots.')
		self.l.info('Starting new plot queue')
		return self._start()
	
	@staticmethod
	def _show():
		plt.show()
		
	def show(self):
		self.l.info('Showing queued plots.')
		self._show()

class ddplot(plotobj):
	def __init__(self,grid=None,file=None):
		if grid is not None:
			self.set_grid(grid)
		self.fname=file
		super(ddplot,self).__init__()
	
	def set_grid(self,grid):
		if grid.u.ndim is not 1:
			raise ValueError(
								'Grid with 1 dimension expected, but input has \
								{} dimensions'.format(grid.u.ndim)
							)
		else:
			self.g=grid

	@staticmethod
	def _plot(x,y, xlabel='',ylabel='',title=''):
		plt.plot(x,y)

class dddplot(plotobj):
	def __init__(self,grid=None):
		if grid is not None:
			self.set_grid(grid)
		super(dddplot,self).__init__()

	def set_grid(self,grid):
		if grid.u.ndim is not 2:
			raise ValueError(
								'Grid with 2 dimensions expected, \
								but input has {} dimensions'.format(grid.u.ndim)
							)
		else:
			self.g=grid



class matplot(plotobj):
	def __init__(self,grid=None):
		if grid is not None:
			self.g=grid
		super(matplot,self).__init__()
	
	def spy(self):
		self._spy(self.g)

	def surf(self):
		self._surf(self.g)
		
	def heatmap(self):
		self._heatmap(self.g)
	
	def heatmap_grayscale(self):
		self._heatmap_grayscale(self.g)
	
		
	@staticmethod
	def _surf(g):
		g.l.info('Matplot: Creating surface plot for grid.')
		f=plt.figure()
		
		ax=f.add_subplot(111,projection='3d')
		try:
			x=g.x
		except AttributeError:
			x=[linspace(0,1,g.u.shape[i],True) for i in range(0,g.u.ndim)]

		g.l.debug('Matplot: Creating mesh grids from x value ranges of length {}.'.format(
																				tuple(len(i) for i in x)
																				)
				)
		xx,yy=tile(x[1],(len(x[0]),1)),tile(transpose([x[0]]),(1,len(x[1])))
		g.l.info('Matplot: Mesh grids created of shape {}.'.format(xx.shape))
		ax.plot_wireframe(yy,xx,g.u)
		return ax

	@staticmethod
	def _spy(g,marker='o'):
		g.l.info('Matplot: Visualizing sparsity pattern for grid.')
		plt.spy(g.u, marker=marker)

	@classmethod
	def _spy_active(cls, mshape, active):
		m=zeros(mshape)
		for ind in active:
			m[tuple(ind)]=1
		cls._spy(trivgrid(m,'matplot.active.tmp'))
#TODO: Make the signatures for active-set viewing agree with other methods

	@staticmethod
	def _print_active(mshape,active):
		m=zeros(mshape)
		for ind in active:
			m[tuple(ind)]=1
		print m

	@staticmethod
	def _heatmap(g):
		g.l.info('Showing heat map for grid with default colormap')
		plt.matshow(g.u)
	
	@staticmethod
	def _heatmap_grayscale(g):
		raise NotImplementedError('_heatmap_grayscale(g) not implemented.')
#		g.l.info('Showing heat map for grid with grayscale colormap')
#		plt.matshow(g.u,cmap=cm.gray)