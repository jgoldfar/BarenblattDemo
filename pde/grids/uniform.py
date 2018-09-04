"""
Collection of "uniform" grids in n-dimensions, with distance in each
direction parameterized by the vector h
"""
import numpy
from .base import grid, active_manip

class uniformgrid(grid,active_manip):
    """Base class for uniform grids"""
    def __init__(self, dims, type='d', x=None, rx=None ):
        """
        Constructor for a uniform grid.
        Creates an solution grid self.u with shape dims, of type type initialized to zero 
        everywhere and an active set self.active of appropriate number of dimensions also 
        initialized to zero. Note that the active set by default limits index values to
        be 16  bit integers, which gives us a large space over which we may calculate the
        solution, but could run out for enormous grid sizes.
        Also, sets self.t to be the type of the solution grid for future reference.
        """
        #TODO: This implementation can be made faster/leaner by not storing all of the x values
        self.nx=list(map(int,dims))
        super(uniformgrid, self).__init__(numpy.zeros(dims, type), numpy.zeros((1,2),'int16'))

        if x is None:
            self.x=map(lambda xi:numpy.linspace(0,1,xi,True),dims)
        else:
            self.x=x
            
        if rx is None:
            self.rx=map(lambda xi:[numpy.amin(xi),numpy.amax(xi)],self.x)
        else:
            self.rx=rx

    def ind_xdist(self,ind1,ind2):
        return self._dist(self.ind_to_x(ind1), self.ind_to_x(ind2))

    def ind_to_x(self, ind):
        """
        input: an index ind
        output: the corresponding x value
        """
        return tuple([self.x[i][ind[i]] for i in range(0,self.u.ndim)])
        
    def x_to_ind(self, x):
        """
        input:  an x value. Note that if the x value is outside the solution domain,
                this function will throw a ValueError
        output: the closest corresponding index in the solution grid
        """
        if len(x) > self.u.ndim:
            raise ValueError('x has {} dimensions, but grid has {}.'.format(len(x),self.u.ndim))
        
        for (xi, rxi) in zip(x, self.rx):
            if xi < rxi[0] or xi > rxi[1]:
                raise ValueError(
                                '{0} exceeds the range of available x values \
                                : [{1},{2}]'.format(xi, rxi[0], rxi[1]))
        return [numpy.argmin(numpy.absolute(x[i]-self.x[i])) for i in range(0,self.u.ndim)]

    def set_bc(self,fn):
        """
        Setting the values on the parabolic boundary (or in the inactive segment
        of the domain) is better done here; the input fn should be defined everywhere,
        or at least in the complement of the active indices in the domain. Setting the
        boundary values in this manner is prefered from a conceptual and practical standpoint
        """
        self.l.info('Setting values on parabolic boundary.')
        self.l.debug('Function docstring: {}'.format(fn.__doc__))
        if not hasattr(fn,'__call__'):
            raise ValueError('set_bc called with a nonfunction argument')
        bcind=numpy.transpose(numpy.nonzero(numpy.ones(self.u.shape)-self.get_active_mat()))
        for ind in bcind:
            self.u[ind[0],ind[1]]=fn(self.ind_to_x(ind))

class twoddir(uniformgrid):
    """Test case for a two dimensional grid"""
    def __init__(self, h=(0.1,0.1), rx=[(0,1),(0,1)], n=None):
        #TODO: More clearly give option to choose number of grid points
        if len(h) != len(rx):
            raise ValueError(
                        'Grid spacing vector h has length {} while grid range matrix xr suggests {}-dimensions\
                        '.format(len(h),len(rx))
                            )
        self.rx = rx
        self.h = h
        super(twoddir, self).__init__( [int((rxi[1]-rxi[0])/hi) for (rxi, hi) in zip(rx, h)])
        self.x=[numpy.linspace(rxi[0], rxi[1], nxi, True) for (rxi, nxi) in zip(rx, self.nx)]
        self._set_timelike_axis(0)

    def _set_bc_left(self,fn):
        """
        Internal: function which sets boundary values according to fn on the left side
        of the boundary (u[:,0])
        """
        for i in range(0,len(self.u[:,0])):
            self.u[i,0]=fn(self.ind_to_x((i,0)))
    def _set_bc_right(self,fn):
        """
        Internal: function which sets boundary values according to fn on the right side
        of the boundary (u[:,-1])
        """
        for i in range(0,len(self.u[:,-1])):
            self.u[i,-1]=fn(self.ind_to_x((i,-1)))
        
    def set_lat_bc(self,*args):
        """
        Set function values on the lateral boundaries (u[0,:] and u[-1,:].)
        Takes a variable number of arguments depending on the preferences of the user.
        If one argument is given, it should be a function defined everywhere in the domain
        (or at least on the complement of the active indices in the domain) which gives
        the boundary value there. If two arguments are given, they should be the left and
        right boundary values, respectively. Mixing these is a very bad idea.
        
        As explained in the documentation for set_bc, the complication and chance of
        making an error is much increased when using set_lat_bc, so it should be avoided.
        """
        for i in range(0,len(args)):
            if not hasattr(args[i],'__call__'):
                raise TypeError('set_lat_bc called with at least one nonfunction argument')
        if len(args) is 2:
            self._set_bc_left(args[0])
            self._set_bc_right(args[1])
        elif len(args) is 1:
            self._set_bc_left(args[0])
            self._set_bc_right(args[0])
        else:
            raise ValueError('1 or 2 arguments expected to set_lat_bc, {} given.'.format(len(args)))
        
    def set_ic(self,fn):
        """Set the initial condition (where row#=0) according to fn"""
        for i in range(0,len(self.u[0,:])):
            self.u[0,i]=fn(self.ind_to_x((0,i)))
