import numpy
from ..log import plog
#TODO: Allow all grid objects to emulate numeric types
#TODO: Allow all grid objects to emulate container types
#TODO: Build reasonable special functions for grid objects (i.e. __str__(), etc)
#TODO: See http://docs.python.org/reference/datamodel.html
class trivgrid(object):
    def __init__(self,u=None,nt=None):
        if u is None:
            u=numpy.zeros((2,2),'d')
        if nt is None:
            nt=self.__class__
        self.u=u
        self.l=plog(str(nt))
        
    def l2_norm(self):
        return self._l2_norm(self.u)

    def max_norm(self):
        return self._max_norm(self.u)

    @staticmethod
    def _l2_norm(u):
        return pow(numpy.sum(numpy.square(u)), 1/2)
        
    @staticmethod
    def _max_norm(u):
        return numpy.amax(numpy.fabs(u))

class grid(trivgrid):
    """ Base grid object """
    def __init__(self, u=numpy.zeros((2, 2), 'd'), active=numpy.zeros((1, 2), 'i')):
        super(grid, self).__init__(u,self.__class__)
        self.l.info(self.l.newinst_info('grid',self.__class__))
        self.active = active
        self.left = None
        self.timelike = None
        self.l.info('Grid has {} dimensions.'.format(self.u.ndim))
        self.l.debug('Grid has dimensions {} for a total size {}.'.format(self.u.shape, self.u.size))
        
    def load(self, fname):
        self.l.info('Loading grid from {}'.format(fname))
        tfile = numpy.load(fname)
        self.u, self.active = tfile['u'], tfile['active']

    def save(self, fname=None, fmt=1):
        if fname is None:
            self.l.warning('Saving to temporary file')
            from tempfile import TemporaryFile
            fname = TemporaryFile()
        
        self.l.info('Saving grid to {}'.format(fname))
        if fmt is 1:
            numpy.savez(fname, u=self.u, active=self.active)
        else:
            numpy.savetxt(fname,self.u)
            self.l.warning('Information about active region will not be exported to text file.')
        return fname
    
    def get_ind_range(self, ax, aset='all'):
        if aset == 'all':
            return range(self.u.shape[ax]);
        elif aset == 'active':
            t = set()
            for i in self.active:
                t.add(i[ax])
            return t

    def _set_active_xv(self, xv):
        """
        Internal:   Set the domain interior (self.active) through 
                    a list of x values xv"""
        self._set_active_ind(map(self.x_to_ind, xv))

    def _set_active_ind(self, indices):
        """Internal:    Set the domain interior (self.active) directly"""
        self.l.info('Setting active area.')
        self.l.debug('Corners of active area bounding box: {} {}'.format(numpy.amax(indices,0),numpy.amin(indices,0)))
        numpy.transpose(indices)
        self.active = indices

    def _timelike_step(self, prevind, delta=1):
        return self._ax_step(prevind, self._get_timelike_axis(), delta)
    
    @classmethod
    def _ax_step(cls, prevind, ax, delta=1):
        #FIXME: Now that this is a class method, there should be a method
        # which does bound checking
        if isinstance(ax, type(1)):
            e = cls._get_basis_elem_static(ax, len(prevind))
        else:
            raise TypeError(plog.type_err('ax',type(ax),'int'))
        return cls._add(prevind, cls._scale(e, delta))

    def _set_timelike_axis(self, axis):
        if isinstance(axis, type(1)):
            self.l.debug('Timelike axis is {}'.format(axis))
            self.timelike = axis
        else:
            raise TypeError(self.l.type_err('axis','int',type(axis)))
        
    def _get_timelike_axis(self):
        if self.timelike is None:
            self.l.warning('Grid has not defined a timelike axis. Guessing zero, which seems likely.')
            return 0
        else:
            return self.timelike

    @staticmethod
    def _scale(ind, scale):
        ind = list(ind)
        for i in range(len(ind)):
            ind[i] = ind[i] * scale
        return tuple(ind)

    @staticmethod
    def _add(ind1, ind2, zerofill=False):
        #TODO: Add capability for adding vectors of different lengths by filling to larger size with zeros
        if len(ind1) == len(ind2):
            ind1 = list(ind1)
            ind2 = list(ind2)
            for i in range(len(ind1)):
                ind1[i] += ind2[i]
            return tuple(ind1)
        else:
            raise ValueError(
                        'Index dimension mismatch: ind1 has {}dimensions, ind2 has {} dimensions'.format(len(ind1), len(ind2))
                        )
    
    @classmethod
    def _sub(cls, ind1, ind2):
        return cls._add(ind1, cls._scale(ind2, -1.0))
        
    @classmethod
    def _dist(cls, ind1, ind2):
        #print numpy.hypot(ind1,ind2)
        #print sum(x**2 for x in cls._sub(ind1,ind2))**(0.5)
        #TODO: Should be able to use (faster?) numpy linalg package command for this.
        return sum(x ** 2 for x in cls._sub(ind1, ind2)) ** (0.5)
    
    def _get_basis_elem(self, axis):
        return self.__get_basis_elem(axis, self.u.ndim)
    
    @staticmethod
    def _get_basis_elem_static(axis, N):
        def l(iv):
            if iv == axis:
                return 1
            else:
                return 0
        return tuple([l(i) for i in range(N)])

    def get_plane_ind(self, val, ax='timelike', aset='active'):
        self.l.debug('Giving plane {} in {} axis from set {}'.format(val,ax,aset))
        if aset == 'active':
            ind = self.active
        elif aset == 'bdy':
            ind = self.get_bdy_ind()
        elif aset == 'left':
            ind = self.left
        else:
            ind = numpy.transpose(
                                numpy.indices(self.u.shape)
                                ).reshape((numpy.prod(self.u.shape), self.u.ndim))
#       print ind
        if ax == 'timelike':
            ax = self._get_timelike_axis()
        elif not isinstance(ax, type(1)):
            raise TypeError(self.l.type_err('ax','int',type(ax)))
        t = []
        for i in ind:
            if i[ax] == val:
                t.append(i)
        if len(t) > 0:
            return map(tuple, numpy.vstack(t[:]))
        else:
            return []

    def _set_left_to_active(self):
        self.l.info('Set of points left replaced by active set.')
        self.left = self.active.copy()

    def get_time_coord(self, ind):
        return ind[self._get_timelike_axis()]

    def set_entrance(self, ind):
        if self._check_entrance(ind):
            self.l.debug('Grid entrance set to {}'.format(ind))
            self.entrance = ind

    def _check_entrance(self, ind):
        #TODO: Determine whether there is a more accurate general test
        return ind in self.active

    def l2_dist(self, v):
        return self._l2_norm(self.u - v)

    def max_dist(self, v):
        return self._max_norm(self.u - v)

    def x_to_u(self, x):
        raise NotImplementedError('Grid does not implement x_to_u.')

    def ind_to_u(self, ind):
        if ind in numpy.indices(self.u.shape):
            return self.u[tuple(ind)]
        raise NotImplementedError(
                                'No interpolation implemented in base grid class: ind_to_u not implemented'
                                )
    
    def ind_to_x(self, ind):
        raise NotImplementedError('Grid does not implement ind_to_x.')

    def x_to_ind(self, x):
        raise NotImplementedError('Grid does not implement x_to_ind.')
    
    def get_active_mat(self):
        """
        Return the active indices in the form of a matrix where the value
        1 is at the active indices and 0 elsewhere. This this matrix
        has the form that the setter function set_active_mat would expect
        """ 
        return self._active_to_mat(self.u.shape, self.active)

    def get_bdy_mat(self):
        return 1 - self.get_active_mat()

    def get_bdy_ind(self):
        return self._active_to_bdy(self.u.shape, self.active)
        
    @classmethod
    def _active_to_bdy(cls, mshape, active):
        return numpy.transpose(numpy.nonzero(1 - cls._active_to_mat(mshape, active)))

    @classmethod
    def _active_to_mat(cls, mshape, active):
        """
        Converts the list of indices active and matrix shape mshape to a 
        matrix of shape mshape with value 1 at all the points referenced
        by active.
            
        input:  shape tuple mshape for the matrix
                list of indices active for the active region
        
        output: a matrix where the value 1 is at the active indices 
                and 0 elsewhere
        """
        m = numpy.zeros(mshape)
        for ind in active:
            m[ind[0], ind[1]] = 1 # Should be able to rewrite with numpy.put
        return m

class active_manip(object):
    def set_active_mat(self, m):
        """
        Set domain interior (self.active)
        input:  N-dimensional array m which has nonzero entries
                in the interior of the domain
        """
        self._set_active_ind(numpy.transpose(numpy.nonzero(m)))
    
    def set_active_all(self):
        """Set domain interior to be whole domain."""
        self.l.info('Setting active area to whole domain.')
        self._set_active_ind(numpy.transpose(
                                            numpy.indices(self.nx)
                                            ).reshape((numpy.prod(self.nx), 2)))
    
    def set_active_int(self):
        """Set domain interior  to be all but extreme points."""
        self.l.info('Setting active area to all but extreme points.')
        def domfn(x):
            """Function which is true so long as the x value is not at the edge of the allowed values."""
            for (xi, rxi) in zip(x, self.rx):
                if xi <= rxi[0] or xi >= rxi[1]:
                    return False
            return True
        self.set_active_eq(domfn)
        
    def set_active_int_plus_top(self):
        self.l.info('Setting active area to be whole domain minus the top')
        ax = self._get_timelike_axis()
        def domfn(x):
            """
            Function which is True on the extreme points of the domain,
            excepting where the timelike axis takes maximum value (the "lid".)
            """
            for (i, rxi) in zip(range(0, self.u.ndim), self.rx):
                if i == ax and x[i] == rxi[0]:
                    return False
                elif x[i] <= rxi[0] or x[i] >= rxi[1]:
                    return False
            return True 
        self.set_active_eq(domfn)

    def set_active_eq(self, domfn):
        """
        Set domain interior (self.active) through domfn
        input:  A function domfn(x) which returns True if x should be in
                active, and False if x should not be in active
        """
        if self.u.ndim == 2:
            row, col = numpy.indices(self.nx)
            td = [self.ind_to_x([ind1, ind2]) for (ind1, ind2) in zip(numpy.concatenate(row), numpy.concatenate(col))]
        else:
            NotImplementedError('set_active_eq not defined for grids with ndim != 2')
            
        tmp = [i for i in td if domfn(i)]
        self._set_active_xv(numpy.vstack(tmp))
        
    def set_active_by_bdy(self, bdy):
            #TODO: Implement the ability to set the active set by its boundary
            raise NotImplementedError('set_active_by_bdy not yet implemented')
