import numpy
import math
#TODO: Make sure these work with grid instances as simple as trivgrid (when necessary)

class int_track(object):
    def __init__(self,grid,tol=10**-10):
        #TODO: Add more interface tracking options and support ndim>2
        if grid.u.ndim is not 2:
            raise ValueError('Interface tracker currently supports 1 spatial dimension.')
        grid.l.info('Tracking interface')
        tax=grid._get_timelike_axis()
        right_edgex=[]
        for ti in grid.get_ind_range(ax=tax,aset='all'):
            vals = [grid.ind_to_x([ti_dup, xi])[1] for (ti_dup, xi) in grid.get_plane_ind(val=ti, ax=tax, aset='all') if abs(grid.u[ti_dup, xi]) > tol]
            if vals:
                right_edgex.append(numpy.amax(vals))
            else:
                right_edgex.append(float('nan'))
        self.rex=right_edgex
