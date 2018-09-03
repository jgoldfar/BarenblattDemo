import numpy
#TODO: Make sure these work with grid instances as simple as trivgrid (when necessary)
class spt(object):
    """
    Define support of a given `grid` with tolerance `tol`
    """
    def __init__(self,grid,tol=10**-10):
        super(spt,self).__init__()
        self.g=grid
        self.s=numpy.argwhere(numpy.absolute(grid.u)>tol)

class int_track(spt):
    def __init__(self,grid,tol=10**-10):
        #TODO: Add more interface tracking options.
        super(int_track,self).__init__(grid,tol)
        grid.l.info('Tracking interface')
        sset=set(map(tuple,self.s))
        if grid.u.ndim is not 2:
            raise ValueError('Interface tracker currently supports 1 spatial dimension.')
        tax=grid._get_timelike_axis()
        oax=[l for l in range(self.g.u.ndim) if l != tax]
        oax=oax[0]
        right_edge=[]
        right_edgex=[]
        for i in grid.get_ind_range(ax=tax,aset='all'):
            curr=set.intersection(sset,self.g.get_plane_ind(val=i,ax=tax,aset='all'))
            #curr=curr.intersection(sset)
            #for ind in sset:
            #   if ind[tax]==i:
            #       print ind[oax]
            try:
                ind_x_m=max(x[oax] for x in curr)
            except ValueError:
                ind_x_m=0
            #print ind_x_m
            #TODO: This hack will not work with nonuniform grids
            tmp=grid.ind_to_x(tuple([i,ind_x_m]))
            right_edge.append(tmp)
            #print tmp
            right_edgex.append(tmp[oax])
        self.rex=right_edgex
        self.re=right_edge
