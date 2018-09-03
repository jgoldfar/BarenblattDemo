"""
This module contains unit tests for the Grids module
"""
from pde.grids import uniform
from pde.bc import hom

def test_twodir_barenblatt():
	grid=uniform.twoddir((0.0005,0.1),((0,0.01),(-2,2)))
	grid.set_active_int_plus_top()
	grid.set_bc(hom.barenblatt(grid, 1, 2))

def test_twodir_xpowalpha():
	grid=uniform.twoddir((0.005,0.1),((0,1),(0,1)))
	grid.set_active_int_plus_top()
	grid.set_bc(hom.dd_xpowalpha(grid, 1.0/3, cutoff=True))
	assert True

# 
import numpy

def test_twodir_set_linorder():
	grid = uniform.twoddir((0.01, 0.1), ((0,1), (0,1)))
	l=map(tuple,
			numpy.transpose(numpy.indices(
										grid.u.shape
										)
							).reshape((numpy.prod(grid.nx),2)))
	j=0
	for i in l:
		grid.u[i]=j
		j=j+1
	assert True

def test_twodir_set_xprod():
    grid = uniform.twoddir((0.01, 0.1), ((0,1), (0,1)))
    l=map(tuple,
        numpy.transpose(
                    numpy.indices(grid.u.shape)
                    ).reshape((numpy.prod(grid.nx),2)))
    for i in l:
        grid.u[i]=numpy.prod(grid.ind_to_x(i))
    assert True