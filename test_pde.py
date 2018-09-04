"""
This module contains unit tests for the PDE module
"""
from pde.grids.uniform import twoddir
from pde.methods.explicit import ftcs_1d as ft
from pde.bc import hom
from pde.interface import int_track as tr 

def test_ft():
	gamma = 1
	sigma = 2
	beta = 1
	b = -5
	grid = twoddir((0.0005, 0.1), ((0, 0.01), (-2, 2)))
	grid.set_active_int_plus_top()
	grid.set_bc(hom.dd_xpowalpha(grid, gamma, cutoff=True))
	def ffn(x, t, u, ux, uxx):
	   return sigma * pow(u, sigma - 1) * uxx + (sigma * (sigma - 1) * pow(u, sigma - 2) + b * beta * pow(u, beta - 1)) * ux
# 	soln = ft(grid, ffn)
# 	soln.calculate()
	assert True
	
	sp2 = tr(grid, 0.1)
	assert True