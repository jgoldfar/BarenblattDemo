"""Completes a sweep of the parameter space 
"""
#TODO: Looks much simpler to use numpy.mgrid to generate parameter combinations
from pde.grids.uniform import twoddir
from pde.methods.explicit import diffusion
from pde.bc import hom
from pde.interface import int_track as tr
from numpy import linspace, linalg

from os import path
import os

def do_sweep(basedir, sigma, nbeta, ngamma, dx = 1e-2, dt = 4e-3):
	dirWithRun = path.join(basedir, 'run'+str(sigma))
	if not path.isdir(dirWithRun):
		os.makedirs(dirWithRun)
	b = 0.01
	threshold1 = 1e-1
	threshold2 = 1e-2
	with open(path.join(dirWithRun, 'results.dat'), 'w') as f:
		f.write('sigma,b,beta,gamma,diff_' + str(threshold1) + ',diff_' + str(threshold2) + '\n')
		for beta in linspace(0,1,nbeta):
			for gamma in linspace(0.9,1.8,ngamma):
				grid = twoddir(h=(dt, dx), rx = [(0, 0.01), (-2, 2)])
				for (j, xj) in enumerate(grid.x[1]):
					grid.u[0, j] = hom._negxpowalpha(xj, gamma)
				
				soln = diffusion(grid,b,beta,sigma)
				
				gridout = soln.calculate()
				sp2 = tr(soln.g, threshold1)
				sp3 = tr(soln.g, threshold2)
				f.write(str(sigma) + "," + str(b) + "," + str(beta) + "," + str(gamma) + "," + str(sp2.rex[-1] - sp2.rex[0]) + "," + str(sp3.rex[-1] - sp3.rex[0]) + "\n")

# CLI
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("basedir", help="Base directory for result output")
parser.add_argument("sigma", help="sigma exponent for PME", type=float)
parser.add_argument("--nbeta", help="Number of beta values to test", type=int, default=2)
parser.add_argument("--ngamma", help="Number of gamma values to test", type=int, default=2)
parser.add_argument("--dx", help="Space grid spacing", type=float, default=1e-2)
parser.add_argument("--tau", help="Time grid spacing", type=float, default=4e-3)

args = parser.parse_args()

do_sweep(args.basedir, args.sigma, args.nbeta, args.ngamma, args.dx, args.tau)