"""Completes a sweep of the parameter space 
"""
#TODO: Looks much simpler to use numpy.mgrid to generate parameter combinations
from pde.grids.uniform import twoddir
from pde.methods.explicit import diffusion
from pde.bc import hom
from pde.interface import int_track as tr
from time import time
from numpy import linspace, linalg

from os import path
import os

def do_sweep(basedir, sigma, nbeta, ngamma):
	dirWithRun = path.join(basedir, 'run'+str(sigma))
	if not path.isdir(dirWithRun):
		os.makedirs(dirWithRun)
	with open(path.join(dirWithRun, 'results' + str(time()) + '.dat'), 'w') as f:
		f.write("sigma,b,beta,gamma,diff_10e-2,diff_10e-3\n")
		for beta in linspace(0,1,nbeta):
			for b in [0.01]:
				for gamma in linspace(0.9,1.8,ngamma):
					grid = twoddir((0.004, 0.1), ((0, 0.01), (-2, 2)))
					grid.set_active_int_plus_top()
					grid.set_bc(hom.dd_xpowalpha(grid, gamma, cutoff=True))
					soln = diffusion(grid,b,beta,sigma)
					try:
						soln.calculate()
					except linalg.LinAlgError:
						f.write(str(sigma) + "," + str(b) + "," + str(beta) + "," + str(gamma) + "," + str(-1000) + "," + str(-1000) + "\n")
						soln.l.warning('Unhandled linear algebra error.')
						continue
					sp2 = tr(soln.g, 0.1)
					sp3 = tr(soln.g, 0.01)
					f.write(str(sigma) + "," + str(b) + "," + str(beta) + "," + str(gamma) + "," + str(sp2.rex[len(sp2.rex) - 1] - sp2.rex[0]) + "," + str(sp3.rex[len(sp3.rex) - 1] - sp3.rex[0]) + "\n")

# CLI
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("basedir", help="Base directory for result output")
parser.add_argument("sigma", help="sigma exponent for PME", type=float)
parser.add_argument("--nbeta", help="Number of beta values to test", type=int, default=2)
parser.add_argument("--ngamma", help="Number of gamma values to test", type=int, default=2)

args = parser.parse_args()

do_sweep(args.basedir, args.sigma, args.nbeta, args.ngamma)