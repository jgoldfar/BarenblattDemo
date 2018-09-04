"""Validates the difference scheme against known solution in
    the case where barenblatt solution applies after translation
    (i.e. linear convection case.)
"""

from pde.grids.uniform import twoddir
from pde.grids.base import trivgrid
from pde.methods.explicit import diffusion
from pde.bc import hom
from numpy import array, savetxt

from os import path
import os

def do_validate(dir_out, sigma = 5, beta = 1, b = 0, t0 = 0.5, h = 0.1, tau = pow(0.1, 2)):
    dir_with_data = path.join(dir_out, 'sigma=' + str(sigma), 'beta=' + str(beta))
    if not path.isdir(dir_with_data):
        os.makedirs(dir_with_data)
    gridp = twoddir(h=(tau, h), rx=((0, 0.2), (-4, 4)))
    gridp.set_active_int_plus_top()
    gridp.set_bc(hom.barenblatt(gridp, t0, sigma))
    soln = diffusion(gridp,b,beta,sigma)
    soln.calculate()
    savetxt(path.join(dir_with_data, 'xv.dat'), gridp.x[1], delimiter=',', fmt='%10.5f')
    savetxt(path.join(dir_with_data, 'tv.dat'), gridp.x[0], delimiter=',', fmt='%10.5f')
    savetxt(path.join(dir_with_data, 'u_approx.dat'), soln.g.u, delimiter=',', fmt='%10.5f')
    grcomp=trivgrid(array([[hom._barenblatt(t0+i, j+b*i, gridp.u.ndim-1, sigma) for j in gridp.x[1]] for i in gridp.x[0]]))
    savetxt(path.join(dir_with_data, 'u_exact.dat'), grcomp.u, delimiter=',', fmt='%10.5f')
    savetxt(path.join(dir_with_data, 'u_diff.dat'), grcomp.u-soln.g.u, delimiter=',', fmt='%10.5f')

# CLI
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("basedir", help="Base directory for result output")
parser.add_argument("sigma", help="sigma exponent for PME", type=float, default = 5)
parser.add_argument("beta", help="beta exponent for NDCE", type=float, default = 1)
parser.add_argument("--t0", help="Initial moment to start Barenblatt solution", type=float, default = 0.5)
parser.add_argument("--h", help="Space grid spacing", type=float, default=0.05)
parser.add_argument("--tau", help="Time grid spacing", type=float, default=pow(0.05, 2))

args = parser.parse_args()

do_validate(args.basedir, args.sigma, args.beta, 0, args.t0, args.h, args.tau)