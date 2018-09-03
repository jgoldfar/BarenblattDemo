"""
This module contains unit tests for the BC module
"""
import pde
from pde import bc

def test_np_safe_abs(NX = 100):
	xzeros = [0.0 for i in range(NX)]
	assert bc.hom.np_safe_abs(xzeros) == 0
	
	xones = [1.0 for i in range(NX)]
	assert bc.hom.np_safe_abs(xones) > 0

import numpy

def test_barenblatt(NX = 100):
	x=numpy.linspace(-2,2,NX,True)
	y = [bc.hom._barenblatt(1,xi, 1, 2) for xi in x]
	assert numpy.all([yi >= 0 for yi in y])