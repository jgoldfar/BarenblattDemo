"""
This module contains unit tests for the PDE module
"""
#from pde.grids.uniform import twoddir
#from pde.plotting import matplot, ddplot
#from pde.grids.base import trivgrid
#from pde.methods.explicit import ftcs1d as ft, lin_diffusion
#from pde.bc import hom
#from pde.interface import spt 
#TODO: Break (all) of these into unit tests to be nice.
#import matplotlib.pyplot as plt
#import numpy
#x=numpy.linspace(-2,2,1000,True)
#y1=[hom._barenblatt(1,xi,1,2,n0=1) for xi in x]
#y2=[hom._barenblatt(1,xi,1,2,n0=2) for xi in x]
#ddplot._plot(x,y1)
#ddplot._plot(x,y2)
#ddplot._show()
#grid1=twoddir((0.0005,0.1),((0,0.01),(-2,2)))
#grid1.set_active_int_plus_top()
#grid1.set_bc(hom.barenblatt(grid1, 1,2))
#grid1.set_bc(hom.dd_xpowalpha(grid1, 1.0/3,cutoff=True,negate=False,poscut=False))
#grid2=twoddir((0.005,0.1),((0,1),(0,1)))
#grid2.set_active_int_plus_top()
#grid2.set_bc(hom.test1(grid2))
#grid3=twoddir((0.005,0.1),((0,1),(0,1)))
#grid3.set_active_int_plus_top()
#grid3.set_bc(hom.test1(grid3))

#try:
#    print grid1.u[11,10]
#except IndexError:
#    print 'Of course you can\'t reference (11,10)'

#m,b,beta=2,0,1
##ux=0.5
#def ffn(x,t,u,ux,uxx):
#    #print m*(u**(m-1))*uxx+b*beta*(u**(beta-1))**ux
#    return m*pow(abs(u),m-1)*uxx+m*(m-1)*pow(abs(u),m-2)*abs(ux)+b*beta*pow(u,beta-1)*ux

#matplot._surf(grid1)
#matplot._show()

#soln=ft(grid1,ffn)
#soln.calculate()
#solnp=lin_diffusion(grid2,alpha=1)
#solnp.calculate()
#print soln.g.u-solnp.g.u
#grid.set_test_linorder()
#print grid.active
#print numpy.amax(grid.active,0)," ",numpy.amin(grid.active,0)

#sp=matplot(grid)
#matplot._surf(soln.g)
#sp=spt(soln.g)
#for ind in sp.s:
#    print ind
#matplot._spy_active(soln.g.u.shape, sp.s)
#matplot._show()
#matplot._surf(soln.g)
#matplot._surf(solnp.g)
#plt.show()
#matplot._surf(trivgrid(soln.g.u-solnp.g.u))
#matplot._show()
#matplot._spy_active(grid.u.shape, grid.active)
#soln=lin_diffusion(grid,alpha=0.5)
#print soln.get_cond()
#soln.calculate()
#matplot._surf(soln.g)
#matplot._show()

#t=timeit.Timer("""
#with open('results' + str(time.time()) + '.dat', 'w') as f:
#    gamma,sigma,beta,b=1,2,1,-5
#    grid = twoddir((0.0005, 0.1), ((0, 0.01), (-2, 2)))
#    grid.set_active_int_plus_top()
#    grid.set_bc(hom.dd_xpowalpha(grid, gamma, cutoff=True, negate=False, poscut=False))
#    def ffn(x, t, u, ux, uxx):
#        return sigma * pow(u, sigma - 1) * uxx + (sigma * (sigma - 1) * pow(u, sigma - 2) + b * beta * pow(u, beta - 1)) * ux
#    soln = ft(grid, ffn)
#    soln.calculate()
#    sp1 = tr(soln.g, 10 ** -10)
#    sp2 = tr(soln.g, 10 ** -8)
#    sp3 = tr(soln.g, 10 ** -6)
#    f.write(str(sigma) + "," + str(b) + "," + str(beta) + "," + str(gamma) + "," + str(sp1.rex[len(sp1.rex) - 1] - sp1.rex[0]) + "," + str(sp2.rex[len(sp2.rex) - 1] - sp2.rex[0]) + "," + str(sp3.rex[len(sp3.rex) - 1] - sp3.rex[0]))""","""gc.enable()
#from pde.grids.uniform import twoddir
#from pde.methods.explicit import ftcs_1d as ft
#from pde.bc import hom
#from pde.interface import int_track as tr""")
#print t.timeit(20)/20
