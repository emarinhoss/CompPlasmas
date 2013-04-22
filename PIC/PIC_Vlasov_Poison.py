# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 18:41:51 2013

@author: sousae
"""
from math import *
from auxiliary_functions import *
from pylab import *
from numpy import *


####################################################################
###
###
###
###
###
####################################################################

# x domain limits and gridpoints
xlower = -pi
xupper = pi
Ngrid  = 32

# v_x limits
vlower = -5.
vupper =  5.

# number of particles
N = 2

# Particle mass and charge
m = 1.0
q = 1.0

# Full width at half maximum
fwhm = 2.

# final time
Tend = 16.*pi
nsteps = 40

# Order of Particle weighting
order = 0


####################################################################
###
###
###
###
###
####################################################################
# Create a grid
grid = createGrid(xlower, xupper, Ngrid)

# Initialize random position and velocity
stdev = fwhm/(2.*sqrt(2.*log(2.)))
#pos, vel = Initialization(N, grid, stdev)

pos = array([-pi/4., pi/4.])
dum = linspace(xlower,xupper,N+2)
neg = dum[1:1+N]
vel = array([0., 0.])

# Plot histogram of velocity distribution
#VelocityHistogram(vel, N, Tend)

## timesteps
dt = Tend/nsteps
time = linspace(0,Tend,nsteps)
KE = zeros(nsteps)
EE = zeros(nsteps)

## Plot Initial conditions
CreatePosVelPlot(pos, vel, 0, grid.xlower, grid.xupper, vlower, vupper)
clf()

## Advance solution
for n in range(0, nsteps):
    
    # Particle Weighting
    if order == 1:
        rho_p = firstOrderParticle(pos, grid)
        rho_n = firstOrderParticle(neg, grid)
    else:
        rho_p = zerothOrderParticle(pos, grid)
        rho_n = zerothOrderParticle(neg, grid)
        
    rho = rho_n-rho_p

    # Field Solve
    Exj = EfieldSolve(rho, grid)
    
    # Field weighting
    if order == 1:
        Exi = firstOrderField(pos, Exj, grid)
    else:
        Exi = zerothOrderField(pos, Exj, grid)
        
    #Particle push
    vel = vel + Exi*dt
    pos = pos + vel*dt
    # Apply BC's
    PeriodicBC(pos, grid)
    
    # calculate the energies in the system
    en = KineticEnergy(vel)
    ee = EfieldEnergy(Exj)
    KE[n] = m*en
    EE[n] = grid.dx*ee
    
    # Plot current solution
    CreatePosVelPlot(pos, vel, n+1, grid.xlower, grid.xupper, vlower, vupper)
    clf()

DataOut = column_stack((time,KE,EE))
KEname = "Energies" + str('%04d' % N) + "_dt" + str('%2.2f' % dt) + ".txt"
savetxt(KEname, DataOut)