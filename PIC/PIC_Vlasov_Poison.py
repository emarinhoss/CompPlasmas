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
N = 64

# Particle mass and charge
m = 1.0
q = 1.0

# Full width at half maximum
fwhm = 2.

# final time
Tend = 20.*pi
nsteps = 200

# Order of Particle weighting
order = 1

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
#stdev = fwhm/(2.*sqrt(2.*log(2.)))
#pos, vel = Initialization(N, grid, stdev)

dum = linspace(xlower,xupper,N+2)
pos = dum[1:1+N]
#pos1 = linspace(xlower,xupper,N)
#pos = concatenate((pos1,pos1))
#pos = array([0.0])
#pos1.append(pos1)

neg = dum[1:1+N]
#neg = linspace(xlower,xupper,N)
#speed = 0.1
#delta = 0.001
#vel1 = speed*ones(N) + SinosoidalVel(delta, pos1, 0.0)
#vel2 = -speed*ones(N) + SinosoidalVel(delta, pos1, pi)
#vel = concatenate((vel1,vel2))
#vel = array([-0. , 0.])
vel = SinosoidalVel(0.005, pos, 0.0)

# Plot histogram of velocity distribution
#VelocityHistogram(vel, N, Tend)

## timesteps
dt = Tend/nsteps
time = linspace(0,Tend,nsteps)
KE = zeros(nsteps)
EE = zeros(nsteps)

## Plot Initial conditions
CreatePosVelPlot(pos, vel , zeros(Ngrid), zeros(Ngrid), 0, grid.xlower, grid.xupper, vlower, vupper)
clf()

if order == 1:
    rho_n = firstOrderParticle(neg, grid)
else:
    rho_n = zerothOrderParticle(neg, grid)

# Advance solution
for n in range(0, nsteps):
    
    # Particle Weighting
    if order == 1:
        rho_p = firstOrderParticle(pos, grid)
    else:
        rho_p = zerothOrderParticle(pos, grid)
        
    rho = -rho_p + rho_n

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
    CreatePosVelPlot(pos, vel, rho, Exj, n+1, grid.xlower, grid.xupper, vlower, vupper)
    clf()

DataOut = column_stack((time,KE,EE))
KEname = "Energies" + str('%04d' % N) + "_dt" + str('%2.2f' % dt) + ".txt"
savetxt(KEname, DataOut)