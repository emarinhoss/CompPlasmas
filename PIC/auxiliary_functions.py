# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 19:24:01 2013

@author: sousae
"""
from numpy import *
from random import *
from pylab import *
from scipy import sparse
from scipy import fftpack

def Initialization(N, grid, stdev):
    
    vel = zeros(N)
    pos = zeros(N)
    
    for n in range(0,N):
        vel[n] = gauss(0.0, stdev)
        pos[n] = uniform(grid.xlower, grid.xupper)
        
    return pos, vel
    
    
def PeriodicBC(pos, grid):
    N = size(pos)
    d = grid.xupper-grid.xlower
    
    for k in range(0,N):
        
        if(pos[k]<grid.xlower):
            pos[k] = pos[k] + d
        elif(pos[k]>grid.xupper):
            pos[k] = pos[k] - d
            
def KineticEnergy(vel):
    N = size(vel)
    KE = 0.0
    
    for k in range(0, N):
        KE += 0.5*vel[k]*vel[k]
    
    return KE
    
def VelocityHistogram(vx, N, Tend):
    clf()
    hist(vx, bins=40)
    xlabel("Velocity")
    ylabel("Number of particles")
    histname = "histogram_N" + str('%04d' % N) + "_T" + str('%2.2f' % Tend) + ".png"
    savefig(histname, dpi=100)
    clf()
    
def CreatePosVelPlot(pos, vel, n, xl, xu, vl, vu):
    plot(pos, vel, '.b'), axis([xl,xu,vl,vu])
    xlabel("Position")
    ylabel("Velocity")
    figname = "tracking" +str('%03d' % n) + '_.png'
    savefig(figname, dpi=100)

def EfieldEnergy(E):
    
    N = size(E)
    EE = 0.0
    for k in range(0,N):
        EE += E[k]*E[k]
        
    return EE
    
class createGrid:
    
    def __init__(self, xl, xu, npoints):
        self.xlower = xl
        self.xupper = xu
        self.gPoints = npoints
        self.grid = linspace(xl,xu,npoints)
        self.dx = self.grid[1]-self.grid[0]
        self.connect = zeros((npoints-1,2))
        
        for k in range(0,npoints-1):
            for n in range(0,2):
                self.connect[k,n] = k+n
        
def EfieldSolve(rho, grid):
    ii = complex(0,1)
    kx = linspace(0,grid.gPoints-1,grid.gPoints)
    #kx = zeros(grid.gPoints)
    #kx[0:grid.gPoints/2-1] = range(0,grid.gPoints/2-1)
    #kx[grid.gPoints/2:grid.gPoints] = range(-grid.gPoints/2,0)
    kx = 2.*pi/grid.xupper*kx
    kx[0] = 1.e-6
    
    rhofft = fftpack.fft(rho)
    phi  = rhofft/(kx*kx)
    Efft = ii*kx*phi
    
    E_x = fftpack.ifft(Efft)
    
    return E_x.real
    
def zerothOrderParticle(pos, grid):
    
    N = size(pos)
    
    dens = zeros(grid.gPoints)
    
    for k in range(0,N):
        el = int((pos[k]+grid.xupper)/grid.dx)
        
        if(pos[k] <= grid.grid[grid.connect[el,0]]+0.5*grid.dx):
            dens[grid.connect[el,0]] += 1.
        else:
            dens[grid.connect[el,1]] += 1.
            
    dens[0] += dens[grid.gPoints-1]
    dens[grid.gPoints-1] = dens[0]
    
    return dens
            
def zerothOrderField(pos, Efield, grid):
    
    N = size(pos)
    
    E_part = zeros(N)
    
    for k in range(0,N):
        el = int((pos[k]+grid.xupper)/grid.dx)
        
        if(pos[k] <= grid.grid[grid.connect[el,0]]+0.5*grid.dx):
            E_part[k] = Efield[grid.connect[el,0]]
        else:
            E_part[k] = Efield[grid.connect[el,0]]
                        
    return E_part
            
def firstOrderParticle(pos, grid):
    
    N = size(pos)
    
    dens = zeros(grid.gPoints)
    
    for k in range(0,N):
        el = int((pos[k]+grid.xupper)/grid.dx)
        
        dens[grid.connect[el,0]] += (grid.grid[grid.connect[el,1]]-pos[k])/grid.dx
        dens[grid.connect[el,1]] += (pos[k]-grid.grid[grid.connect[el,0]])/grid.dx
        
def firstOrderField(pos, Efield, grid):
    
    N = size(pos)
    
    E_part = zeros(N)
    
    for k in range(0,N):
        el = int((pos[k]+grid.xupper)/grid.dx)
        E_part[k] = (grid.grid[grid.connect[el,1]]-pos[k])/grid.dx * \
                    Efield[grid.connect[el,0]] + \
                    (pos[k]-grid.grid[grid.connect[el,0]])/grid.dx * \
                    Efield[grid.connect[el,1]]