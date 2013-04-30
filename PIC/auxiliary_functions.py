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
    
def CreatePosVelPlot(pos, vel, E, r, n, xl, xu, vl, vu):
    subplot(2,2,1)
    plot(pos, vel, '.b'), axis([xl,xu,vl,vu])
    xlabel("Position")
    ylabel("Velocity")
    subplot(2,2,2)
    plot(r,'-*')
    subplot(2,2,4)
    plot(E,'-o')
    figname = "tracking" +str('%03d' % n) + '_.png'
    savefig(figname, dpi=100)

def EfieldEnergy(E):
    
    N = size(E)
    EE = 0.0
    for k in range(0,N-1):
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
        
        self.connect[0,0] = npoints-1
        
def EfieldSolve(rho, grid):
    ii = complex(0,1)
#    Efield = zeros(grid.gPoints)
#    kx = linspace(-grid.gPoints/2,grid.gPoints/2-1,grid.gPoints)
    kx = linspace(0,grid.gPoints-1,grid.gPoints)
    kx = 2.*pi/(grid.xupper-grid.xlower)*kx
    kx[kx==0] = 1.e-6
    
    rhofft = fftpack.fft(rho)
    phi  = -rhofft/(kx*kx)
    Efft = -ii*kx*phi
    
    E_x = fftpack.ifft(Efft)
    
#    Efield[0:grid.gPoints-1] = E_x.real
#    Efield[grid.gPoints-1] = Efield[0]
    
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
        #print el
        dens[grid.connect[el,0]] += abs(grid.grid[grid.connect[el,1]]-pos[k])/grid.dx
        dens[grid.connect[el,1]] += abs(pos[k]-grid.grid[grid.connect[el,0]])/grid.dx
        
    #dens[0] += dens[grid.gPoints-1]
    #dens[grid.gPoints-1] = dens[0]
        
    return dens
        
def firstOrderField(pos, Efield, grid):
    
    N = size(pos)
    
    E_part = zeros(N)
    
    for k in range(0,N):
        el = int((pos[k]+grid.xupper)/grid.dx)
        E_part[k] = (grid.grid[grid.connect[el,1]]-pos[k])/grid.dx * \
                    Efield[grid.connect[el,0]] + \
                    (pos[k]-grid.grid[grid.connect[el,0]])/grid.dx * \
                    Efield[grid.connect[el,1]]
                    
    return E_part
    
def SinosoidalVel(Amp, pos, phase):
    N = size(pos)
    v = zeros(N)
    
    for k in range(0,N):
        v[k] = Amp*sin(pos[k]+phase)
        
    return v