# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:09:46 2013

@author: sousae
"""

def createFigures(arr, comp):
    clf()
    plot(arr[:,0],arr[:,comp], linewidth=2.0)
    xlabel(r'time')
    ylabel(r'Field Energy')

from numpy import *
from pylab import *
from math import *

part_a = loadtxt('Energies0002_dt0.31.txt')
part_b = loadtxt('Energies0002_b.txt')
part_b2= loadtxt('Energies0002_b2.txt')
part_c = loadtxt('Energies0002_c.txt')

createFigures(part_a,2)
savefig("part_a.png", dpi=100)

createFigures(part_b,1)
ylabel(r'Kinetic Energy')
savefig("part_b.png", dpi=100)

createFigures(part_b2,1)
ylabel(r'Kinetic Energy')
savefig("part_b2.png", dpi=100)

createFigures(part_c,2)
savefig("part_c.png", dpi=100)

clf()
dt = pi/10.
w  = array([2.780,2.618,2.474,2.212,2.020,1.939,1.424,1.025])
wpe= array([5.917,5.781,5.642,5.276,4.886,4.514,3.192,2.257])

yy = w*dt/2.
xx = wpe*dt/2.
th = zeros(size(xx))
for n in range(0,size(xx)):
	th[n] = asin(xx[n])
	
plot(xx,yy, 'ob', label='data', linewidth=2.0)
plot(xx,th, '-m', label='analytical', linewidth=1.0)
xlabel(r'$\omega_{pe}\Delta t /2$')
ylabel(r'$\omega\Delta t/2$')
savefig("part_c2.png",dpi=100)
