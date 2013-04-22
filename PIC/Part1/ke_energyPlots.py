# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 17:45:39 2013

@author: sousae
"""

from pylab import *
from numpy import *


N = ["0128", "0512", "2048", "8192"]
cr= ["b", "r", "g", "k"]

dt= ["1.57", "0.79", "0.39"]
sh= ["-", "--", ":"]

for k in range(0,4):
    for n in range(0,3):
        fiName = "../kineticEnergy"+N[k]+"_dt"+dt[n]+".txt"
        DataIn = loadtxt(fiName)
        plot(DataIn[:,0],DataIn[:,1], color=cr[k], linestyle=sh[n], linewidth=1.0)

axis([0, 8.*pi, 0., 3100])
xlabel("time")
ylabel("Kinetic Energy")
figname = 'KE_history.png'
savefig(figname, dpi=200)