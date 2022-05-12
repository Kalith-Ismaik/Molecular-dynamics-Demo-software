#Initial velocity and momentum

import random
from input import *
from parameters import *

vel0 = []
momentum = [0.0,0.0,0.0]
for k in range(N_atom):
    vel0.append([random.gauss(0,v0), random.gauss(0,v0), random.gauss(0,v0)])
    for kk in range(3):
        momentum[kk] += mass*vel0[k][kk] 

#Removal of residual momentum        
residual_mom = []        
for c in range(3):
    residual_mom.append(momentum[c]/N_atom)
for i in range(N_atom):
    for c in range(3):
        vel0[i][c] -= residual_mom[c]/mass

force = []
for i in range(N_atom):
    force.append([0.0, 0.0, 0.0])
    
#Author      : Kalith M Ismail.
#Objective   : Conditioning atoms in the crystal with initial velocity and momentum correspondint to temperature. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
