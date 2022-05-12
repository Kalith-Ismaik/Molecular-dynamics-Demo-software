#Estimating the interatomic potential and force

import numpy as np
import pandas as pd
import math

from input import *
from parameters import *

def lennard_jones(pos,neigh,N):
    distt2 = []
    for x in range(N):
        dist = []
        for y in range(len(neigh[x])):
            distt = pos.iloc[x,0]*pos.iloc[neigh[x][y],0] + pos.iloc[x,1]*pos.iloc[neigh[x][y],1]+pos.iloc[x,2]*pos.iloc[neigh[x][y],2]
            #distt = distt**2
            #distt = math.sqrt(distt)
            dist.append(distt)
        distt2.append(dist)    
#    print(distt2)
    U = []
    F = []
    for i in range(N):
        UU = []
        FF = []
        for j in range(len(distt2[i])): 
            dist6   = distt2[i][j]**3
            dist12  = dist6**2
            dist13  = dist12*math.sqrt(distt2[i][j])
            dist7   = dist6*math.sqrt(distt2[i][j])
    
            u  = 4*lj_e*(lj_s12/dist12 - lj_s6/dist6)                            # Potential energy calculation
            f  = 4*lj_e*(((12 * lj_s12) / dist13) - ((6 * lj_s6) / dist7))       # Force calculation
            
            UU.append(u)
            FF.append(f)
        if len(UU) > 0:    
            U_avg = sum(UU)/len(UU)
            #F_avg = sum(FF)/len(FF)
            UU    = []
            #FF    = []
            UU.append(U_avg)
            #FF.append(F_avg)
        U.append(UU)
        F.append(FF)
        
    return (U, F)
    
#Author      : Kalith M Ismail.
#Objective   : Estimating the interatomic potential and force corresponding to lennard jones condition.
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
