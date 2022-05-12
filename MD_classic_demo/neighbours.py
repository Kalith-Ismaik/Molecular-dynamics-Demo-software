#Neighbouring atoms estimation

import  numpy as np
import pandas as pd

def Skin(pos,N,cutoff):
    r = cutoff + 4.0
    neigh = []
    for x in range(N):
        nh = []
        for y in range(x+1,N):
            if pos.iloc[x,0]-pos.iloc[y,0]<=r and pos.iloc[x,1]-pos.iloc[y,1]<=r and pos.iloc[x,2]-pos.iloc[y,2]<=r:
                nh.append(y)
        for y in range(x-1,0,-1):
            if pos.iloc[x,0]-pos.iloc[y,0]<=r and pos.iloc[x,1]-pos.iloc[y,1]<=r and pos.iloc[x,2]-pos.iloc[y,2]<=r:
                nh.append(y)
        neigh.append(nh)
    return neigh
    
def Neighbours(pos,skin,N,cutoff):
    r = cutoff + 2.0
    neighh = []
    for x in range(N):
        nhh = []
        for y in range(len(skin[x])):
            if pos.iloc[x,0]-pos.iloc[skin[x][y],0]<=r and pos.iloc[x,1]-pos.iloc[skin[x][y],1]<=r and pos.iloc[x,2]-pos.iloc[skin[x][y],2]<=r:
                nhh.append(skin[x][y])
        neighh.append(nhh)
    return neighh
    
#Author      : Kalith M Ismail.
#Objective   : Estimating the neighbouring atoms corresponding to their cutoff distance. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
