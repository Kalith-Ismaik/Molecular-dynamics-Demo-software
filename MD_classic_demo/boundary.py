#Periodic Boundary Condition

import pandas as pd

def periodic_cond(pos,dim,N,bound):
    
    for x in range(N):        
        for y in range(3):
            if bound[y] == 1:
                if pos.iloc[x,y] > dim[y]:
                    pos.iloc[x,y] = -(dim[y]) + pos.iloc[x,y]
                elif pos.iloc[x,y] < 0:
                    pos.iloc[x,y] = (dim[y]) + pos.iloc[x,y]
                else:
                    pos.iloc[x,y] = pos.iloc[x,y]
            else:
                pos.iloc[x,y] = pos.iloc[x,y]
    return pos
    
#Author      : Kalith M Ismail.
#Objective   : Modification positive periodic boundary condition for MD simulation. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
