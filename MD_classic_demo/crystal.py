#create the experiment crystal.

import random
import pandas as pd
import numpy as np

from input import *

x_coords = []
y_coords = []
z_coords = []

for c in range(N_atom):
    x = random.random()
    y = random.random() 
    z = random.random() 
    pos = [x,y,z]
    
    for i in range(3):
        if pos[i] > (box[i]/100):
            pos[i] -= (box[i]/100)
            
        else:
            pos[i] = pos[i]
        pos[i] = round(pos[i] * (10**5)) / 1000
        
    x_coords.append(pos[0])
    y_coords.append(pos[1])
    z_coords.append(pos[2])
    
df = pd.DataFrame(list(zip(x_coords,y_coords,z_coords)),
                  columns =['x_coords','y_coords','z_coords'])
coords = df.to_numpy()
np.savetxt('Initial_conf.out', coords)
#dff = df.copy()

#Author      : Kalith M Ismail.
#Objective   : Create the experiment crystal with random arrangement of atoms. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
