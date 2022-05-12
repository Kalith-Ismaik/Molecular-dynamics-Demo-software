# MAIN SIMULATION LOOP

import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from input import *
from boundary import *
from cutoff import *
from files import *
from LJ import *
from neighbours import *
from parameters import *
from vel import *
from VV_algo import *

df = read_out("Initial_conf.out")

co = cutoff()

for step in range(N_step + 1):
    
    xx_coords = []
    yy_coords = []
    zz_coords = []
    
    # Zero the forces
    for i in range(N_atom):
        for c in range(3):
            force[i][c] = 0.0
           
    # Calculate interaction forces
    if step%10 == 0:
        sk = Skin(df,N_atom,co)
    
    neig     = Neighbours(df,sk,N_atom,co)
    
    U_E, I_F =lennard_jones(df,neig,N_atom)
    
    ener_pot = 0.0
    for kk in range(N_atom):
        if len(U_E[kk]) > 0:
            ener_pot += U_E[kk][0]
    
    force_new = force
    for i in range(N_atom):
        if len(I_F[i]) > 0:
            for j in range(len(I_F[i])):
                for k in range(3):
                    force_new[i][k]          += I_F[i][j]  
                    force_new[neig[i][j]][k] -= I_F[i][j]
    force = force_new
    
    pos_new = []
    vel_new = []
    for i in range(N_atom):
        p = []
        v = []
        for j in range(3):
            ijk = new_pos(df.iloc[i,j],vel0[i][j],force_new[i][j],mass,dt)
            vvv = new_vel(force[i][j],force_new[i][j],vel0[i][j],dt,mass)
            p.append(ijk)
            v.append(vvv)
        vel_new.append(v)
        pos_new.append(p)
            
        xx_coords.append(p[0])
        yy_coords.append(p[1])
        zz_coords.append(p[2])
        
    # Calculate the kinetic energy
    ener_kin = 0.0
    for i in range(N_atom):
        ener_kin += 0.5*mass*(vel_new[i][0]**2 + vel_new[i][1]**2 + vel_new[i][2]**2)
    ener_total = ener_pot + ener_kin
    mean_temp = 2.0*ener_kin/(3*bltz_const*(N_atom-1))
    print("step %9d  ener_total %9.4f  ener_pot %9.4f  ener_kin %9.4f mean_temp %8.3f" % (step,ener_total,ener_pot,ener_kin,mean_temp))
    
    n_df = pd.DataFrame(list(zip(xx_coords,yy_coords,zz_coords)),
                    columns =['xx_coords','yy_coords','zz_coords'])
    periodic_cond(n_df,box,N_atom,bound)
    df = n_df
    vel0 = vel_new
    if step%N_write == 0:
        crds = df.to_numpy()
        name = str(step) + '.out'
        np.savetxt(name, crds)
        
#Author      : Kalith M Ismail.
#Objective   : Molecular Dynamics simulation of noble gas atoms in the box. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
