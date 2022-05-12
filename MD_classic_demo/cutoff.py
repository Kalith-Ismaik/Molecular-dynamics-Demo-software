#Interatomic potential calculation

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from input import *
from parameters import *

def cutoff():
    U    = []
    X_np = np.arange(0.001, 10, 0.001)
    X    = list(X_np)
    F    = []
    for x in X:
        dist2   = x*x
        dist6   = dist2**3
        dist12  = dist6**2
        dist13  = dist12*x
        dist7   = dist6*x
    
        u  = 4*lj_e*(lj_s12/dist12 - lj_s6/dist6)                            # Potential energy calculation
        f  = 4*lj_e*(((12 * lj_s12) / dist13) - ((6 * lj_s6) / dist7))       # Force calculation 
        
        U.append(u)
        F.append(f)
        
    min_index = U.index(min(U))
    print("Minima of interatomic potential with He-He:", U[min_index])
    print("Equilibrium interatomic distance of He-He:", X[min_index])
    
    PC = pd.DataFrame(list(zip(X,U,F)),
                      columns =['Int_dist','E_Pot','Force'])
    coords = PC.to_numpy()
    np.savetxt('LJ_Curve.out', coords)
    
    fig, ax = plt.subplots()
    ax.plot(X[min_index - 800 : min_index + 500] , U[min_index - 800 : min_index + 500],label = "E_pot")
    ax.plot(X[min_index - 800 : min_index + 500] , F[min_index - 800 : min_index + 500],label = "Force")


    ax.set(xlabel='Interatomic_Distance(r)', ylabel='Force(F) & E_pot(U)',
           title='Lennard-jones potential curve')
    ax.legend()
    
    fig.savefig("LJ_pot.png")
    plt.show()
    
    return X[min_index]
    
#Author      : Kalith M Ismail.
#Objective   : Estimate the equilibrium interatomic distance of atoms. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
