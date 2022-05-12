#Read and write files

import pandas as pd
import numpy as np

def read_out(name):
    df = pd.read_csv(name, sep=" " , header = None)
    df.columns = ['x_coords','y_coords','z_coords']
    return df

def write_out(dataframe,name):
    xx = dataframe.to_numpy()
    np.savetxt(name, x)
    
#Author      : Kalith M Ismail.
#Objective   : Access, read and write coordinate files. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
