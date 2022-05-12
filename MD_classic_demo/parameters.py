#Pre-defined factors

import math
from input import *

dt          = Delta/tf
dt2         = dt * dt

v0          = math.sqrt(bltz_const * temp/mass)

lj_s6       = lj_s**6
lj_s12      = lj_s**12

#Author      : Kalith M Ismail.
#Objective   : Parameters estimation for simulation. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
