#Velocity verlet algorithm
def new_pos(X_old,V_old,F_old,m,t):                                    #Algorithm For new position
    X_new = X_old + (t)*V_old + ((t**2)*F_old)/(2*m)
    return(X_new)

def new_vel(F_old,F_new,V_old,t,m):                                    #Algorithm for new velocity
    V_new = (F_old + F_new)
    V_new = (V_new * t)/(2*m)
    V_new = V_old + V_new
    return(V_new)
    
#Author      : Kalith M Ismail.
#Objective   : Velocity verlet algorithm.  
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
