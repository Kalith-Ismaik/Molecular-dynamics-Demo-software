#Parameters are provided for Helium(He) atom.

box         = (50.0,50.0,50.0)   # Dimensions             [In Angstrom]
bound       = (1,1,1)
N_atom      = 50                  # No of atoms
N_step      = 500                 # No of Steps
N_write     = 100                 # Steps at which co-ordinate to be saved
Delta       = 10.0               # Time steps             [In ferrosecond]

temp        = 200.0;              # Temperature            [In Kelvin]

lj_s        = 2.5238;            # Lennard_jones_sigma    [In Angstrom]
lj_e        = 0.01962;           # Lennard_jones_epsilon  [In Kcal/mol]
mass        = 4.002602;          # Mass of atom           [In amu] 
bltz_const  = 0.001987191;       # Boltzmann constant     [In Kcal/mol/K]

tf          = 1**(-15)           # Time factor            [Convert second to ferrosecond]

#Author      : Kalith M Ismail.
#Objective   : Input conditions for simulation. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 12/05/2022.
