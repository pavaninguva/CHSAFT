# Importing relevant packages 
# random is needed for noise generation in initialisation of system
# dolfin is the finite element solver package within the fenics environment
# numpy is for general operations
import random
from dolfin import *
import numpy as np

# Defining constants for system
#chi_ab is the flory-huggins binary interaction parameter
chi_ab = 0.006
#N_A / N_B are the polymer chain lengths in terms of monomer units
N_A = 1000
N_B = 1000
# D_AB is the diffusion coefficient for the polymeric species
D_AB = 1e-11 
