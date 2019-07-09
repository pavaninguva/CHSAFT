# Importing relevant packages 
# random is needed for noise generation in initialisation of system
# dolfin is the finite element solver package within the fenics environment
# numpy is for general operations
import random
from dolfin import *
import numpy as np

#chi_ab is the flory-huggins binary interaction parameter
chi_ab = 0.006
#N_A / N_B are the polymer chain lengths in terms of monomer units
N_A = 1000
N_B = 1000
# D_AB is the diffusion coefficient for the polymeric species
D_AB = 1e-11 

# Intial mole fraction of species A
A_RAW = 0.5

#Numerics 
DT = 0.05
TIME_MAX = 20
N_CELLS = 80
DOMAIN_LENGTH = 40
theta_ch = 0.5


class InitialConditions(Expression):
    def __init__(self, **kwargs):
        random.seed(1234)

    def eval(self, values, x):
        values[0] = A_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[2] = 0.0

    def value_shape(self):
        return (2,)

class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity)
        self.reset_sparsity = False

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"

mesh = RectangleMesh(
    Point(0.0, 0.0), Point(DOMAIN_LENGTH, DOMAIN_LENGTH), N_CELLS, N_CELLS
)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)

CH = FunctionSpace(mesh, MixedElement([P1, P1, P1, P1, P1]))

dch = TrialFunction(CH)
h_1, j_1 = TestFunctions(CH)

ch = Function(CH)
ch0 = Function(CH)

# Split mixed functions
da, dmu_AB = split(dch)
a, mu_AB = split(ch)
a0, mu0_AB = split(ch0)

a = variable(a)

ch_init = InitialConditions(degree=1)
ch.interpolate(ch_init)
ch0.interpolate(ch_init)

kappa = (2.0/3.0)*chi_AB

#Flory-Huggins Expression
g = ( a * ln(a) / N_A ) + ((1.0-a)*ln(1-a)/ N_B) + a*(1.0-a)*chi_AB

# Using the fenics autodifferentiation toolkit 
dgda = diff(g,a)

# Introduce an expression for mu_{n+theta}
mu_AB_mid = (1.0 - theta_ch) * mu0_AB + theta_ch * mu_AB

dt = DT

F_a = (
    a * h_1 * dx
    - a0 * h_1 * dx
    + dt * a * (1.0 - a) * D_AB_ * dot(grad(mu_AB_mid), grad(h_1)) * dx
)

F_mu_AB = (
    mu_AB * j_1 * dx
    - dgda * j_1 * dx
    - kappa * dot(grad(a), grad(j_1)) * dx
)

F = F_a + F_mu_AB

#Compute directional derivative about u in the direction of du
a = derivative(F, ch, dch)

problem = CahnHilliardEquation(a, F)

solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6

file_a = File("concentration_A.pvd", "compressed")

t = 0.0
timestep = 0

# Output intial conditions
file_a << (ch.split()[0], t)

space = FunctionSpace(mesh, P1)

while t < TIME_MAX:


    timestep += 1
    t += dt

    if MPI.rank(mpi_comm_world()) == 0:

        
        print "Timestep", timestep, "Time", t

    proj_a = project(ch.split()[0], FunctionSpace(mesh, P1))

    gather_a = Vector()

    proj_a.vector().gather(gather_a, np.array(range(space.dim()), "intc"))

    if MPI.rank(mpi_comm_world()) == 0:


        print gather_a.array().shape
        print gather_a.array().min()
        print gather_a.array().max()