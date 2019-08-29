# Importing relevant packages 
# random is needed for noise generation in initialisation of system
# dolfin is the finite element solver package within the fenics environment
import random
from dolfin import *
# from cbcpost import *



#chi_ab is the flory-huggins binary interaction parameter
chi_AB = 6
#N_A / N_B are the polymer chain lengths in terms of monomer units
N_A = 1000
N_B = 1000
# D_AB is the diffusion coefficient for the polymeric species
D_AB = 1e-11 
# Intial mole fraction of species A
A_RAW = 0.5

#Numerics 
DT = 0.1
TIME_MAX = 20
N_CELLS = 80
DOMAIN_LENGTH = 40
theta_ch = 0.5
NOISE_MAGNITUDE = 0.03
N=int(N_CELLS)
MESH_TYPE = "structured"

class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        # self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A)#, reset_sparsity=self.reset_sparsity)
        # self.reset_sparsity = False

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True

# INITIAL AND BOUNDARY CONDITIONS OF THE PROBLEM #

# Initial conditions which include random noise as a perturbation


class InitialConditions(UserExpression):
    def __init__(self, **kwargs):
        random.seed(2 + MPI.rank(MPI.comm_world))
        super().__init__(**kwargs)
    def eval(self, values, x):
        # [0] corresponds to the concentration field for species A 
        # [1] coresponds to the \mu_{AB} field
        values[0] = A_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[1] = 0.0
    def value_shape(self):
        return (2,)

# Setting up Mesh and Periodic Boundary Conditions #

## Setting up periodic boundary conditions. 

## We are creating classes for defining the various parts of the boundaries (top, bottom, left and right)

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)


class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], DOMAIN_LENGTH)


class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0)


class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], DOMAIN_LENGTH)


class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0.0))

    # Map RightBoundary to LeftBoundary
    def map(self, x, y):
        y[0] = x[0] - DOMAIN_LENGTH
        y[1] = x[1]


N = int(N_CELLS)

if MESH_TYPE == "structured":
    mesh = RectangleMesh(
        Point(0.0, 0.0), Point(DOMAIN_LENGTH, DOMAIN_LENGTH), N, N
    )
else:
    domain_vertices = [
        Point(0.0, 0.0),
        Point(DOMAIN_LENGTH, 0.0),
        Point(DOMAIN_LENGTH, DOMAIN_LENGTH),
        Point(0.0, DOMAIN_LENGTH),
    ]

    domain = Polygon(domain_vertices)
    mesh = generate_mesh(domain, N)


# CG stands for continuous galerkin can is a lagrange type element. 
# The "1" corresponds to the order. So this is a linear lagrange element. 
# CH = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())
# P1 = FiniteElement("Lagrange", interval, 1)

# # # This function has been deprecated for 2016.1.0, but still works. 
# # # The requirement is to created a combined concentration and chemical potential function space
# # # Since we are solving for two variables. 
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
CH = FunctionSpace(mesh, P1*P1)

dch = TrialFunction(CH)
h_1, j_1 = TestFunctions(CH)

ch = Function(CH)
ch0 = Function(CH)

# Split mixed functions
da, dmu_AB = split(dch)
x_a, mu_AB = split (ch)
x_a0, mu0_AB = split(ch0)

x_a = variable(x_a)

ch_init = InitialConditions(degree=1)
ch.interpolate(ch_init)
ch0.interpolate(ch_init)

kappa = (2.0/3.0)*chi_AB

#Flory-Huggins Expression
g = ( x_a * ln(x_a) / N_A ) + ((1.0-x_a)*ln(1-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB #TO DO: Shift this to a master script

# Using the fenics autodifferentiation toolkit 
dgdx_a = diff(g,x_a)

mu_AB_mid = (1.0 - theta_ch) * mu0_AB + theta_ch * mu_AB

dt = DT

F_a = (
    x_a * h_1 * dx
    - x_a0 * h_1 * dx
    + dt * x_a * (1.0 - x_a) * dot(grad(mu_AB_mid), grad(h_1)) * dx
)

F_mu_AB = (
    mu_AB * j_1 * dx
    - dgdx_a * j_1 * dx
    - kappa * dot(grad(x_a), grad(j_1)) * dx
)

F = F_a + F_mu_AB

#Compute directional derivative about u in the direction of du
a = derivative(F, ch, dch)

problem = CahnHilliardEquation(a, F)
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "residual"
solver.parameters["relative_tolerance"] = 1e-6

# Setting up post processor

# casedir = "output_RESULTS" 
# postprocessor = PostProcessor({"casedir": casedir})


# postprocessor.add_field(SolutionField("Concentration_A", dict(save=True,
#                                 save_as=["hdf5", "xdmf"],
#                                 plot=False,
#                                 stride_timestep=20
#                                 )))

# Setting up time stepping and solving

# t = 0.0
# timestep = 0

# while (t < TIME_MAX):
#     print (t)
#     t += dt
#     timestep += 1
#     ch0.vector()[:] = ch.vector()
#     solver.solve(problem, ch.vector())
#     # It seems this is the line to link up the relevant field to the postprocesser field
#     postprocessor.update_all ({"Concentration_A"= lambda: ch.split()[0]}, t, timestep)
    
# postprocessor.finalize_all()

### IF ALL GOES WRONG JUST USE THE DAMN PVD FILE OUTPUT
# file = XDMFFile("output.xdmf")

# # Step in time
# t = 0.0
# while (t < TIME_MAX):
#     t += dt
#     ch0.vector()[:] = ch.vector()
#     solver.solve(problem, ch.vector())
#     file.write (ch.split()[0], t)
file = File("output.pvd", "compressed")

t = 0.0
while (t < TIME_MAX):
    t += dt
    ch0.vector()[:] = ch.vector()
    solver.solve(problem, ch.vector())
    file << (ch.split()[0], t)

