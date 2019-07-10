# Importing relevant packages 
# random is needed for noise generation in initialisation of system
# dolfin is the finite element solver package within the fenics environment
# numpy is for general operations
import random
from dolfin import *
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


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
DT = 0.05
TIME_MAX = 20
N_CELLS = 80
DOMAIN_LENGTH = 40
theta_ch = 0.5
NOISE_MAGNITUDE = 0.03
N=int(N_CELLS)


class InitialConditions(Expression):
    def __init__(self, **kwargs):
        random.seed(1234)

    def eval(self, values, x):
        values[0] = A_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[1] = 0.0

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
        assemble(self.a, tensor=A)#, reset_sparsity=self.reset_sparsity)
        self.reset_sparsity = False

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
# parameters["form_compiler"]["representation"] = "quadrature"

# This can be used to generate a structured mesh
mesh = RectangleMesh(
    Point(0.0, 0.0), Point(DOMAIN_LENGTH, DOMAIN_LENGTH), N_CELLS, N_CELLS
)
# This is used to generated an unstructured mesh
# mesh = Mesh()

# domain_vertices = [
#                 Point(0.0, 0.0),
#                 Point(DOMAIN_LENGTH, 0.0),
#                 Point(DOMAIN_LENGTH, DOMAIN_LENGTH),
#                 Point(0.0, DOMAIN_LENGTH),
#             ]

# domain = Polygon(domain_vertices)

# mesh = generate_mesh(domain, N)


# plot(mesh, interactive=True)

P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)

# CH = FunctionSpace(mesh, MixedElement([P1, P1, P1, P1, P1]))

CH = FunctionSpace(mesh, P1*P1)

dch = TrialFunction(CH)
h_1, j_1 = TestFunctions(CH)

ch = Function(CH)
ch0 = Function(CH)

# Split mixed functions
da, dmu_AB = split(dch)
x_a, mu_AB = split(ch)
x_a0, mu0_AB = split(ch0)

x_a = variable(x_a)

ch_init = InitialConditions(degree=1)
ch.interpolate(ch_init)
ch0.interpolate(ch_init)

kappa = (2.0/3.0)*chi_AB

#Flory-Huggins Expression
g = ( x_a * ln(x_a) / N_A ) + ((1.0-x_a)*ln(1-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB

# Using the fenics autodifferentiation toolkit 
dgdx_a = diff(g,x_a)

# Plot Gibbs energy functional and derivative

# molfracs = np.arange(0.0, 2.0, 0.01)
# fig, ax = plt.subplots()
# line1, = ax.plot (molfracs, g, label ='Gibbs energy of mixing')
# # line2, = ax.plot(x, dgdx_a, label='Derivative')
# ax.set(xlabel='mole fraction of A', ylabel='Energy')
# fig.savefig('gibbs.png')
# plt.show()

# Introduce an expression for mu_{n+theta}
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
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6

file_a = File("concentration_A.pvd", "compressed")

t = 0.0

while (t < TIME_MAX):
    print (t)
    t += dt
    ch0.vector()[:] = ch.vector()
    solver.solve(problem, ch.vector())
    file_a << (ch.split()[0], t)

plot(ch.split()[0])
interactive()