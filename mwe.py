import random
from dolfin import *
import numpy as np

# from common.constants import (
chi_AB = 0.005
chi_AC = 0.005
chi_BC = 0.005
D_AB = 1e-11
D_AC = 1e-11
D_BC = 1e-11
N_A = 1000
N_B = 1000
N_C = 1000
# kappa_AA 
# kappa_BB,
# kappa_AB,
theta_ch = 0.5
maximum_iterations = 30
# relaxation_parameter,
# convergence_criterion,
# absolute_tolerance,
# relative_tolerance,
N_SCALE_OPTION = N_A
D_SCALE_OPTION = D_AB
# use_ramp,
# step_duration,
# n_steps,
# mult_start,
# mult_end,
A_RAW = 0.3
B_RAW = 0.3
NOISE_MAGNITUDE = 0.03
DT = 0.1
N_CELLS = 80
DOMAIN_LENGTH = 40
TIME_MAX = 5
# )


# Class representing the intial conditions
class InitialConditions(Expression):
    def __init__(self, **kwargs):
        random.seed(1234)

    def eval(self, values, x):
        values[0] = A_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[1] = B_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[2] = 0.0
        values[3] = 0.0
        values[4] = 0.0

    def value_shape(self):
        return (5,)


# Class for interfacing with the Newton solver
class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a

    def F(self, b, x):
        assemble(self.L, tensor=b)

    def J(self, A, x):
        assemble(self.a, tensor=A)


# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"

# Create mesh and build function space

mesh = RectangleMesh(
    Point(0.0, 0.0), Point(DOMAIN_LENGTH, DOMAIN_LENGTH), N_CELLS, N_CELLS
)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)

CH = FunctionSpace(mesh, MixedElement([P1, P1, P1, P1, P1]))

# Define trial and test functions
dch = TrialFunction(CH)
h_1, h_2, j_1, j_2, j_3 = TestFunctions(CH)

ch = Function(CH)
ch0 = Function(CH)

a, b, N_mu_AB, N_mu_AC, N_mu_BC = split(ch)
a0, b0, N_mu0_AB, N_mu0_AC, N_mu0_BC = split(ch0)

a = variable(a)
b = variable(b)

# Create intial conditions and interpolate
ch_init = InitialConditions(degree=1)
ch.interpolate(ch_init)
ch0.interpolate(ch_init)

N_scale_options = {"N_A": N_A, "N_B": N_B, "N_C": N_C}
N_SCALE = N_scale_options[N_SCALE_OPTION]

D_scale_options = {"D_AB": D_AB, "D_AC": D_AC, "D_BC": D_BC}
D_SCALE = D_scale_options[D_SCALE_OPTION]

kappa_AA = (1.0 / 3.0) * chi_AC
kappa_BB = (1.0 / 3.0) * chi_BC
kappa_AB = (1.0 / 6.0) * (chi_AC +  chi_BC - 2.0 * chi_AB)

N_kappa_AA = N_SCALE * kappa_AA
N_kappa_BB = N_SCALE * kappa_BB
N_kappa_AB = N_SCALE * kappa_AB

# "ordinary free energy" gradients
N_dgda = N_SCALE * (
    (1.0 / N_A) * (1.0 + ln(a)) + chi_AB * b + chi_AC * (1.0 - a - b)
)
N_dgdb = N_SCALE * (
    (1.0 / N_B) * (1.0 + ln(b)) + chi_AB * a + chi_BC * (1.0 - a - b)
)
N_dgdc = N_SCALE * (
    (1.0 / N_C) * (1.0 + ln(1.0 - a - b)) + chi_AC * a +  chi_BC * b
)

N_mu_AB_mid = (1.0 - theta_ch) * N_mu0_AB + theta_ch * N_mu_AB
N_mu_AC_mid = (1.0 - theta_ch) * N_mu0_AC + theta_ch * N_mu_AC
N_mu_BC_mid = (1.0 - theta_ch) * N_mu0_BC + theta_ch * N_mu_BC

# scale diffusivity
D_AB_ = D_AB / D_SCALE
D_AC_ = D_AC / D_SCALE
D_BC_ = D_BC / D_SCALE


dt = DT

# transport equations
F_a = (
    a * h_1 * dx
    - a0 * h_1 * dx
    + dt * a * b * D_AB_ * dot(grad(N_mu_AB_mid), grad(h_1)) * dx
    + dt * a * (1.0 - a - b) * D_AC_ * dot(grad(N_mu_AC_mid), grad(h_1)) * dx
)

F_b = (
    b * h_2 * dx
    - b0 * h_2 * dx
    - dt * a * b * D_AB_ * dot(grad(N_mu_AB_mid), grad(h_2)) * dx
    + dt * b * (1.0 - a - b) * D_BC_ * dot(grad(N_mu_BC_mid), grad(h_2)) * dx
)

# chemical potential equations
F_N_mu_AB = (
    N_mu_AB * j_1 * dx
    - N_dgda * j_1 * dx
    + N_dgdb * j_1 * dx
    - (N_kappa_AA - N_kappa_AB) * dot(grad(a), grad(j_1)) * dx
    + (N_kappa_BB - N_kappa_AB) * dot(grad(b), grad(j_1)) * dx
)

F_N_mu_AC = (
    N_mu_AC * j_2 * dx
    - N_dgda * j_2 * dx
    + N_dgdc * j_2 * dx
    - N_kappa_AA * dot(grad(a), grad(j_2)) * dx
    - N_kappa_AB * dot(grad(b), grad(j_2)) * dx
)

F_N_mu_BC = (
    N_mu_BC * j_3 * dx
    - N_dgdb * j_3 * dx
    + N_dgdc * j_3 * dx
    - N_kappa_BB * dot(grad(b), grad(j_3)) * dx
    - N_kappa_AB * dot(grad(a), grad(j_3)) * dx
)

F = F_a + F_b + F_N_mu_AB + F_N_mu_AC + F_N_mu_BC

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(F, ch, dch)

# # Create nonlinear problem and Newton solver
problem = CahnHilliardEquation(a, F)

solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6


# def Max(a, b): return (a+b+abs(a-b))/Constant(2)
# def Min(a, b): return (a+b-abs(a-b))/Constant(2)

# Output file
file_a = File("concentration_A.pvd", "compressed")
file_b = File("concentration_B.pvd", "compressed")

# Step in time
t = 0.0
timestep = 0

# output initial condition
file_a << (ch.split()[0], t)
file_b << (ch.split()[1], t)

space = FunctionSpace(mesh, P1)

while t < TIME_MAX:


    timestep += 1
    t += dt

    if MPI.rank(mpi_comm_world()) == 0:
        print ("Timestep", timestep, "Time", t)


    ch0.vector()[:] = ch.vector()
    solver.solve(problem, ch.vector())
    
    # if timestep % 10 == 0:
    file_a << (ch.split()[0], t)
    file_b << (ch.split()[1], t)

    proj_a = project(ch.split()[0], FunctionSpace(mesh, P1))
    proj_b = project(ch.split()[1], FunctionSpace(mesh, P1))

    gather_a = Vector()
    gather_b = Vector()
    
    proj_a.vector().gather(gather_a, np.array(range(space.dim()), "intc"))
    proj_b.vector().gather(gather_b, np.array(range(space.dim()), "intc"))

    if MPI.rank(mpi_comm_world()) == 0:
    

        print gather_a.array().shape
        print gather_a.array().min()
        print gather_a.array().max()

        print gather_b.array().shape
        print gather_b.array().min()
        print gather_b.array().max()
