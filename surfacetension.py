import random
from dolfin import *
import csv
import os
import time
import sys
from FHTaylor import taylorapprox_fullFH, taylorapprox_logonlyFH, heatofmixing, symquarticdoublewell, polynomialfit_fullFH
from Thermo.Thermo import RK, ThermoMix

from parameters.params import (
    A_RAW,
    NOISE_MAGNITUDE,
    TIME_MAX,
    DT,
    N_CELLS,
    DOMAIN_LENGTH,
    theta_ch,
    MESH_TYPE,
    TIME_STRIDE,
    SPECIES,
    chi_AB,
    N_A,
    N_B,
    GIBBS,
    TEMPERATURE,
    FINITE_ELEMENT,
    FINITE_ELEMENT_ORDER,
    SOLVER_CONFIG,
    MOBILITY_MODEL,
    SIZE_DISPARITY,
    A_SYM,
    MOBILITY_MODEL
)

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
parameters["form_compiler"]["quadrature_degree"] = 2

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


N = int(N_CELLS)

mesh = IntervalMesh(N, 0.0,DOMAIN_LENGTH)

P1 = FiniteElement(FINITE_ELEMENT, mesh.ufl_cell(), FINITE_ELEMENT_ORDER)
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

if GIBBS == "FH":
    # r = RK("FH",["PB","PS"],[N_A,N_B])
    g = ( x_a * ln(x_a) )/N_A + ((1.0-x_a)*ln(1-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB 
    # g = r.G_RK(x_a)
    print("Vanilla FH")
if GIBBS == "UNIFAC":
    r = RK("UNIFAC",SPECIES,[N_A,N_B],[TEMPERATURE])
    g = r.G_RK(x_a)
    print("Redlich-Kister UNIFAC")
elif GIBBS == "PCSAFT_CR":
    r = RK("PCSAFT",SPECIES,[N_A,N_B],[TEMPERATURE],CR="On")
    g = r.G_RK(x_a)
    print("Redlich-Kister PCSAFT")
elif GIBBS == "PCSAFT_Fit":
    r = RK("PCSAFT",SPECIES,[N_A,N_B],[TEMPERATURE],CR="Off")
    g = r.G_RK(x_a)
    print("Redlich-Kister PCSAFT")
elif GIBBS == "TaylorApproxFullFH":
    g = taylorapprox_fullFH(N_A, N_B, chi_AB, x_a)
    print("full taylor approx of FH")
elif GIBBS == "TaylorApproxLogOnlyFH":
    g = taylorapprox_logonlyFH(N_A, N_B, chi_AB, x_a)
    print ("Taylor approx of log term in FH only")
elif GIBBS == "polynomialfit_fullFH":
    g = polynomialfit_fullFH(N_A, N_B, chi_AB, x_a, 4)
    print("Full curve fit ahaha")
elif GIBBS == "symquarticdoublewell_fullFH":
    g = symquarticdoublewell(x_a, A_SYM)
    print ("BS quartic polynomial")
elif GIBBS == "Heatofmixing":
    g = heatofmixing(chi_AB, x_a)
else: 
    print ("work harder")

if SIZE_DISPARITY == "SMALL":
    if GIBBS=="FH":
        kappa = (2.0/3.0)*chi_AB
    else:
        chi_AB = r.RK(0)/(N_A*N_B)**0.5
        kappa = (2.0/3.0)*chi_AB[0]
    print ("about the same size")
elif SIZE_DISPARITY == "LARGE": 
    if GIBBS=="FH":
        kappa = (1.0/3.0)*chi_AB
    else:
        chi_AB = r.RK(0)/(N_A*N_B)**0.5
        kappa = (1.0/3.0)*chi_AB[0]
    print ("big size difference")

print(chi_AB)

# Using the fenics autodifferentiation toolkit 
dgdx_a = diff(g,x_a)

mu_AB_mid = (1.0 - theta_ch) * mu0_AB + theta_ch * mu_AB

dt = DT

F_a = (
    x_a * h_1 * dx
    - x_a0 * h_1 * dx
    + dt * x_a * (1.0 - x_a) * dot(grad(mu_AB_mid), grad(h_1)) * dx
)

F_a_constant = (
    x_a * h_1 * dx
    - x_a0 * h_1 * dx
    + dt * dot(grad(mu_AB_mid), grad(h_1)) * dx
)

F_mu_AB = (
    mu_AB * j_1 * dx
    - dgdx_a * j_1 * dx
    - kappa * dot(grad(x_a), grad(j_1)) * dx
)

if MOBILITY_MODEL == "Variable":
    F = F_a + F_mu_AB
    print ("work hard model")
elif MOBILITY_MODEL == "Constant":
    F = F_a_constant + F_mu_AB
    print("less work hehe")
else:
    print("wrong model implemented")
    sys.exit()

a = derivative(F, ch, dch)

if SOLVER_CONFIG == "LU":

    problem = CahnHilliardEquation(a, F)
    # problem.set_bounds(x_a_min, x_a_max)
    solver = NewtonSolver()
    solver.parameters["linear_solver"] = "lu"
    # solver.parameters["linear_solver"] = "gmres"
    # solver.parameters["preconditioner"] = "ilu"
    solver.parameters["convergence_criterion"] = "residual"
    solver.parameters["relative_tolerance"] = 1e-10
    solver.parameters["absolute_tolerance"] = 1e-16
    solver.parameters["relaxation_parameter"] = 0.5

elif SOLVER_CONFIG == "KRYLOV":
    class CustomSolver(NewtonSolver):

        def __init__(self):
            NewtonSolver.__init__(self, mesh.mpi_comm(),
                                PETScKrylovSolver(), PETScFactory.instance())

        def solver_setup(self, A, P, problem, iteration):
            self.linear_solver().set_operator(A)

            PETScOptions.set("ksp_type", "gmres")
            PETScOptions.set("ksp_monitor")
            PETScOptions.set("pc_type", "hypre")
            PETScOptions.set("pc_hypre_type", "euclid")
            PETScOptions.set("ksp_rtol", "1.0e-8")
            PETScOptions.set("ksp_atol", "1.0e-16")
            PETScOptions.set('ksp_max_it', '1000')

            self.linear_solver().set_from_options()

    problem = CahnHilliardEquation(a, F)
    solver = CustomSolver()


# Initialising the output files
gibbs_list = []
surface_tension = []
# file = XDMFFile("output.xdmf")
file = File("output.pvd", "compressed")

start = time.time()

t = 0.0
time_stride = TIME_STRIDE
timestep = 0
while (t < TIME_MAX):
    t += dt
    print (f"Time = {t}")
    ch0.vector()[:] = ch.vector()
    solver.solve(problem, ch.vector())
    timestep += 1

    for i in ch.vector()[::2]:
        if i < 0.0:
            i == 0.0 + 1e-16
        elif i > 1.0:
            i == 1.0 - 1e-16

    # Assembling the various terms of the Landau-Ginzburg free energy functional
    homogenous_energy = assemble(g * dx())
    gradient_energy = assemble(Constant(0.5)*kappa * dot(grad(x_a), grad(x_a)) * dx())
    gibbs= homogenous_energy + gradient_energy
    gibbs_list.append(gibbs)

    # Assembling the suface tension contribution
    Gamma = assemble(kappa*dot(grad(x_a), grad(x_a))*dx())
    surface_tension.append(Gamma)

    # fpath1 = "./output_gibbs.csv"
    fpath2 = "./output_surfacetension.csv"
    # headers1 = ["time", "gibbs"]
    headers2 = ["time", "Gamma"]
    end = time.time()
    print(end-start)
    # Write header row (for first timestep)
    if not os.path.exists(fpath2):
        with open(fpath2,"w") as f:
            w = csv.DictWriter(f,headers2)
            w.writeheader()

    if (timestep % time_stride ==0):
        # file.write (ch.split()[0], t)
        file << (ch.split()[0], t)
        # Appending case data
        with open(fpath2, "a") as f:
            w = csv.DictWriter(f, headers2)
            w.writerow({"time": float(t), "Gamma": float(Gamma)})