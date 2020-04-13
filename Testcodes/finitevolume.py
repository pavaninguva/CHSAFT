from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver, parallelComm, ExplicitDiffusionTerm
from fipy.tools import numerix
import time
from fipy import PeriodicGrid2D


from params_fipy import (
    A_RAW,
    NOISE_MAGNITUDE,
    TIME_MAX,
    DT,
    N_CELLS,
    DOMAIN_LENGTH,
    TIME_STRIDE,
    chi_AB,
    N_A,
    N_B,
    GIBBS,
    DESIRED_RESIDUAL
)

print ("Yay")

# # Define mesh
# mesh = Grid2D(dx=dx, dy=dx, nx=N_CELLS, ny=N_CELLS)
mesh = PeriodicGrid2D(nx=50.0, ny=50.0, dx=0.1, dy=0.1)
print ("mesh loaded")

# We need to define the relevant variables: 
x_a = CellVariable(name=r"x_a", mesh = mesh, hasOld=1)
mu_AB = CellVariable(name=r"mu_AB", mesh = mesh, hasOld=1)


# We need to introduce the noise
noise = GaussianNoiseVariable(mesh=mesh,
                              mean = A_RAW,
                              variance = NOISE_MAGNITUDE).value

x_a[:] = noise

# x_a.setValue(GaussianNoiseVariable(mesh=mesh,
#                                    mean=A_RAW,
#                                    variance=NOISE_MAGNITUDE)
# )

# def g(x_a):
#     if GIBBS == "FH":
#         #Flory-Huggins Expression
#         # The numerix.log is needed since x_a is a CellVariable and not a simple array of size 1. 
#         # The standard math.log cannot cope with this. 
#         return ( x_a * log(x_a) / N_A ) + ((1.0-x_a)*log(1-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB 
#     elif GIBBS != "FH":
#         print ("Other free energy functions yet to be implemented")

# # Use automatic differentiation to evaluate the derivatives of the free energy
# print (g(0.5))
# dgdx_a = grad(g)
# d2gdx_a2 = grad(dgdx_a)
# print (dgdx_a(0.4))
if GIBBS == "FH":
    dgdx_a = ((1.0/N_A) - (1.0/N_B)) + (1.0/N_A)*numerix.log(x_a) - (1.0/N_B)*numerix.log(1.0 - x_a) + chi_AB*(1.0 - 2*x_a)
    d2gdx_a2 = (1.0/(N_A*x_a)) + (1.0/(N_B*(1.0 - x_a))) - 2*chi_AB
elif GIBBS != "FH": 
    print("Implent more stuff you lazy fuck")

# Define the equations

# evaluating kappa
kappa = (2.0/3.0)*chi_AB

# eqn 1 is the 2nd order transport equation
eq1 = (TransientTerm(var=x_a)) == DiffusionTerm(coeff = x_a * (1 - x_a), var=mu_AB)

# Try constant mobility
eq0 = (TransientTerm(var=x_a)) == DiffusionTerm(coeff = 1, var=mu_AB)

# eqn 2 is the chemical potential definition
eq2 = (ImplicitSourceTerm(coeff=1. , var=mu_AB)) == ImplicitSourceTerm(coeff=d2gdx_a2, var=x_a) - d2gdx_a2 * x_a + dgdx_a - DiffusionTerm(coeff=kappa, var=x_a)

# write eq2 without the fipy trick:
eq3 = (ImplicitSourceTerm(coeff=1, var=mu_AB)) == dgdx_a - DiffusionTerm(coeff= kappa, var = x_a)

# Adding the equations together
eq = eq1 & eq2

elapsed = 0.
dt = DT
if __name__ == "__main__":
    duration = TIME_MAX


time_stride = TIME_STRIDE
timestep = 0

# Defining the solver to improve numerical stabilty
solver = LinearLUSolver(tolerance=1e-9, iterations=50, precon= "lu")
# solver = PETSc.KSP().create()
start = time.time()

while elapsed < duration: 
    if (timestep == 0):
        vw = VTKCellViewer(vars=(x_a, mu_AB))
        vw.plot(filename="0_output.%d.vtk" %(parallelComm.procID))
    elapsed += dt
    timestep += 1
    x_a.updateOld()
    mu_AB.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        x_a.setValue(0.0, where=x_a < 0.0)
        x_a.setValue(1.0, where=x_a > 1.0 )
        print ("sweep!")
    print (elapsed)
    end = time.time()
    print(end-start)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(x_a, mu_AB))
        vw.plot(filename="%s_output.%d.vtk" %(elapsed, parallelComm.procID))









