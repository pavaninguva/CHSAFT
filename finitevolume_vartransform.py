from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver, parallelComm
from fipy.tools import numerix
import time


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
mesh = Grid2D(nx=10.0, ny=10.0, dx=1.0, dy=1.0)
print ("mesh loaded")

# We need to define the relevant variables: 
a = CellVariable(name=r"a", mesh = mesh, hasOld=1)
exp_a = CellVariable(name=r"exp_a", mesh = mesh, hasOld=1)
mu_AB = CellVariable(name=r"mu_AB", mesh = mesh, hasOld=1)


# We need to introduce the noise
noise = GaussianNoiseVariable(mesh=mesh,
                              mean = numerix.log(A_RAW),
                              variance = abs(numerix.log(NOISE_MAGNITUDE))).value

a[:] = noise

if GIBBS == "FH":
    dgda = (1+a)*(numerix.exp(a)/N_A) - (numerix.exp(a)/ N_B)*(1+ numerix.log(1.0 - numerix.exp(a))) + chi_AB*numerix.exp(a) -2.0*chi_AB*numerix.exp(2.0*a)
elif GIBBS != "FH": 
    print("Implent more stuff you lazy fuck")

# Define the equations

# evaluating kappa
kappa = (1.0/6.0)*chi_AB

# # eqn 1 is the 2nd order transport equation
# eq1 = (TransientTerm(var=x_a)) == DiffusionTerm(coeff = x_a * (1 - x_a), var=mu_AB)

# # eqn 2 is the chemical potential definition
# eq2 = (ImplicitSourceTerm(coeff=1. , var=mu_AB)) == ImplicitSourceTerm(coeff=d2gdx_a2, var=x_a) - d2gdx_a2 * x_a + dgdx_a - DiffusionTerm(coeff=kappa, var=x_a)

# eq1 is the transport equation
eq1 = (TransientTerm(coeff=numerix.exp(a), var=a)) == DiffusionTerm(coeff = numerix.exp(a)*(1.0- numerix.exp(a)), var=mu_AB)

#eq2 is the chemical potential
eq2 = (ImplicitSourceTerm(coeff=1. , var=mu_AB)) == dgda*(1.0/numerix.exp(a)) - DiffusionTerm(coeff=kappa, var=exp_a)

eq3 = (ImplicitSourceTerm(coeff=1, var = exp_a)) == numerix.exp(a)
# Adding the equations together
eq = eq1 & eq2 & eq3

elapsed = 0.
dt = DT
if __name__ == "__main__":
    duration = TIME_MAX


time_stride = TIME_STRIDE
timestep = 0

# Defining the solver to improve numerical stabilty
solver = LinearLUSolver(tolerance=1e-10, iterations=25)
# solver = PETSc.KSP().create()
start = time.time()

while elapsed < duration: 
    if (timestep == 0):
        vw = VTKCellViewer(vars=(a, mu_AB))
        vw.plot(filename="0_output.%d.vtk" %(parallelComm.procID))
    elapsed += dt
    timestep += 1
    a.updateOld()
    mu_AB.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        print ("sweep!")
    print (elapsed)
    end = time.time()
    print(end-start)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(a, mu_AB))
        vw.plot(filename="%s_output.%d.vtk" %(elapsed, parallelComm.procID))