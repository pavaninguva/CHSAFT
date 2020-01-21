from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver
from fipy.tools import numerix

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
# dx = DOMAIN_LENGTH / N_CELLS
# dy = DOMAIN_LENGTH / N_CELLS

# # Define mesh
# mesh = Grid2D(dx=dx, dy=dx, nx=N_CELLS, ny=N_CELLS)
mesh = Grid2D(nx=100.0, ny=100.0, dx=0.5, dy=0.5)

# We need to define the relevant variables: 
x_a = CellVariable(name=r"x_a", mesh = mesh, hasOld=1)
mu_AB = CellVariable(name=r"mu_AB", mesh = mesh, hasOld=1)

# We need to introduce the noise
noise = GaussianNoiseVariable(mesh=mesh,
                              mean = A_RAW,
                              variance = NOISE_MAGNITUDE).value

x_a[:] = noise

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
kappa = (1.0/6.0)*chi_AB

# eqn 1 is the 2nd order transport equation
eq1 = (TransientTerm(var=x_a)) == DiffusionTerm(coeff = x_a * (1 - x_a), var=mu_AB)

# eqn 2 is the chemical potential definition
eq2 = (ImplicitSourceTerm(coeff=1. , var=mu_AB)) == ImplicitSourceTerm(coeff=d2gdx_a2, var=x_a) - d2gdx_a2 * x_a + dgdx_a - DiffusionTerm(coeff=kappa, var=x_a)

# Adding the equations together
eq = eq1 & eq2

elapsed = 0.
dt = DT
if __name__ == "__main__":
    duration = TIME_MAX


time_stride = TIME_STRIDE
timestep = 0

# Defining the solver to improve numerical stabilty
solver = LinearLUSolver(tolerance=1e-10, iterations=25, maxIterations=50)

while elapsed < duration: 
    elapsed += dt
    timestep += 1
    x_a.updateOld()
    mu_AB.updateOld()
    res = 1e+10
    # while res > DESIRED_RESIDUAL:
    # eq.solve(dt=dt)
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        print ("sweep!")
    print (elapsed)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(x_a, mu_AB))
        vw.plot(filename="%s_output.vtk" %elapsed)









