from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix
from autograd import grad
from math import log

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
    GIBBS
)

# Define mesh
mesh = Grid2D(Lx=DOMAIN_LENGTH, Ly=DOMAIN_LENGTH, nx=N_CELLS, ny=N_CELLS)

# We need to define the relevant variables: 
x_a = CellVariable(name=r"$\x_a$", mesh = mesh)
mu_AB = CellVariable(name=r"$\mu_AB$", mesh = mesh)

# We need to introduce the noise
noise = GaussianNoiseVariable(mesh=mesh,
                              mean = A_RAW,
                              variance = NOISE_MAGNITUDE).value

x_a[:] = noise

if __name__ == "__main__":
    viewer = Viewer(vars=(x_a,), datamin=0., datamax=1.)

print (type(x_a))
if GIBBS == "FH":
    #Flory-Huggins Expression
    g = ( x_a * log(x_a) / N_A ) + ((1.0-x_a)*log(1-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB 
elif GIBBS != "FH":
    print ("Other free energy functions yet to be implemented")

# Use automatic differentiation to evaluate the derivatives of the free energy
dgdx_a = grad(g, x_a)
d2gdx_a2 = grad(dgdx_a, x_a)

# Define the equations

# evaluating kappa
kappa = (2.0/3.0)*chi_AB

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


while elapsed < duration: 
    elapsed += dt
    eq.solve(dt=dt)

    if __name__ == "__main__":
        viewer.plot()

from builtins import input

if __name__ == '__main__':
    input("Coupled equations. Press <return> to proceed...")




