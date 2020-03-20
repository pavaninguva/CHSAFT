from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix
from math import log

if __name__ == "__main__":
    nx = ny = 20
else:
        nx = ny = 10
mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
phi = CellVariable(name=r"$\phi$", mesh=mesh)
psi = CellVariable(name=r"$\psi$", mesh=mesh)

noise = GaussianNoiseVariable(mesh=mesh, mean=0.5, variance=0.01).value
phi[:] = noise

if __name__ == "__main__":
    viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

D = a = epsilon = 1.
# kappa = log(phi)
N_A = 1000
N_B = 1000
chi_AB = 0.03
# dfdphi = a**2 * phi * (1 - phi) * (1 - 2 * phi)
# dfdphi_ = a**2 * (1 - phi) * (1 - 2 * phi)
# d2fdphi2 = a**2 * (1 - 6 * phi * (1 - phi))

dfdphi = ((1.0/N_A) - (1.0/N_B)) + (1.0/N_A)*numerix.log(phi) - (1.0/N_B)*numerix.log(1.0 - phi) + chi_AB*(1.0 - 2*phi)
d2fdphi2 = (1.0/(N_A*phi)) + (1.0/(N_B*(1.0 - phi))) - 2*chi_AB
eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
eq2 = (ImplicitSourceTerm(coeff=1., var=psi)
        == ImplicitSourceTerm(coeff=d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi
        - DiffusionTerm(coeff=epsilon**2, var=phi))


eq = eq1 & eq2      

elapsed = 0.
dt = 0.5
if __name__ == "__main__":
    duration = 100.


while elapsed < duration: 
    elapsed += dt
    eq.solve(dt=dt)

    if __name__ == "__main__":
        viewer.plot()

from builtins import input

if __name__ == '__main__':
    input("Coupled equations. Press <return> to proceed...")