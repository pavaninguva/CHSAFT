from fipy import CellVariable, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver,  PeriodicGrid1D, Viewer, UniformNoiseVariable, Grid1D
from fipy.tools import numerix
from math import log

nx = 100
dx = 0.3

mesh = Grid1D(dx=dx, nx=nx)
phi = CellVariable(name=r"$\phi$", mesh=mesh)
psi = CellVariable(name=r"$\psi$", mesh=mesh)

# noise = GaussianNoiseVariable(mesh=mesh, mean=0.8, variance=0.002).value
# phi[:] = noise

x = mesh.cellCenters[0]
phi.value = 0.1
phi.setValue(0.2, where=x < 0.3*(nx*dx))
phi.setValue(0.6, where=x > 0.3*(nx*dx))
phi.setValue(0.2, where=x > 0.32*(nx*dx))
phi.setValue(0.6, where=x > 0.7*(nx*dx))
phi.setValue(0.2, where=x > 0.72*(nx*dx))



D = a = epsilon = 1.
N_A = 500
N_B = 500
chi_AB = 0.015
kappa = 2*chi_AB / 3.0

dfdphi = ((1.0/N_A) - (1.0/N_B)) + (1.0/N_A)*numerix.log(phi) - (1.0/N_B)*numerix.log(1.0 - phi) + chi_AB*(1.0 - 2*phi)
d2fdphi2 = (1.0/(N_A*phi)) + (1.0/(N_B*(1.0 - phi))) - 2*chi_AB
eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
eq2 = (ImplicitSourceTerm(coeff=1., var=psi)
        == ImplicitSourceTerm(coeff=d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi
        - DiffusionTerm(coeff=kappa, var=phi))

eq = eq1 & eq2      

elapsed = 0.
dt = 1.0
if __name__ == "__main__":
    duration = 1000.

solver = LinearLUSolver(tolerance=1e-9, iterations=500)

# Create viewer
if __name__ == "__main__":
    viewer = Viewer(name=r"$\phi$", vars=phi, datamin=0.,datamax=1.)

while elapsed < duration: 
    elapsed += dt
    eq.solve(dt=dt, solver=solver)
    # phi.setValue(0.0, where=phi < 0.0 )
    # phi.setValue(1.0, where=phi > 1.0 )
    print (elapsed)
    if (elapsed % 1 ==0):
        if __name__ == "__main__":
            viewer.plot()