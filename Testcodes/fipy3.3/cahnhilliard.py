from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver,  PeriodicGrid2D
from fipy.tools import numerix
from math import log


mesh =  PeriodicGrid2D(nx=128, ny=128, dx=1.0, dy=1.0)
phi = CellVariable(name=r"$\phi$", mesh=mesh)
psi = CellVariable(name=r"$\psi$", mesh=mesh)

noise = GaussianNoiseVariable(mesh=mesh, mean=0.5, variance=0.01).value
phi[:] = noise



D = a = epsilon = 1.
# kappa = log(phi)
N_A = 400
N_B = 400
chi_AB = 0.0075
kappa = chi_AB / 6.0
# dfdphi = a**2 * phi * (1 - phi) * (1 - 2 * phi)
# dfdphi_ = a**2 * (1 - phi) * (1 - 2 * phi)
# d2fdphi2 = a**2 * (1 - 6 * phi * (1 - phi))

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

while elapsed < duration: 
    elapsed += dt
    eq.solve(dt=dt, solver=solver)
    phi.setValue(0.0, where=phi < 0.0 )
    phi.setValue(1.0, where=phi > 1.0 )
    print (elapsed)
    if (elapsed % 50 ==0):
        vw = VTKCellViewer(vars=(phi, psi))
        vw.plot(filename="%s_vashista_output.vtk" %(elapsed)) 

    