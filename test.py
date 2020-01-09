from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, Viewer, VTKCellViewer
from fipy.tools import numerix
from autograd import grad
from autograd.numpy import log

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
# def g(phi):
#     return phi**2 *(1.0-phi)**2
# dfdphi = grad(g)
# d2fdphi2 = grad(dfdphi)

dfdphi = a**2 * phi * (1 - phi) * (1 - 2 * phi)
print (type (dfdphi))
dfdphi_ = a**2 * (1 - phi) * (1 - 2 * phi)
d2fdphi2 = a**2 * (1 - 6 * phi * (1 - phi))
eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
eq2 = (ImplicitSourceTerm(coeff=1., var=psi)
        == ImplicitSourceTerm(coeff=d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi
        - DiffusionTerm(coeff=epsilon**2, var=phi))

eq = eq1 & eq2      

elapsed = 0.
dt = 0.5
if __name__ == "__main__":
    duration = 100.



time_stride = 100
timestep = 0

while elapsed < duration: 
    elapsed += dt
    timestep += 1
    eq.solve(dt=dt)
    print (elapsed)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(phi, psi))
        vw.plot(filename="%s_output.vtk" %elapsed)
    if __name__ == "__main__":
        viewer.plot()

    
from builtins import input

if __name__ == '__main__':
    input("Coupled equations. Press <return> to proceed...")