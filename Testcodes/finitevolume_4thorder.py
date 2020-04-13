from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver, parallelComm, ExplicitDiffusionTerm, ConvectionTerm
from fipy.tools import numerix
import time
from fipy import NthOrderBoundaryCondition, PeriodicGrid2D



A_RAW = 0.5
NOISE_MAGNITUDE = 0.01
TIME_MAX = 100
DT = 0.5
TIME_STRIDE = 10
chi_AB = 0.08
N_A = 1000
N_B = 1000


print ("Yay")

# # Define mesh
mesh = PeriodicGrid2D(nx=10.0, ny=10.0, dx=1.0, dy=1.0)
print ("mesh loaded")

x_a = CellVariable(name=r"x_a", mesh = mesh, hasOld=1)
xi = CellVariable(name=r"xi", mesh = mesh, hasOld=1)

x_a.constrain(x_a.faceValue, mesh.facesLeft)
x_a.constrain(x_a.faceValue, mesh.facesRight)

xi.constrain(xi.faceValue, mesh.facesLeft)
xi.constrain(xi.faceValue, mesh.facesRight)

# We need to introduce the noise
noise = GaussianNoiseVariable(mesh=mesh,
                              mean = A_RAW,
                              variance = NOISE_MAGNITUDE).value

x_a[:] = noise


dgdx_a = ((1.0/N_A) - (1.0/N_B)) + (1.0/N_A)*numerix.log(x_a) - (1.0/N_B)*numerix.log(1.0 - x_a) + chi_AB*(1.0 - 2*x_a)


kappa = (1.0/6.0)*chi_AB

# eq1 = TransientTerm(var=x_a) == ConvectionTerm(coeff= ((x_a*(1.0 - x_a))*[[1]]).rank *xi.faceGrad) + DiffusionTerm((x_a *(1.0 - x_a), kappa), var = x_a)
eq1 = TransientTerm(var=x_a) == DiffusionTerm(coeff= x_a*(1.0 - x_a), var = xi) + DiffusionTerm((x_a *(1.0 - x_a), kappa), var = x_a)

eq2 = ImplicitSourceTerm(coeff=1, var= xi) == dgdx_a

eq = eq1 & eq2


elapsed = 0.
dt = DT
if __name__ == "__main__":
    duration = TIME_MAX


time_stride = TIME_STRIDE
timestep = 0

# Defining the solver
solver = LinearLUSolver(tolerance=1e-10, iterations=50)

start = time.time()

while elapsed < duration: 
    if (timestep == 0):
        vw = VTKCellViewer(vars=(x_a, xi))
        vw.plot(filename="0_output.%d.vtk" %(parallelComm.procID))
    elapsed += dt
    timestep += 1
    x_a.updateOld()
    xi.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        print ("sweep!")
    print (elapsed)
    end = time.time()
    print(end-start)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(x_a, xi))
        vw.plot(filename="%s_output.%d.vtk" %(elapsed, parallelComm.procID))