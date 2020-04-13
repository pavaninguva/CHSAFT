from pystencils.session import *

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
    chi_AB,
    N_A,
    N_B,
    GIBBS
)



# Defining the mesh and field
dh = ps.create_data_handling(domain_size=(50, 50), periodicity=True)
μ_field = dh.add_array('mu', latex_name='μ')
c_field = dh.add_array('c')

# Definiting the free energy functional
κ, chi = sp.symbols("κ chi")

c = c_field.center
μ = μ_field.center



def f(c):
    return (c*sp.ln(c)/ N_A) + ((1-c)*sp.ln(1-c)/ N_B) + chi*c*(1-c)
 
bulk_free_energy_density = f(c)
grad_sq = sum(ps.fd.diff(c, i)**2 for i in range(dh.dim))
interfacial_free_energy_density = κ/2 * grad_sq

free_energy_density = bulk_free_energy_density + interfacial_free_energy_density
free_energy_density

# Discretization

discretize = ps.fd.Discretization2ndOrder(dx=10, dt=0.01)

μ_update_eq = ps.fd.functional_derivative(free_energy_density, c)
μ_update_eq = ps.fd.expand_diff_linear(μ_update_eq, constants=[κ])  # pull constant κ in front of the derivatives
μ_update_eq_discretized = discretize(μ_update_eq)
μ_update_eq_discretized


μ_kernel = ps.create_kernel([ps.Assignment(μ_field.center,
                                           μ_update_eq_discretized.subs(chi, chi_AB).subs(κ, chi_AB/3))]
                           ).compile()

M = sp.Symbol("M")
cahn_hilliard = ps.fd.transient(c) - ps.fd.diffusion(μ, M)
cahn_hilliard


c_update = discretize(cahn_hilliard)
c_update

c_kernel = ps.create_kernel([ps.Assignment(c_field.center,
                                           c_update.subs(M, (c*(1-c))))]
                           ).compile()

def init(value=0.4, noise=NOISE_MAGNITUDE):
    for b in dh.iterate():
        b['c'].fill(value)
        np.add(b['c'], noise*np.random.rand(*b['c'].shape), out=b['c'])


def timeloop(steps=100000):
    c_sync = dh.synchronization_function(['c'])
    μ_sync = dh.synchronization_function(['mu'])
    for t in range(steps):
        c_sync()
        dh.run_kernel(μ_kernel)
        μ_sync()
        dh.run_kernel(c_kernel)
    return dh.gather_array('c')
init()


if 'is_test_run' in globals():
    timeloop(10)
    result = None
else:
    ani = ps.plot.scalar_field_animation(timeloop, rescale=True, frames=600, cmap="RdBu")
    plt.show()
    # result = ps.jupyter.display_as_html_video(ani)
result