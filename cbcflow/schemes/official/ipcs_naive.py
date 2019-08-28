# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
r"""
This incremental pressure correction scheme (IPCS) is an operator splitting scheme that
follows the idea of Goda [1]_.
This scheme preserves the exact same stability properties
as Navier-Stokes and hence does not introduce additional dissipation in the flow.

The idea is to replace the unknown pressure with an approximation. This is chosen as
the pressure solution from the previous solution.

The time discretization is done using backward Euler, the diffusion term is handled with Crank-Nicholson, and the convection is handled explicitly, making the
equations completely linear. Thus, we have a discretized version of the Navier-Stokes equations as

.. math:: \frac{1}{\Delta t}\left( u^{n+1}-u^{n} \right)-\nabla\cdot\nu\nabla u^{n+\frac{1}{2}}+u^n\cdot\nabla u^{n}+\frac{1}{\rho}\nabla p^{n+1}=f^{n+1}, \\
    \nabla \cdot u^{n+1} = 0,

where :math:`u^{n+\frac{1}{2}} = \frac{1}{2}u^{n+1}+\frac{1}{2}u^n.`

For the operator splitting, we use the pressure solution from the previous timestep as an estimation, giving an equation for a tentative velocity, :math:`\tilde{u}^{n+1}`:

.. math:: \frac{1}{\Delta t}\left( \tilde{u}^{n+1}-u^{n} \right)-\nabla\cdot\nu\nabla \tilde{u}^{n+\frac{1}{2}}+u^n\cdot\nabla u^{n}+\frac{1}{\rho}\nabla p^{n}=f^{n+1}.

This tenative velocity is not divergence free, and thus we define a velocity correction :math:`u^c=u^{n+1}-\tilde{u}^{n+1}`. Substracting the second equation from the first, we see that

.. math::
    \frac{1}{\Delta t}u^c-\frac{1}{2}\nabla\cdot\nu\nabla u^c+\frac{1}{\rho}\nabla\left( p^{n+1} - p^n\right)=0, \\
    \nabla \cdot u^c = -\nabla \cdot \tilde{u}^{n+1}.

The operator splitting is a first order approximation, :math:`O(\Delta t)`, so we can, without reducing the order of the approximation simplify the above to

.. math::
    \frac{1}{\Delta t}u^c+\frac{1}{\rho}\nabla\left( p^{n+1} - p^n\right)=0, \\
    \nabla \cdot u^c = -\nabla \cdot \tilde{u}^{n+1},

which is reducible to a Poisson problem:

.. math::
   \Delta p^{n+1} = \Delta p^n+\frac{\rho}{\Delta t}\nabla \cdot \tilde{u}^{n+1}.

The corrected velocity is then easily calculated from

.. math::
    u^{n+1} = \tilde{u}^{n+1}-\frac{\Delta t}{\rho}\nabla\left(p^{n+1}-p^n\right)

The scheme can be summarized in the following steps:
    #. Replace the pressure with a known approximation and solve for a tenative velocity :math:`\tilde{u}^{n+1}`.

    #. Solve a Poisson equation for the pressure, :math:`p^{n+1}`

    #. Use the corrected pressure to find the velocity correction and calculate :math:`u^{n+1}`

    #. Update t, and repeat.

.. [1] Goda, Katuhiko. *A multistep technique with implicit difference schemes for calculating two-or three-dimensional cavity flows.* Journal of Computational Physics 30.1 (1979): 76-95.

"""

from __future__ import division

from cbcflow.core.nsscheme import *

from cbcflow.schemes.utils import (compute_regular_timesteps,
                                   assign_ics_split,
                                   make_velocity_bcs,
                                   make_pressure_bcs,
                                   NSSpacePoolSplit)

def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, mu):
    "Return stress tensor."
    return 2.0*mu*epsilon(u) - p*Identity(len(u))


class IPCS_Naive(NSScheme):
    "Incremental pressure-correction scheme, naive implementation."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            #theta = 0.5,
            )
        return params

    def solve(self, problem, timer):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = NSSpacePoolSplit(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = Function(V, name="u0")
        u1 = Function(V, name="u1")
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        u1.assign(u0)
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        rho = float(problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step
        u_mean = 0.5 * (u + u0)
        u_diff = (u - u0)
        F_u_tent = ((1/k) * inner(v, u_diff) * dx()
                    + inner(v, grad(u0)*u0) * dx()
                    + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx()
                    - nu * inner(grad(u_mean).T*n, v) * ds()
                    + inner(v, p0*n) * ds()
                    - inner(v, f) * dx())

        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(u1)*dx()

        # Velocity correction
        a_u_corr = inner(v, u)*dx()
        L_u_corr = inner(v, u1)*dx() - k*inner(v, grad(p1-p0))*dx()

        # Assemble matrices
        A_u_tent = assemble(a_u_tent)
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        # Yield initial data for postprocessing
        yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                        t=float(t), timestep=start_timestep, u=u0, p=p0)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            t.assign(timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Scale to solver pressure
            p0.vector()[:] *= 1.0/rho

            # Compute tentative velocity step
            b = assemble(L_u_tent)
            for bc in bcu:
                bc.apply(A_u_tent, b)
            A_u_tent.apply("insert")
            b.apply("insert")
            timer.completed("u1 construct rhs")

            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)
            timer.completed("u1 solve (%s, %d, %d)"%(', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            b = assemble(L_p_corr)
            if len(bcp) == 0:
                normalize(b)
            else:
                # Scale to physical pressure
                b *= rho
                for bc in bcp:
                    bc.apply(A_p_corr, b)
                A_p_corr.apply("insert")
                b.apply("insert")
                # ... and back to solver pressure
                b *= 1.0/rho
            timer.completed("p construct rhs")

            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0:
                normalize(p1.vector())
            timer.completed("p solve (%s, %d, %d)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            b = assemble(L_u_corr)
            for bc in bcu:
                bc.apply(A_u_corr, b)
            A_u_corr.apply("insert")
            b.apply("insert")
            timer.completed("u2 construct rhs")

            solver_params = self.params.solver_u_corr
            iter = solve(A_u_corr, u1.vector(), b, *solver_params)
            timer.completed("u2 solve (%s, %d)"%(', '.join(solver_params), b.size()),{"iter": iter})

            # Rotate functions for next timestep
            u0.assign(u1)
            p0.assign(p1)

            # Scale to physical pressure
            p0.vector()[:] *= rho

            # Yield data for postprocessing
            yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                            t=float(t), timestep=timestep, u=u0, p=p0, state=(u1,p1))
