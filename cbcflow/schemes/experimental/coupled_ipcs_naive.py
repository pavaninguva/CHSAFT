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

from cbcflow.core.nschscheme import *

from cbcflow.schemes.utils import (
    compute_regular_timesteps,
    assign_ics_nsch_split,
    make_velocity_bcs,
    make_pressure_bcs,
    NSCHSpacePoolSplit,
)


def epsilon(u):
    "Return symmetric gradient."
    return 0.5 * (grad(u) + grad(u).T)


def sigma(u, p, mu):
    "Return stress tensor."
    return 2.0 * mu * epsilon(u) - p * Identity(len(u))


class CoupledIPCSNaive(NSCHScheme):
    "Incremental pressure-correction scheme, naive implementation."

    def __init__(self, params=None, constrained_domain=None):
        self.constrained_domain = constrained_domain
        NSCHScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSCHScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree=1,
            p_degree=1,
            c_degree=1,
            mu_degree=1,
            # theta = 0.5,
        )
        return params

    def solve(self, problem, timer):
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = NSCHSpacePoolSplit(
            mesh,
            self.params.u_degree,
            self.params.p_degree,
            self.params.c_degree,
            self.params.mu_degree,
            constrained_domain=self.constrained_domain,
        )

        V = spaces.V
        Q = spaces.Q
        CH = spaces.CH  # mixed concentratino-potential space (Cahn-Hilliard space)

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

        # Cahn-Hilliard specific details
        theta_ch = (
            1.0
        )  # time stepping family, e.g. theta_ch=1 -> backward Euler, theta_ch=0.5 -> Crank-Nicolson

        Pe = float(problem.params.Pe)
        Sc = float(problem.params.Sc)
        n_chain = float(problem.params.n_chain)
        coupling_multiplier = float(problem.params.coupling_multiplier)

        dch = TrialFunction(CH)
        h_1, h_2, j_1, j_2, j_3 = TestFunctions(CH)

        ch = Function(CH)
        ch0 = Function(CH)

        a, b, mu_AB, mu_AC, mu_BC = split(ch)
        a0, b0, mu0_AB, mu0_AC, mu0_BC = split(ch0)

        a = variable(a)
        b = variable(b)

        n_chi_AB = 6.0
        n_chi_AC = 6.0
        n_chi_BC = 6.0

        kappa_1 = (2.0 / 3.0) * n_chi_AC
        kappa_2 = (2.0 / 3.0) * n_chi_BC
        kappa_12 = (1.0 / 3.0) * (n_chi_AC + n_chi_BC - n_chi_AB)

        # "ordinary free energy" gradients
        dgda = 1.0 + ln(a) + n_chi_AB * b + n_chi_AC * (1.0 - a - b)
        dgdb = 1.0 + ln(b) + n_chi_AB * a + n_chi_BC * (1.0 - a - b)
        dgdc = 1.0 + ln(1.0 - a - b) + n_chi_AC * a + n_chi_BC * b

        # omega_(n+theta)
        mu_AB_mid = (1.0 - theta_ch) * mu0_AB + theta_ch * mu_AB
        mu_AC_mid = (1.0 - theta_ch) * mu0_AC + theta_ch * mu_AC
        mu_BC_mid = (1.0 - theta_ch) * mu0_BC + theta_ch * mu_BC

        # Use explicit velocity in advection term
        F_a = (
            a * h_1 * dx
            - a0 * h_1 * dx
            - dt * dot(a * u0, grad(h_1)) * dx
            + dt * a * b * dot(grad(mu_AB_mid), grad(h_1)) * dx
            + dt * a * (1.0 - a - b) * dot(grad(mu_AC_mid), grad(h_1)) * dx
        )

        F_b = (
            b * h_2 * dx
            - b0 * h_2 * dx
            - dt * dot(b * u0, grad(h_2)) * dx
            - dt * a * b * dot(grad(mu_AB_mid), grad(h_2)) * dx
            + dt * b * (1.0 - a - b) * dot(grad(mu_BC_mid), grad(h_2)) * dx
        )

        F_mu_AB = (
            mu_AB * j_1 * dx
            - dgda * j_1 * dx
            + dgdb * j_1 * dx
            - (kappa_1 - kappa_12) * dot(grad(a), grad(j_1)) * dx
            + (kappa_2 - kappa_12) * dot(grad(b), grad(j_1)) * dx
        )

        F_mu_AC = (
            mu_AC * j_2 * dx
            - dgda * j_2 * dx
            + dgdc * j_2 * dx
            - kappa_1 * dot(grad(a), grad(j_2)) * dx
            - kappa_12 * dot(grad(b), grad(j_2)) * dx
        )

        F_mu_BC = (
            mu_BC * j_3 * dx
            - dgdb * j_3 * dx
            + dgdc * j_3 * dx
            - kappa_2 * dot(grad(b), grad(j_3)) * dx
            - kappa_12 * dot(grad(a), grad(j_3)) * dx
        )

        F = F_a + F_b + F_mu_AB + F_mu_AC + F_mu_BC

        J_ch = derivative(F, ch, dch)

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_nsch_split(u0, p0, ch0, ch, spaces, ics)

        u1.assign(u0)
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Problem coefficients
        nu = Constant(problem.params.mu / problem.params.rho)
        rho = float(problem.params.rho)
        k = Constant(dt)
        f = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step
        u_mean = 0.5 * (u + u0)
        u_diff = u - u0
        F_u_tent = (
            (1 / k) * inner(v, u_diff) * dx()
            + inner(v, grad(u0) * u0) * dx()
            + inner(epsilon(v), sigma(u_mean, p0, nu)) * dx()
            - nu * inner(grad(u_mean).T * n, v) * ds()
            + inner(v, p0 * n) * ds()
            - inner(v, f) * dx()
        )

        a_u_tent = lhs(F_u_tent)
        L_u_tent = rhs(F_u_tent)

        # Pressure correction
        a_p_corr = inner(grad(q), grad(p)) * dx()
        L_p_corr = inner(grad(q), grad(p0)) * dx() - (1 / k) * q * div(u1) * dx()

        # Velocity correction
        a_u_corr = inner(v, u) * dx()
        L_u_corr = inner(v, u1) * dx() - k * inner(v, grad(p1 - p0)) * dx()

        # Assemble matrices
        A_u_tent = assemble(a_u_tent)
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)

        # Define non-linear solver for Cahn-Hilliard equation
        problem_ch = NonlinearVariationalProblem(F, ch, bcs=[], J=J_ch)
        solver_ch = NonlinearVariationalSolver(problem_ch)
        prm = solver_ch.parameters
        # TODO check that `newton_solver` is active

        prm["newton_solver"]["absolute_tolerance"] = 1e-5
        # prm["newton_solver"]["relative_tolerance"] = 1E-7
        # prm["newton_solver"]["maximum_iterations"] = 50
        # prm["newton_solver"]["relaxation_parameter"] = 1.0

        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        # Yield initial data for postprocessing
        yield ParamDict(
            spaces=spaces,
            observations=observations,
            controls=controls,
            t=float(t),
            timestep=start_timestep,
            u=u0,
            p=p0,
            ch=ch0,
        )

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep + 1, len(timesteps)):
            t.assign(timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, ch, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Solve Cahn-Hilliard equation
            iter = solver_ch.solve()

            # Scale to solver pressure
            p0.vector()[:] *= 1.0 / rho

            # Compute tentative velocity step
            b = assemble(L_u_tent)
            for bc in bcu:
                bc.apply(A_u_tent, b)
            A_u_tent.apply("insert")
            b.apply("insert")
            timer.completed("u1 construct rhs")

            iter = solve(A_u_tent, u1.vector(), b, *self.params.solver_u_tent)
            timer.completed(
                "u1 solve (%s, %d, %d)"
                % (", ".join(self.params.solver_u_tent), b.size(), iter)
            )

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
                b *= 1.0 / rho
            timer.completed("p construct rhs")

            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0:
                normalize(p1.vector())
            timer.completed(
                "p solve (%s, %d, %d)" % (", ".join(solver_p_params), b.size(), iter)
            )

            # Velocity correction
            b = assemble(L_u_corr)
            for bc in bcu:
                bc.apply(A_u_corr, b)
            A_u_corr.apply("insert")
            b.apply("insert")
            timer.completed("u2 construct rhs")

            solver_params = self.params.solver_u_corr
            iter = solve(A_u_corr, u1.vector(), b, *solver_params)
            timer.completed(
                "u2 solve (%s, %d)" % (", ".join(solver_params), b.size()),
                {"iter": iter},
            )

            # Rotate functions for next timestep
            u0.assign(u1)
            p0.assign(p1)
            ch0.vector()[:] = ch.vector()

            # Scale to physical pressure
            p0.vector()[:] *= rho

            # Yield data for postprocessing
            yield ParamDict(
                spaces=spaces,
                observations=observations,
                controls=controls,
                t=float(t),
                timestep=timestep,
                u=u0,
                p=p0,
                ch=ch0,
                state=(u1, p1),
            )

