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
This scheme follows the same logic as in :class:`.IPCS_Naive`, but with a few notable exceptions.

A parameter :math:`\theta` is added to the diffusion and convection terms,
allowing for different evaluation of these, and the convection is handled semi-implicitly:

.. math::
    \frac{1}{\Delta t}\left( \tilde{u}^{n+1}-u^{n} \right)-
    \nabla\cdot\nu\nabla \tilde{u}^{n+\theta}+
    u^*\cdot\nabla \tilde{u}^{n+\theta}+\nabla p^{n}=f^{n+1},

where

.. math::
    u^* = \frac{3}{2}u^n - \frac{1}{2}u^{n-1}, \\
    \tilde{u}^{n+\theta} = \theta \tilde{u}^{n+1}+\left(1-\theta\right)u^n.

This convection term is unconditionally stable, and with :math:`\theta=0.5`,
this equation is second order in time and space [1]_.


In addition, the solution process is significantly faster by solving for each of the
velocity components separately, making for D number of smaller linear systems compared
to a large system D times the size.



.. [1] Simo, J. C., and F. Armero. *Unconditional stability and long-term behavior
    of transient algorithms for the incompressible Navier-Stokes and Euler equations.*
    Computer Methods in Applied Mechanics and Engineering 111.1 (1994): 111-154.

"""

from __future__ import division

from cbcpost.utils import cbc_log

from cbcflow.core.nsscheme import *

from cbcflow.schemes.utils import (compute_regular_timesteps,
                                   assign_ics_segregated,
                                   make_segregated_velocity_bcs,
                                   make_pressure_bcs,
                                   NSSpacePoolSegregated,
                                   RhsGenerator, create_solver)

def _get_weighted_gradient(mesh, dims, v,p):
    from fenicstools.WeightedGradient import compiled_gradient_module, weighted_gradient_matrix
    DG = FunctionSpace(mesh, "DG", 0)
    q = TestFunction(DG)
    A = assemble(TrialFunction(DG)*v*dx)
    dg = Function(DG)

    dPdX = []
    for d in dims:
        dP = assemble(p.dx(d)*q*dx)
        dPmat = as_backend_type(dP).mat()
        compiled_gradient_module.compute_DG0_to_CG_weight_matrix(A, dg)
        Amat = as_backend_type(A).mat()

        Cpmat = Amat.matMultSymbolic(dPmat)
        Amat.matMultNumeric(dPmat, Cpmat)

        # Perform some strange copies that apparently saves memory
        Cpmat2 = Cpmat.copy()
        Cpmat.destroy()
        Cp = PETScMatrix(Cpmat2)
        Cp = PETScMatrix(Cp)

        dPdX.append(Cp)

        MPI.barrier(mpi_comm_world())
        # Destroy petsc4py objects
        dPmat.destroy()
        Amat.destroy()
        Cpmat2.destroy()

    return dPdX

def update_extrapolation(u_est, u1, u0, theta):
    theta = float(theta)
    assert theta <= 0, "Non-linear schemes not supported"
    for d in range(len(u1)):
        u_est[d].vector().zero()
        u_est[d].vector().axpy(1-theta, u1[d].vector())
        u_est[d].vector().axpy(theta, u0[d].vector())

    return u_est

class IPCS(NSScheme):
    "Incremental pressure-correction scheme, fast and stable version."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            alpha = 0.5,
            theta = -0.5,
            rebuild_prec_frequency = 1e16,
            u_tent_prec_structure = "same_nonzero_pattern",
            u_tent_solver_parameters = {},
            p_corr_solver_parameters = {},
            u_corr_solver_parameters = {},
            low_memory_version = False,
            store_rhs_matrix_p_corr = True,
            store_rhs_matrix_u_tent = True,
            store_rhs_matrix_u_corr = True,
            store_stiffness_matrix = True,
            isotropic_diffusion=0.,
            streamline_diffusion=0.,
            crosswind_diffusion=0.,
            #assemble_convection = "standard", # unassembled, debug
            )
        return params

    def solve(self, problem, timer):
        #parameters["krylov_solver"]["monitor_convergence"] = True

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())
        alpha = self.params.alpha

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        u = TrialFunction(U)
        p = TrialFunction(Q)

        # Functions
        u0 = as_vector([Function(U, name="u0_%d"%d) for d in dims]) # u^n
        u1 = as_vector([Function(U, name="u1_%d"%d) for d in dims]) # u^{n+1}
        u_est = as_vector([Function(U, name="u_est_%d"%d) for d in dims]) # Velocity interpolation used in convection

        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims:
            u1[d].assign(u0[d])

        # Update extrapolation term for first timestep
        update_extrapolation(u_est, u1, u0, self.params.theta)

        #for d in dims: u2[d].assign(u1[d])
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        rho = float(problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        timer.completed("create function spaces, functions and boundary conditions")

        # Create forms for LHS of tenative velocity
        # Convection linearized as in Simo/Armero (1994)
        # Will set u_est = (1-theta)*u1[r] + theta*u0[r]
        a1 = (1/k) * inner(v, u) * dx()
        a2 = inner(grad(v), nu*grad(u)) * dx()
        a_conv  = inner(v, dot(u_est, nabla_grad(u)))*dx()

        # Optional stabilization
        h = CellSize(mesh)
        #common_stab = h*sqrt(inner(u_est,u_est))
        isotropic_diffusion = h*sqrt(inner(u_est, u_est))*inner(grad(u), grad(v))*dx()
        streamline_diffusion = h/(sqrt(inner(u_est, u_est))+h)*inner(dot(u_est, grad(u)), dot(u_est, grad(v)))*dx()
        crosswind_diffusion = isotropic_diffusion-streamline_diffusion

        a_stab = int(self.params["isotropic_diffusion"]!=0)*Constant(self.params["isotropic_diffusion"])*isotropic_diffusion
        a_stab += int(self.params["streamline_diffusion"]!=0)*Constant(self.params["streamline_diffusion"])*streamline_diffusion
        a_stab += int(self.params["crosswind_diffusion"]!=0)*Constant(self.params["crosswind_diffusion"])*crosswind_diffusion

        # Create the static part of the coefficient matrix for the tentative
        # velocity step. The convection is added in the time loop. We use a
        # single matrix for all dimensions, which means that the BCs must match
        # (they must apply to the full vector space at the vertex, not a
        # subspace.
        A_u_tent = assemble(a1)
        if not self.params.low_memory_version and self.params.store_stiffness_matrix:
            A = assemble(a2)
            K_conv = assemble(a_conv+a_stab)
        else:
            A = assemble(a2+a_conv+a_stab)

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*len(dims)
        if not self.params.low_memory_version and self.params.store_rhs_matrix_u_tent:
            for d in dims:
                C = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())
                rhs_u_tent[d] = RhsGenerator(U)
                rhs_u_tent[d] += A_u_tent, u0[d]
                rhs_u_tent[d] += C, p0
        else:
            rhs_u_tent = [lambda: A_u_tent*u0[d].vector()+assemble(-v*p0*n[d]*ds() + v.dx(d)*p0*dx()) for d in dims]

        # Tentative velocity solver
        solver_u_tent = create_solver(*self.params.solver_u_tent)
        if 'preconditioner' in solver_u_tent.parameters:
            solver_u_tent.parameters['preconditioner']['structure'] = 'same'

        timer.completed("create tenative velocity solver")

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        if not self.params.low_memory_version and self.params.store_rhs_matrix_p_corr:
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            Ku = [None]*len(dims)
            for d in dims:
                Ku[d] = assemble(-(1/k)*q*u.dx(d)*dx()) # TODO: Store forms in list, this is copied below
                rhs_p_corr += Ku[d], u1[d]
        else:
            rhs_p_corr = lambda: A_p_corr*p0.vector() + assemble(-(1/k)*q*sum(u1[d].dx(d) for d in dims)*dx())

        # Pressure correction solver
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        for bc in bcp:
            bc.apply(A_p_corr)
        A_p_corr.apply("insert")

        solver_p_corr = create_solver(*solver_p_params)
        solver_p_corr.set_operator(A_p_corr)
        if 'preconditioner' in solver_p_corr.parameters:
                solver_p_corr.parameters['preconditioner']['structure'] = 'same'
        solver_p_corr.parameters.update(self.params.p_corr_solver_parameters)

        timer.completed("create pressure correction solver")

        # Velocity correction solver
        if self.params.solver_u_corr not in ["WeightedGradient"]:
            # Velocity correction. Like for the tentative velocity, a single LHS is used.
            M = assemble(inner(u,v)*dx())
            A_u_corr = assemble(inner(u,v)*dx())
            if not self.params.low_memory_version and self.params.store_rhs_matrix_u_corr:
                rhs_u_corr = [None]*len(dims)
                Kp = [None]*len(dims)
                for d in dims:
                    Kp[d] = assemble(-k*inner(v, grad(p)[d])*dx())
                    rhs_u_corr[d] = RhsGenerator(U)
                    rhs_u_corr[d] += M, u1[d]
                    rhs_u_corr[d] += Kp[d], p1
                    rhs_u_corr[d] -= Kp[d], p0
            else:
                rhs_u_corr = [lambda: A_u_corr*u1[d].vector()+assemble(-k*inner(v, grad(p1-p0)[d])*dx()) for d in dims]

            # Apply BCs to LHS
            for bc in bcu:
                bc[0].apply(A_u_corr)
            A_u_corr.apply("insert")

            solver_u_corr = create_solver(*self.params.solver_u_corr)
            solver_u_corr.set_operator(A_u_corr)
            if 'preconditioner' in solver_u_corr.parameters:
                solver_u_corr.parameters['preconditioner']['structure'] = 'same'
            solver_u_corr.parameters.update(self.params.u_corr_solver_parameters)

        elif self.params.solver_u_corr == "WeightedGradient":
            assert self.params.p_degree == 1
            dPdX = _get_weighted_gradient(mesh, dims, v,p)

        timer.completed("create velocity correction solver")

        # Yield initial conditions
        yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                        t=float(t), timestep=start_timestep, u=u1, p=p1)
        timer.completed("initial postprocessor update")

        # Time loop
        for timestep in xrange(start_timestep+1,len(timesteps)):
            t.assign(timesteps[timestep])

            # Update various functions
            problem.update(spaces, u1, p1, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            p0.vector()[:] *= 1.0/rho
            p1.vector()[:] *= 1.0/rho

            # Assemble convection and stabilization
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                # Assemble only convection matrix
                K_conv.zero()
                assemble(a_conv+a_stab, tensor=K_conv)
                A_u_tent.axpy(-(1.0-alpha), K_conv, True)
            else:
                # Assemble convection, stabilization, and diffusion matrix in one
                assemble(a2+a_conv+a_stab, tensor=A)

            A_u_tent.axpy(-(1.0-alpha), A, True)
            timer.completed("built A_u_tent for rhs")

            # Use A_u_tent in current form to create rhs
            # Note: No need to apply bcs to A_u_tent (this is set directly on b)
            b = [None]*len(dims)
            for d in dims:
                b[d] = rhs_u_tent[d]()
                for bc in bcu:
                    bc[d].apply(b[d])
            for d in dims: b[d].apply("insert")
            timer.completed("built tentative velocity rhs")

            # Construct lhs for tentative velocity
            A_u_tent.axpy(1.0, A, True)
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                A_u_tent.axpy(1.0, K_conv, True)
            for bc in bcu:
                bc[0].apply(A_u_tent)
            A_u_tent.apply("insert")

            timer.completed("u_tent construct lhs")

            # Check if preconditioner is to be rebuilt
            if timestep % self.params.rebuild_prec_frequency == 0 and 'preconditioner' in solver_u_tent.parameters:
                solver_u_tent.parameters['preconditioner']['structure'] = self.params.u_tent_prec_structure

            # Compute tentative velocity step
            for d in dims:
                iter = solver_u_tent.solve(A_u_tent, u1[d].vector(), b[d])

                # Preconditioner is the same for all three components, so don't rebuild several times
                if 'preconditioner' in solver_u_tent.parameters:
                    solver_u_tent.parameters['preconditioner']['structure'] = "same"

                timer.completed("u_tent solve (%s, %d dofs)"%(', '.join(self.params.solver_u_tent), b[d].size()), {"iter": iter})



            # Reset A_u_tent to mass matrix
            A_u_tent.axpy(-alpha, A, True)
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                A_u_tent.axpy(-alpha, K_conv, True)

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0:
                normalize(b)
            b *= rho
            for bc in bcp:
                # Restore physical pressure and apply bcs
                bc.apply(b)

            # Rescale to solver pressure
            b.apply("insert")
            b *= 1.0/rho

            timer.completed("p_corr construct rhs")

            # Solve p_corr
            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0:
                normalize(p1.vector())

            timer.completed("p_corr solve (%s, %d dofs)"%(', '.join(solver_p_params), b.size()), {"iter": iter})

            if self.params.solver_u_corr not in ["WeightedGradient"]:
                # Velocity correction
                for d in dims:
                    b = rhs_u_corr[d]()
                    for bc in bcu: bc[d].apply(b)
                    b.apply("insert")
                    timer.completed("u_corr construct rhs")

                    iter = solver_u_corr.solve(u1[d].vector(), b)
                    timer.completed("u_corr solve (%s, %d dofs)"%(', '.join(self.params.solver_u_corr), b.size()),{"iter": iter})
            elif self.params.solver_u_corr == "WeightedGradient":
                for d in dims:
                    u1[d].vector().axpy(-dt, dPdX[d]*(p1.vector()-p0.vector()))
                    for bc in bcu:
                        bc[d].apply(u1[d].vector())
                    u1[d].vector().apply("insert")
                    timer.completed("u_corr solve (weighted_gradient, %d dofs)" % u1[d].vector().size())

            # Update extrapolation term for next timestep
            update_extrapolation(u_est, u1, u0, self.params.theta)

             # Rotate functions for next timestep
            for d in dims:
                u0[d].assign(u1[d])
            p0.assign(p1)

            p0.vector()[:] *= rho
            p1.vector()[:] *= rho

            # Yield data for postprocessing
            yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                            t=float(t), timestep=timestep, u=u1, p=p1, state=(u1,p1))
            timer.completed("updated postprocessing (completed timestep)")
