from __future__ import division

import csv
import os

from cbcflow.dol import solve as dolfin_solve

from cbcpost.utils import cbc_log
from cbcflow.core.chscheme import *
from cbcflow.schemes.utils import (
    compute_regular_timesteps,
    assign_ics_ch,
    make_pressure_bcs,
    CHSpacePool,
)

from common.constants import (
    chi_AB,
    chi_AC,
    chi_BC,
    D_AB,
    D_AC,
    D_BC,
    N_A,
    N_B,
    N_C,
    kappa_AA,
    kappa_BB,
    kappa_AB,
    theta_ch,
    maximum_iterations,
    relaxation_parameter,
    convergence_criterion,
    absolute_tolerance,
    relative_tolerance,
    N_SCALE_OPTION,
    D_SCALE_OPTION,
    MOLECULAR_SCALE,
    PARTICLE_SCALE,
)


class VanillaCH(CHScheme):
    "Vanilla Cahn-Hilliard formulation."

    def __init__(self, params=None, constrained_domain=None):
        self.constrained_domain = constrained_domain
        CHScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = CHScheme.default_params()
        params.update(c_degree=1, mu_degree=1, theta=-0.5)  # TODO check role of theta
        return params

    def solve(self, problem, timer):
        # parameters["krylov_solver"]["monitor_convergence"] = True

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n = FacetNormal(mesh)
        dims = range(mesh.topology().dim())

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = CHSpacePool(
            mesh,
            self.params.c_degree,
            self.params.mu_degree,
            constrained_domain=self.constrained_domain,
            configuration="c2-mu3",
        )

        CH = spaces.CH  # mixed concentration-potential space (Cahn-Hilliard space)

        # Cahn-Hilliard specific details

        # time stepping family, e.g. theta_ch=1 -> backward Euler, theta_ch=0.5 -> Crank-Nicolson

        dch = TrialFunction(CH)
        h_1, h_2, j_1, j_2, j_3 = TestFunctions(CH)

        ch = Function(CH)
        ch0 = Function(CH)

        a, b, N_mu_AB, N_mu_AC, N_mu_BC = split(ch)
        a0, b0, N_mu0_AB, N_mu0_AC, N_mu0_BC = split(ch0)

        a = variable(a)
        b = variable(b)

        N_scale_options = {"N_A": N_A, "N_B": N_B, "N_C": N_C}
        N_SCALE = N_scale_options[N_SCALE_OPTION]

        D_scale_options = {"D_AB": D_AB, "D_AC": D_AC, "D_BC": D_BC}
        D_SCALE = D_scale_options[D_SCALE_OPTION]

        MP_SCALE = MOLECULAR_SCALE / PARTICLE_SCALE

        N_MP_kappa_AA = N_SCALE * MP_SCALE ** 2 * kappa_AA
        N_MP_kappa_BB = N_SCALE * MP_SCALE ** 2 * kappa_BB
        N_MP_kappa_AB = N_SCALE * MP_SCALE ** 2 * kappa_AB

        # "ordinary free energy" gradients
        N_dgda = N_SCALE * (
            (1.0 / N_A) * (1.0 + ln(a)) + chi_AB * b + chi_AC * (1.0 - a - b)
        )
        N_dgdb = N_SCALE * (
            (1.0 / N_B) * (1.0 + ln(b)) + chi_AB * a + chi_BC * (1.0 - a - b)
        )
        N_dgdc = N_SCALE * (
            (1.0 / N_C) * (1.0 + ln(1.0 - a - b)) + chi_AC * a + chi_BC * b
        )

        N_mu_AB_mid = (1.0 - theta_ch) * N_mu0_AB + theta_ch * N_mu_AB
        N_mu_AC_mid = (1.0 - theta_ch) * N_mu0_AC + theta_ch * N_mu_AC
        N_mu_BC_mid = (1.0 - theta_ch) * N_mu0_BC + theta_ch * N_mu_BC

        # scale diffusivity
        D_AB_ = D_AB / D_SCALE
        D_AC_ = D_AC / D_SCALE
        D_BC_ = D_BC / D_SCALE

        # transport equations
        F_a = (
            a * h_1 * dx
            - a0 * h_1 * dx
            + dt * a * b * D_AB_ * dot(grad(N_mu_AB_mid), grad(h_1)) * dx
            + dt * a * (1.0 - a - b) * D_AC_ * dot(grad(N_mu_AC_mid), grad(h_1)) * dx
        )

        F_b = (
            b * h_2 * dx
            - b0 * h_2 * dx
            - dt * a * b * D_AB_ * dot(grad(N_mu_AB_mid), grad(h_2)) * dx
            + dt * b * (1.0 - a - b) * D_BC_ * dot(grad(N_mu_BC_mid), grad(h_2)) * dx
        )

        # chemical potential equations
        F_N_mu_AB = (
            N_mu_AB * j_1 * dx
            - N_dgda * j_1 * dx
            + N_dgdb * j_1 * dx
            - (N_MP_kappa_AA - N_MP_kappa_AB) * dot(grad(a), grad(j_1)) * dx
            + (N_MP_kappa_BB - N_MP_kappa_AB) * dot(grad(b), grad(j_1)) * dx
        )

        F_N_mu_AC = (
            N_mu_AC * j_2 * dx
            - N_dgda * j_2 * dx
            + N_dgdc * j_2 * dx
            - N_MP_kappa_AA * dot(grad(a), grad(j_2)) * dx
            - N_MP_kappa_AB * dot(grad(b), grad(j_2)) * dx
        )

        F_N_mu_BC = (
            N_mu_BC * j_3 * dx
            - N_dgdb * j_3 * dx
            + N_dgdc * j_3 * dx
            - N_MP_kappa_BB * dot(grad(b), grad(j_3)) * dx
            - N_MP_kappa_AB * dot(grad(a), grad(j_3)) * dx
        )

        F = F_a + F_b + F_N_mu_AB + F_N_mu_AC + F_N_mu_BC

        J_ch = derivative(F, ch, dch)

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_ch(ch0, ch, spaces, ics)

        # Define non-linear solver for Cahn-Hilliard equation
        problem_ch = NonlinearVariationalProblem(F, ch, bcs=[], J=J_ch)
        solver_ch = NonlinearVariationalSolver(problem_ch)
        parameters = solver_ch.parameters

        parameters["nonlinear_solver"] = "newton"
        parameters["newton_solver"]["absolute_tolerance"] = absolute_tolerance
        parameters["newton_solver"]["relative_tolerance"] = relative_tolerance
        parameters["newton_solver"]["maximum_iterations"] = maximum_iterations
        parameters["newton_solver"]["relaxation_parameter"] = relaxation_parameter
        parameters["newton_solver"]["convergence_criterion"] = convergence_criterion

        # dump configuration to stdout (blocks execution)
        # info(parameters, True)

        k = Constant(dt)
        f = as_vector(problem.body_force(spaces, t))

        timer.completed("create function spaces, functions and boundary conditions")

        # Yield initial conditions
        # TODO set ch instead of ch0
        yield ParamDict(
            spaces=spaces,
            observations=observations,
            controls=controls,
            t=float(t),
            timestep=start_timestep,
            ch=ch0,
        )
        timer.completed("initial postprocessor update")

        # Time loop
        gibbs_list = []
        for timestep in xrange(start_timestep + 1, len(timesteps)):

            t.assign(timesteps[timestep])

            # Solve Cahn-Hilliard equation
            iter_ = solver_ch.solve()
            # dolfin_solve(F == 0, ch, bcs=[])

            # Assembling Gradient Energy Contributions
            term_1 = assemble(Constant(0.5) * kappa_AA * dot(grad(a), grad(a)) * dx())
            term_2 = assemble(kappa_AB * dot(grad(a), grad(b)) * dx())
            term_3 = assemble(Constant(0.5) * kappa_BB * dot(grad(b), grad(b)) * dx())

            # Assembling Flory Huggins Contributions
            # Entropic Terms
            term_4 = assemble(((a * ln(a)) / N_A) * dx)
            term_5 = assemble(((b * ln(b)) / N_B) * dx)
            term_6 = assemble((((1.0 - a - b) * ln(1.0 - a - b)) / N_C) * dx)

            # Enthalpic Terms
            term_7 = assemble((chi_AB * a * b) * dx)
            term_8 = assemble((chi_AC * a * (1.0 - a - b)) * dx)
            term_9 = assemble((chi_BC * b * (1.0 - a - b)) * dx)

            gibbs = (
                term_1
                + term_2
                + term_3
                + term_4
                + term_5
                + term_6
                + term_7
                + term_8
                + term_9
            )
            gibbs_list.append(gibbs)

            fpath = "./output_RESULTS/gibbs.csv"
            headers = ["time", "gibbs"]

            # write the header row only (for the first timestep)
            if not os.path.exists(fpath):
                with open(fpath, "w") as f:
                    w = csv.DictWriter(f, headers)
                    w.writeheader()

            # append the case data
            with open(fpath, "a") as f:
                w = csv.DictWriter(f, headers)
                w.writerow({"time": float(t), "gibbs": float(gibbs)})

            # Rotate functions for next timestep
            ch0.vector()[:] = ch.vector()

            # Yield data for postprocessing
            yield ParamDict(
                spaces=spaces,
                observations=observations,
                controls=controls,
                t=float(t),
                timestep=timestep,
                ch=ch0,
            )
            timer.completed("updated postprocessing (completed timestep)")
