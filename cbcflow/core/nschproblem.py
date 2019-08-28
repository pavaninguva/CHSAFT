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
from __future__ import division

import ufl

from cbcpost import ParamDict, Parameterized
from cbcflow.dol import Constant, MeshFunction

from common.constants import PECLET, SCHMIDT


class NSCHProblem(Parameterized):
    """Base class for all Navier-Stokes/Cahn-Hilliard problems."""

    def __init__(self, params):
        """Initialize problem instance.

        Params will be taken from default_params and overridden
        by the values given to this constructor.
        """
        Parameterized.__init__(self, params)

        # Optionally set end time based on period and number of periods
        if self.params.T is None:
            if self.params.num_periods is None or self.params.period is None:
                raise ValueError(
                    "You must provide parameter values for either end time T, or period and num_periods."
                )
            self.params.T = (
                self.params.T0 + self.params.period * self.params.num_periods
            )
        else:
            if self.params.num_periods is not None:
                raise ValueError(
                    "Ambiguous time period, cannot set both T and num_periods."
                )

    @classmethod
    def default_params(cls):
        """Returns the default parameters for a problem.

        Explanation of parameters:

        Time parameters:

          - start_timestep: int, initial time step number
          - dt: float, time discretization value
          - T0: float, initial time
          - T: float, end time
          - period: float, length of period
          - num_periods: float, number of periods to run

        Either T or period and num_period must be set.
        If T is not set, T=T0+period*num_periods is used.

        Physical parameters:

          - mu: float, kinematic viscosity
          - rho: float, mass density

        Space discretization parameters:

          - mesh_file: str, filename to load mesh from (if any)

        """
        params = ParamDict(
            # Physical parameters:
            mu=None,
            rho=None,
            # scaling parameters
            Sc=SCHMIDT,
            Pe=PECLET,
            n_chain=1000,
            coupling_multiplier=0.1,
            # Time parameters:
            start_timestep=0,
            dt=None,
            T0=0.0,
            T=None,
            period=None,
            num_periods=None,
            # Spatial discretization parameters:
            mesh_file=None,
        )
        return params

    def initialize_geometry(self, mesh, facet_domains=None, cell_domains=None):
        """Stores mesh, domains and related quantities in a canonical member naming.

        Creates attributes on self:

            - mesh
            - facet_domains
            - cell_domains
            - ds
            - dS
            - dx

        """
        # Store geometry properties
        self.mesh = mesh
        self.facet_domains = facet_domains
        self.cell_domains = cell_domains

        # Fetch domains from mesh if necessary and avaiable
        domains = mesh.domains()
        if domains is not None:
            dim = mesh.geometry().dim()
            if self.facet_domains is None:
                self.facet_domains = MeshFunction("size_t", mesh, dim - 1, domains)
            if self.cell_domains is None:
                self.cell_domains = MeshFunction("size_t", mesh, dim, domains)

        # Attach domains to measures for convenience
        self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_domains)
        self.dS = ufl.Measure("dS", domain=self.mesh, subdomain_data=self.facet_domains)
        self.dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_domains)

    def observations(self, spaces, t):
        """Return observations of velocity for optimization problem.

        Optimization problem support is currently experimental.
        Can be ignored for non-control problems.

        TODO: Document expected observations behaviour here.
        """
        return []

    def controls(self, spaces):
        """Return controls for optimization problem.

        Optimization problem support is currently experimental.
        Can be ignored for non-control problems.

        TODO: Document expected controls behaviour here.
        """
        return []

    def cost_functionals(self, spaces, t, observations, controls):
        """Return cost functionals for optimization problem.

        Optimization problem support is currently experimental.
        Can be ignored for non-control problems.

        TODO: Document expected cost functionals behaviour here.
        """
        return []

    def __kinematic_viscosity(self, controls):  # TODO: Enable and use this
        """Return the kinematic viscosity nu=mu/rho."""
        return Constant(self.params.mu / self.params.rho, name="kinematic_viscosity")

    def __dynamic_viscosity(self, controls):  # TODO: Enable and use this
        """Return the dynamic viscosity mu."""
        return Constant(self.params.mu, name="dynamic_viscosity")

    def __density(self, controls):  # TODO: Enable and use this
        """Return the density rho."""
        return Constant(self.params.rho, name="density")

    def initial_conditions(self, spaces, controls):
        """Return initial conditions.

        The initial conditions should be specified as follows: ::
            # Return u=(x,y,0) and p=0 as initial conditions
            u0 = [Expression("x[0]"), Expression("x[1]"), Constant(0)]
            p0 = Constant(0)
            return u0, p0

        Note that the velocity is specified as a list of scalars instead of
        vector expressions.

        This function must be overridden py subclass.

        Returns: u, p
        """
        raise NotImplementedError("initial_conditions must be overridden in subclass")

    def boundary_conditions(self, spaces, u, p, t, controls):
        """Return boundary conditions in raw format.

        Boundary conditions should . The boundary conditions
        can be specified as follows: ::

            # Specify u=(0,0,0) on mesh domain 0 and u=(x,y,z) on mesh domain 1
            bcu = [
                ([Constant(0), Constant(0), Constant(0)], 0),
                ([Expression("x[0]"), Expression("x[1]"), Expression("x[2]")], 1)
                ]

            # Specify p=x^2+y^2 on mesh domain 2 and p=0 on mesh domain 3
            bcp = [
                (Expression("x[0]*x[0]+x[1]*x[1]"), 2),
                (Constant(0), 3)
            ]

            return bcu, bcp

        Note that the velocity is specified as a list of scalars instead of
        vector expressions.

        For schemes applying Dirichlet boundary conditions, the domain
        argument(s) are parsed to DirichletBC and can be specified in a matter
        that matches the signature of this class.

        This function must be overridden py subclass.

        Returns: a tuple with boundary conditions for velocity and pressure

        """
        raise NotImplementedError("boundary_conditions must be overridden in subclass")

    def body_force(self, spaces, t):
        """ Return body force, defaults to 0.

        If not overridden by subclass this function will return zero.

        Returns: list of scalars.
        """
        d = self.mesh.geometry().dim()
        return [Constant(0.0, name="body_force%d" % i) for i in range(d)]

    # TODO: Add body_force here, maybe also viscosity? Maybe time to generalize?
    def update(
        self,
        spaces,
        u,
        p,
        ch,
        t,
        timestep,
        boundary_conditions,
        observations=None,
        controls=None,
        cost_functionals=None,
    ):
        """Update functions previously returned to new timestep.

        This function is called before computing the solution at a new timestep.

        The arguments boundary_conditions, observations, controls should be the
        exact lists of objects returned by boundary_conditions, observations, controls.

        Typical usage of this function would be to update time-dependent boundary
        conditions: ::

            bcu, bcp = boundary_conditions
            for bc, _ in bcu:
                bc.t = t

            for bc, _ in bcp:
                bc.t = t

        returns None

        """
        pass

    def analytical_solution(self, spaces, t):
        """Return analytical solution.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        TODO: Document expected analytical_solution behaviour here.

        Returns: u, p
        """
        raise NotImplementedError(
            "analytical_solution must be overridden in problem subclass to use analytical solution fields"
        )

    def test_functionals(self, spaces):
        """Return fields to be used by regression tests.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        Returns: list of fields.
        """
        return []

    def test_references(self):
        """Return reference values corresponding to test_functionals to be used by regression tests.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        Returns: list of reference values.
        """
        return []
