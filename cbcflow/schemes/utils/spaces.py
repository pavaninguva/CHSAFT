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


from cbcflow.dol import (
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    BoundaryMesh,
    dolfin_version,
    MixedElement,
)
from cbcpost import SpacePool
from distutils.version import LooseVersion


def galerkin_family(degree):
    return "CG" if degree > 0 else "DG"


def decide_family(family, degree):
    return galerkin_family(degree) if family == "auto" else family


class NSSpacePool:
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."

    def __init__(self, mesh, u_degree, p_degree, u_family="auto", p_family="auto"):
        self.spacepool = SpacePool(mesh)
        assert isinstance(u_degree, int)
        assert isinstance(p_degree, int)
        assert isinstance(u_family, str)
        assert isinstance(p_family, str)
        self.u_degree = u_degree
        self.p_degree = p_degree
        self.u_family = u_family
        self.p_family = p_family
        self._spaces = {}

        # Get dimensions for convenience
        cell = mesh.ufl_cell()
        self.gdim = cell.geometric_dimension()
        self.tdim = cell.topological_dimension()
        self.gdims = range(self.gdim)
        self.tdims = range(self.tdim)

        # For compatibility, remove when code has been converted
        self.d = self.gdim
        self.dims = self.gdims

    @property
    def U(self):
        "Scalar valued space for velocity components."
        return self.spacepool.get_space(self.u_degree, 0, family=self.u_family)

    @property
    def V(self):
        "Vector valued space for velocity vector."
        return self.spacepool.get_space(self.u_degree, 1, family=self.u_family)

    @property
    def Q(self):
        "Scalar valued space for pressure."
        return self.spacepool.get_space(self.p_degree, 0, family=self.p_family)

    @property
    def DU0(self):
        "Scalar valued space for gradient component of single velocity component."
        return self.spacepool.get_space(self.u_degree - 1, 0, family=self.u_family)

    @property
    def DU(self):
        "Vector valued space for gradients of single velocity components."
        return self.spacepool.get_space(self.u_degree - 1, 1, family=self.u_family)

    @property
    def DV(self):
        "Tensor valued space for gradients of velocity vector."
        return self.spacepool.get_space(self.u_degree - 1, 2, family=self.u_family)

    @property
    def DQ0(self):
        "Scalar valued space for pressure gradient component."
        return self.spacepool.get_space(self.p_degree - 1, 0, family=self.p_family)

    @property
    def DQ(self):
        "Vector valued space for pressure gradient."
        return self.spacepool.get_space(self.p_degree - 1, 1, family=self.p_family)

    @property
    def W(self):
        "Mixed velocity-pressure space."
        space = self._spaces.get("W")
        if space is None:
            if LooseVersion(dolfin_version()) > LooseVersion("1.6.0"):
                space = FunctionSpace(
                    self.spacepool.mesh,
                    MixedElement(self.V.ufl_element(), self.Q.ufl_element()),
                )
            else:
                space = self.V * self.Q
            self._spaces["W"] = space
        return space


class NSSpacePoolMixed(NSSpacePool):
    "A function space pool with custom named spaces for use with mixed Navier-Stokes schemes."
    pass


class NSSpacePoolSplit(NSSpacePool):
    "A function space pool with custom named spaces for use with split Navier-Stokes schemes."
    pass


class NSSpacePoolSegregated(NSSpacePool):
    "A function space pool with custom named spaces for use with segregated Navier-Stokes schemes."
    pass


class NSCHSpacePool:
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."

    def __init__(
        self,
        mesh,
        u_degree,
        p_degree,
        c_degree,
        mu_degree,
        u_family="auto",
        p_family="auto",
        c_family="auto",
        mu_family="auto",
        constrained_domain=None,
    ):

        self.spacepool = SpacePool(mesh)
        assert isinstance(u_degree, int)
        assert isinstance(p_degree, int)
        assert isinstance(c_degree, int)
        assert isinstance(mu_degree, int)
        assert isinstance(u_family, str)
        assert isinstance(p_family, str)
        assert isinstance(c_family, str)
        assert isinstance(mu_family, str)
        self.u_degree = u_degree
        self.p_degree = p_degree
        self.c_degree = c_degree
        self.mu_degree = mu_degree
        self.u_family = u_family
        self.p_family = p_family
        self.c_family = c_family
        self.mu_family = mu_family
        self._spaces = {}

        # Get dimensions for convenience
        cell = mesh.ufl_cell()
        self.gdim = cell.geometric_dimension()
        self.tdim = cell.topological_dimension()
        self.cdim = 1  # dimension of concentration vector

        self.gdims = range(self.gdim)
        self.tdims = range(self.tdim)
        self.cdims = range(self.cdim)

        # For compatibility, remove when code has been converted
        self.d = self.gdim
        self.dims = self.gdims

        # For periodic boundary initial_conditions
        self.constrained_domain = constrained_domain

    @property
    def U(self):
        "Scalar valued space for velocity components."
        return self.spacepool.get_space(self.u_degree, 0, family=self.u_family)

    @property
    def V(self):
        "Vector valued space for velocity vector."
        return self.spacepool.get_space(self.u_degree, 1, family=self.u_family)

    @property
    def Q(self):
        "Scalar valued space for pressure."
        return self.spacepool.get_space(self.p_degree, 0, family=self.p_family)

    @property
    def C(self):
        "Scalar valued space for concentration."
        return self.spacepool.get_space(self.c_degree, 0, family=self.c_family)

    @property
    def M(self):
        "Scalar valued space for potential."
        return self.spacepool.get_space(self.mu_degree, 0, family=self.mu_family)

    @property
    def DU0(self):
        "Scalar valued space for gradient component of single velocity component."
        return self.spacepool.get_space(self.u_degree - 1, 0, family=self.u_family)

    @property
    def DU(self):
        "Vector valued space for gradients of single velocity components."
        return self.spacepool.get_space(self.u_degree - 1, 1, family=self.u_family)

    @property
    def DV(self):
        "Tensor valued space for gradients of velocity vector."
        return self.spacepool.get_space(self.u_degree - 1, 2, family=self.u_family)

    @property
    def DQ0(self):
        "Scalar valued space for pressure gradient component."
        return self.spacepool.get_space(self.p_degree - 1, 0, family=self.p_family)

    @property
    def DC0(self):
        "Scalar valued space for concentration gradient component."
        return self.spacepool.get_space(self.c_degree - 1, 0, family=self.c_family)

    @property
    def DM0(self):
        "Scalar valued space for potential gradient component."
        return self.spacepool.get_space(self.mu_degree - 1, 0, family=self.mu_family)

    @property
    def DQ(self):
        "Vector valued space for pressure gradient."
        return self.spacepool.get_space(self.p_degree - 1, 1, family=self.p_family)

    @property
    def DC(self):
        "Vector valued space for concentration gradient."
        return self.spacepool.get_space(self.c_degree - 1, 1, family=self.c_family)

    @property
    def DM(self):
        "Vector valued space for potential gradient."
        return self.spacepool.get_space(self.mu_degree - 1, 1, family=self.mu_family)

    @property
    def W(self):
        "Mixed velocity-pressure space."
        space = self._spaces.get("W")
        if space is None:
            if LooseVersion(dolfin_version()) > LooseVersion("1.6.0"):
                space = FunctionSpace(
                    self.spacepool.mesh,
                    MixedElement(self.V.ufl_element(), self.Q.ufl_element()),
                )
            else:
                space = self.V * self.Q
            self._spaces["W"] = space
        return space

    @property
    def CH(self):
        "Mixed contentration-potential space."
        space = self._spaces.get("CH")

        if space is None:
            if LooseVersion(dolfin_version()) > LooseVersion("1.6.0"):
                space = FunctionSpace(
                    self.spacepool.mesh,
                    MixedElement(
                        self.C.ufl_element(),
                        self.C.ufl_element(),
                        self.M.ufl_element(),
                        self.M.ufl_element(),
                        self.M.ufl_element(),
                    ),
                    constrained_domain=self.constrained_domain,
                )
            else:
                # space = self.C*self.M
                raise NotImplementedError("Use dolfin version > 1.6.0")
            self._spaces["CH"] = space
        return space


class NSCHSpacePoolMixed(NSCHSpacePool):
    "A function space pool with custom named spaces for use with mixed Navier-Stokes schemes."
    pass


class NSCHSpacePoolSplit(NSCHSpacePool):
    "A function space pool with custom named spaces for use with split Navier-Stokes schemes."
    pass


class NSCHSpacePoolSegregated(NSCHSpacePool):
    "A function space pool with custom named spaces for use with segregated Navier-Stokes schemes."
    pass


class CHSpacePool:
    "A function space pool with custom named spaces for use with Cahn-Hillard schemes."

    def __init__(
        self,
        mesh,
        c_degree,
        mu_degree,
        c_family="auto",
        mu_family="auto",
        constrained_domain=None,
        configuration="c2-mu3",
    ):

        self.spacepool = SpacePool(mesh)
        assert isinstance(c_degree, int)
        assert isinstance(mu_degree, int)
        assert isinstance(c_family, str)
        assert isinstance(mu_family, str)
        self.c_degree = c_degree
        self.mu_degree = mu_degree
        self.c_family = c_family
        self.mu_family = mu_family
        self._spaces = {}
        self.configuration = configuration

        # Get dimensions for convenience
        cell = mesh.ufl_cell()
        self.gdim = cell.geometric_dimension()
        self.tdim = cell.topological_dimension()
        self.cdim = 1  # dimension of concentration vector

        self.gdims = range(self.gdim)
        self.tdims = range(self.tdim)
        self.cdims = range(self.cdim)

        # For compatibility, remove when code has been converted
        self.d = self.gdim
        self.dims = self.gdims

        # For periodic boundary initial_conditions
        self.constrained_domain = constrained_domain

    @property
    def C(self):
        "Scalar valued space for concentration."
        return self.spacepool.get_space(self.c_degree, 0, family=self.c_family)

    @property
    def M(self):
        "Scalar valued space for potential."
        return self.spacepool.get_space(self.mu_degree, 0, family=self.mu_family)

    @property
    def DC0(self):
        "Scalar valued space for concentration gradient component."
        return self.spacepool.get_space(self.c_degree - 1, 0, family=self.c_family)

    @property
    def DM0(self):
        "Scalar valued space for potential gradient component."
        return self.spacepool.get_space(self.mu_degree - 1, 0, family=self.mu_family)

    @property
    def DC(self):
        "Vector valued space for concentration gradient."
        return self.spacepool.get_space(self.c_degree - 1, 1, family=self.c_family)

    @property
    def DM(self):
        "Vector valued space for potential gradient."
        return self.spacepool.get_space(self.mu_degree - 1, 1, family=self.mu_family)

    @property
    def CH(self):
        "Mixed contentration-potential space."
        space = self._spaces.get("CH")

        if space is None:
            if LooseVersion(dolfin_version()) > LooseVersion("1.6.0"):
                if self.configuration == "c2-mu2":
                    space = FunctionSpace(
                        self.spacepool.mesh,
                        MixedElement(
                            self.C.ufl_element(),
                            self.C.ufl_element(),
                            self.M.ufl_element(),
                            self.M.ufl_element(),
                        ),
                        constrained_domain=self.constrained_domain,
                    )
                elif self.configuration == "c2-mu3":
                    space = FunctionSpace(
                        self.spacepool.mesh,
                        MixedElement(
                            self.C.ufl_element(),
                            self.C.ufl_element(),
                            self.M.ufl_element(),
                            self.M.ufl_element(),
                            self.M.ufl_element(),
                        ),
                        constrained_domain=self.constrained_domain,
                    )
                else:
                    raise NotImplementedError(
                        "Configuration {} not supported.".format(configuration)
                    )
            else:
                # space = self.C*self.M
                raise NotImplementedError("Use dolfin version > 1.6.0")
            self._spaces["CH"] = space
        return space

