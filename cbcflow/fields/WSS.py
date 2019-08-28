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

from dolfin import (TestFunction, Function,  FacetNormal,
                    Constant, dot, grad, ds, assemble, inner, dx,
                    TrialFunction, LinearSolver)

from cbcpost import Field, SpacePool
from cbcpost.utils import (mesh_to_boundarymesh_dofmap, cbc_warning,
                           get_set_vector)
import numpy as np
from cbcflow.fields.DynamicViscosity import DynamicViscosity
from cbcflow.schemes.utils.solver_creation import create_solver

class WSS(Field):

    def add_fields(self):
        return [DynamicViscosity()]

    def before_first_compute(self, get):
        u = get("Velocity")
        V = u.function_space()

        spaces = SpacePool(V.mesh())
        degree = V.ufl_element().degree()
        
        if degree <= 1:
            Q = spaces.get_grad_space(V, shape=(spaces.d,))
        else:
            if degree > 2:
                cbc_warning("Unable to handle higher order WSS space. Using CG1.")
            Q = spaces.get_space(1,1)

        Q_boundary = spaces.get_space(Q.ufl_element().degree(), 1, boundary=True)

        self.v = TestFunction(Q)
        self.tau = Function(Q, name="WSS_full")
        self.tau_boundary = Function(Q_boundary, name="WSS")

        local_dofmapping = mesh_to_boundarymesh_dofmap(spaces.BoundaryMesh, Q, Q_boundary)
        self._keys = np.array(local_dofmapping.keys(), dtype=np.intc)
        self._values = np.array(local_dofmapping.values(), dtype=np.intc)
        self._temp_array = np.zeros(len(self._keys), dtype=np.float_)

        Mb = assemble(inner(TestFunction(Q_boundary), TrialFunction(Q_boundary))*dx)
        self.solver = create_solver("gmres", "jacobi")
        self.solver.set_operator(Mb)

        self.b = Function(Q_boundary).vector()

        self._n = FacetNormal(V.mesh())

    def compute(self, get):
        u = get("Velocity")
        mu = get("DynamicViscosity")
        if isinstance(mu, (float, int)):
            mu = Constant(mu)

        n = self._n
        T = -mu*dot((grad(u) + grad(u).T), n)
        Tn = dot(T, n)
        Tt = T - Tn*n

        tau_form = dot(self.v, Tt)*ds()
        assemble(tau_form, tensor=self.tau.vector())

        #self.b[self._keys] = self.tau.vector()[self._values] # FIXME: This is not safe!!!
        get_set_vector(self.b, self._keys, self.tau.vector(), self._values, self._temp_array)

        # Ensure proper scaling
        self.solver.solve(self.tau_boundary.vector(), self.b)

        return self.tau_boundary
