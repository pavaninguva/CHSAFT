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

from cbcpost import Field, SpacePool
from dolfin import Function, grad, Constant

class Q(Field):
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.replace(
            expr2function="project", # "assemble" | "project" | "interpolate"
            )
        return params

    def before_first_compute(self, get):
        u = get("Velocity")
        U = u.function_space()
        spaces = SpacePool(U.mesh())

        if self.params.expr2function == "assemble":
            degree = 0
        elif True:
            degree = 1
        else:
            # TODO: Use accurate degree? Plotting will be projected to CG1 anyway...
            degree = 2 * (U.ufl_element().degree() - 1)

        V = spaces.get_space(degree, 0)

        self._function = Function(V, name=self.name)

    def compute(self, get):
        u = get("Velocity")

        S = (grad(u) + grad(u).T)/2
        Omega = (grad(u) - grad(u).T)/2
        expr = 0.5*(Omega**2 - S**2)

        return self.expr2function(expr, self._function)
