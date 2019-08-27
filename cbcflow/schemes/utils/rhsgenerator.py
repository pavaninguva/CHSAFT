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

from cbcflow.dol import (GenericMatrix, GenericVector, Expression, Function,
                         assemble, inner, TestFunction, dx, DirichletBC)

import ufl

class RhsGenerator(object):
    """Class for storing the instructions to create the RHS vector b.
    The two main purposes of this class are:

    - make it easy to define the LHS and RHS in the same place

    - make it easy to generate RHS from matrix-XXX products, where XXX may be either
      * a Constant (which can be projected to a vector at once)
      * an Expression (which must be projected each time, because its parameters may change)
      * a Function
      """

    def __init__(self, space):
        self.space = space
        self.matvecs = []
        self.form = None
        self.vecs = []
        f = Function(self.space)
        self.b = f.vector().copy()

    def __iadd__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            assert isinstance(A, GenericMatrix)
            self.matvecs.append((A, self._as_vector_or_timedep(x), 1))
        elif isinstance(ins, GenericVector):
            self.vecs.append(ins)
        elif isinstance(ins, ufl.Form):
            if self.form is None:
                self.form = ins
            else:
                self.form += ins
        else:
            raise RuntimeError("Unknown RHS generator "+str(type(ins)))
        return self

    def __isub__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            if isinstance(A, GenericMatrix):
                self.matvecs.append((A, self._as_vector_or_timedep(x), -1))
                return self
        raise RuntimeError("Try '+=' instead")

    def _as_vector_or_timedep(self, x):
        if isinstance(x, (GenericVector, Expression, Function)):
            return x
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def _as_vector(self, x):
        if isinstance(x, GenericVector):
            return x
        if isinstance(x, Function):
            return x.vector()
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def __call__(self, bcs=None, symmetric_mod=None):
        b = self.b
        b.zero()
        for mat, x, alpha in self.matvecs:
            b_ = mat * self._as_vector(x)
            if alpha != 1:
                b_ *= alpha
            b += b_
        for vec in self.vecs:
            b += vec
        if self.form is not None:
            assemble(self.form, tensor=b, add_values=True)
        for bc in self._wrap_in_list(bcs, "bcs", DirichletBC):
            bc.apply(b)
        if symmetric_mod:
            b -= symmetric_mod*b
        return b

    def _wrap_in_list(self, obj, name, types=type):
        if obj is None:
            lst = []
        elif hasattr(obj, '__iter__'):
            lst = list(obj)
        else:
            lst = [obj]
        for obj in lst:
            if not isinstance(obj, types):
                raise TypeError("expected a (list of) %s as '%s' argument" % (str(types),name))
        return lst
