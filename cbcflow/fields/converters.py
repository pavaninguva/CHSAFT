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

from dolfin import Function, FunctionAssigner
from cbcflow.schemes.utils import NSSpacePoolMixed, NSSpacePoolSegregated


class PressureConverter:
    def __call__(self, p, spaces):
        if not isinstance(p, Function):
            if not hasattr(self, "_p"):
                self._p = Function(spaces.Q)
                assert isinstance(spaces, NSSpacePoolMixed)
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))

            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

            p = self._p

        assert isinstance(p, Function)
        return p


class VelocityConverter:
    def __call__(self, u, spaces):
        if not isinstance(u, Function):
            d = spaces.d
            if not hasattr(self, "_u"):
                self._u = Function(spaces.V)

                if isinstance(spaces, NSSpacePoolMixed):
                    self._assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
                elif isinstance(spaces, NSSpacePoolSegregated):
                    self._assigner = FunctionAssigner(spaces.V, [spaces.U] * d)
                else:
                    error(
                        "It doesnt make sense to create a function assigner for a split space."
                    )

            if isinstance(spaces, NSSpacePoolMixed):
                # Hack: u is a ListTensor([Indexed(Coefficient()),...]),
                # get the underlying mixed function
                w = u.operands()[0].operands()[0]
                assert w.shape() == (d + 1,)
                us = w.sub(0)

            elif isinstance(spaces, NSSpacePoolSegregated):
                us = [u[i] for i in range(d)]

            self._assigner.assign(self._u, us)
            u = self._u

        assert isinstance(u, Function)
        return u
