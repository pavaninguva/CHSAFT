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

from dolfin import Function, FunctionAssigner, Constant
from cbcflow.schemes.utils import NSCHSpacePoolMixed, NSCHSpacePoolSegregated


# REFACTOR use config "c2-mu2" parameter for setting appropriate sub(n) values


class ConcentrationConverter_A:
    def __call__(self, ch, spaces):
        # extract `a` field from `ch` Cahn-Hilliard object

        self._c = Function(spaces.C)
        self._assigner = FunctionAssigner(spaces.C, spaces.CH.sub(0))
        self._assigner.assign(self._c, ch.sub(0))

        c = self._c
        assert isinstance(c, Function)
        return c


class ConcentrationConverter_B:
    def __call__(self, ch, spaces):
        # extract `b` field from `ch` Cahn-Hilliard object

        self._c = Function(spaces.C)
        self._assigner = FunctionAssigner(spaces.C, spaces.CH.sub(1))
        self._assigner.assign(self._c, ch.sub(1))

        c = self._c
        assert isinstance(c, Function)
        return c


class ConcentrationConverter_C:
    def __call__(self, ch, spaces):
        # extract `c` field from `ch` Cahn-Hilliard object

        self._c = Function(spaces.C)
        self._assigner = FunctionAssigner(spaces.C, spaces.CH.sub(2))
        self._assigner.assign(self._c, ch.sub(2))

        c = self._c
        assert isinstance(c, Function)
        return c


class PotentialConverter_AB:
    def __call__(self, ch, spaces):
        # extract `mu_AB` field from `ch` Cahn-Hilliard object

        self._mu = Function(spaces.M)
        self._assigner = FunctionAssigner(spaces.M, spaces.CH.sub(2))
        self._assigner.assign(self._mu, ch.sub(2))

        mu = self._mu
        assert isinstance(mu, Function)
        return mu


class PotentialConverter_AC:
    def __call__(self, ch, spaces):
        # extract `mu_AC` field from `ch` Cahn-Hilliard object

        self._mu = Function(spaces.M)
        self._assigner = FunctionAssigner(spaces.M, spaces.CH.sub(3))
        self._assigner.assign(self._mu, ch.sub(3))

        mu = self._mu
        assert isinstance(mu, Function)
        return mu


class PotentialConverter_BC:
    def __call__(self, ch, spaces):
        # extract `mu_BC` field from `ch` Cahn-Hilliard object

        self._mu = Function(spaces.M)
        self._assigner = FunctionAssigner(spaces.M, spaces.CH.sub(4))
        self._assigner.assign(self._mu, ch.sub(4))

        mu = self._mu
        assert isinstance(mu, Function)
        return mu


class PotentialConverter_A:
    def __call__(self, ch, spaces):
        # extract `mu_A` field from `ch` Cahn-Hilliard object

        self._mu = Function(spaces.M)
        self._assigner = FunctionAssigner(spaces.M, spaces.CH.sub(2))
        self._assigner.assign(self._mu, ch.sub(2))

        mu = self._mu
        assert isinstance(mu, Function)
        return mu


class PotentialConverter_B:
    def __call__(self, ch, spaces):
        # extract `mu_B` field from `ch` Cahn-Hilliard object

        self._mu = Function(spaces.M)
        self._assigner = FunctionAssigner(spaces.M, spaces.CH.sub(3))
        self._assigner.assign(self._mu, ch.sub(3))

        mu = self._mu
        assert isinstance(mu, Function)
        return mu


class PressureConverter:
    def __call__(self, p, spaces):

        if not isinstance(p, Function):

            if not hasattr(self, "_p"):
                self._p = Function(spaces.Q)
                assert isinstance(spaces, NSCHSpacePoolMixed)
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

                if isinstance(spaces, NSCHSpacePoolMixed):
                    self._assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
                elif isinstance(spaces, NSCHSpacePoolSegregated):
                    self._assigner = FunctionAssigner(spaces.V, [spaces.U] * d)
                else:
                    error(
                        "It doesnt make sense to create a function assigner for a split space."
                    )

            if isinstance(spaces, NSCHSpacePoolMixed):
                # Hack: u is a ListTensor([Indexed(Coefficient()),...]),
                # get the underlying mixed function
                w = u.operands()[0].operands()[0]
                assert w.shape() == (d + 1,)
                us = w.sub(0)

            elif isinstance(spaces, NSCHSpacePoolSegregated):
                us = [u[i] for i in range(d)]

            self._assigner.assign(self._u, us)
            u = self._u

        assert isinstance(u, Function)
        return u
