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

from cbcflow.dol import *

from cbcpost import ParamDict, Parameterized


class CHScheme(Parameterized):
    """Base class for Cahn-Hilliard schemes."""

    def __init__(self, params=None):
        Parameterized.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Discretization parameters
            # Required to be set explicitly by scheme subclass!
            c_degree=None,
            mu_degree=None,
        )
        return params

    def solve(self, problem, timer):
        """Solve Cahn-Hilliard problem by executing scheme."""
        raise NotImplementedError("Scheme must implement solve method!")
