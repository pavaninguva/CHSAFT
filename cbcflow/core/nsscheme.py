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

class NSScheme(Parameterized):
    """Base class for all Navier-Stokes schemes.

    TODO: Clean up and document new interface.
    """

    def __init__(self, params=None):
        Parameterized.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Discretization parameters
            u_degree = None, # Require this to be set explicitly by scheme subclass!
            p_degree = None, # Require this to be set explicitly by scheme subclass!

            # TODO: These ipcs solvers are scheme specific and should maybe not be here,
            #       however they are used by most of the splitting schemes...
            # TODO: Split these into separate parameters for solver/preconditioner?
            solver_u_tent=("gmres", "default"),
            solver_p_neumann=("gmres", "amg"),
            solver_p_dirichlet=("gmres", "amg"),
            solver_p=None, # overrides neumann/dirichlet if given
            solver_u_corr=("bicgstab", "default"),

            # Timestepping method
            # FIXME: Adaptive timestepping still unsupported
            # adaptive_timestepping=False,
            )
        return params

    def solve(self, problem, timer):
        """Solve Navier-Stokes problem by executing scheme."""
        raise NotImplementedError("Scheme must implement solve method!")
