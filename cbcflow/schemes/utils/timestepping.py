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

from numpy import arange
from cbcpost.utils import cbc_warning
from dolfin import error

# --- Functions for timestepping

def compute_regular_timesteps(problem):
    """Compute fixed timesteps for problem.

    The first timestep will be T0 while the last timestep will be in the interval [T, T+dt).

    Returns (dt, timesteps, start_timestep).
    """
    # Get the time range and timestep from the problem
    T0 = problem.params.T0
    T = problem.params.T
    dt = problem.params.dt
    start_timestep = problem.params.start_timestep

    # Compute regular timesteps, including T0 and T
    timesteps = arange(T0, T+dt, dt)

    if abs(dt - (timesteps[1]-timesteps[0])) > 1e-8:
        error("Computed timestep size does not match specified dt.")

    if timesteps[-1] < T - dt*1e-6:
        error("Computed timesteps range does not include end time.")

    if timesteps[-1] > T + dt*1e-6:
        cbc_warning("End time for simulation does not match end time set for problem (T-T0 not a multiple of dt).")

    if start_timestep < 0 or start_timestep >= len(timesteps):
        error("start_timestep is beyond the computed timesteps.")

    return dt, timesteps, start_timestep
