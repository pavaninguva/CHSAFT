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

from cbcflow.dol import as_vector, project

# --- Initial condition helper functions for schemes


def assign_ics_mixed(up0, spaces, ics):
    """Assign initial conditions from ics to up0.

    up0 is a mixed function in spaces.W = spaces.V * spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    up = as_vector(list(ics[0]) + [ics[1]])
    # project(up, spaces.W, function=up0) # TODO: Can do this in fenics dev
    upp = project(up, spaces.W)
    upp.rename("icup0_projection", "icup0_projection")
    up0.assign(upp)


def assign_ics_split(u0, p0, spaces, ics):
    """Assign initial conditions from ics to u0, p0.

    u0 is a vector valued function in spaces.V and p0 is a scalar function in spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    u = as_vector(list(ics[0]))
    p = ics[1]
    u0.assign(project(u, spaces.V))  # , name="u0_init_projection"))
    p0.assign(project(p, spaces.Q))  # , name="p0_init_projection"))
    # project(u, spaces.V, function=u0) # TODO: Can do this in fenics dev
    # project(p, spaces.Q, function=p0) # TODO: Can do this in fenics dev


def assign_ics_segregated(u0, p0, spaces, ics):
    """Assign initial conditions from ics to u0[:], p0.

    u0 is a list of scalar functions each in spaces.U and p0 is a scalar function in spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """
    for d in spaces.dims:
        u0[d].assign(project(ics[0][d], spaces.U))  # , name="u0_%d_init_projection"%d))
        # project(ics[0][d], spaces.U, function=u0[d]) # TODO: Can do this in fenics dev
    p0.assign(project(ics[1], spaces.Q))  # , name="p0_init_projection"))
    # project(ics[1], spaces.Q, function=p0) # TODO: Can do this in fenics dev


def assign_ics_ch(ch0, ch, spaces, ics):
    """Assign initial conditions from ics
    """
    # [0] --> a0
    # [1] --> b0
    # [2] --> _mu0_AB
    # [3] --> _mu0_AC
    # [4] --> _mu0_BC

    ch_init = as_vector([ics[0], ics[1], ics[2], ics[3], ics[4]])
    ch0.assign(project(ch_init, spaces.CH))
    ch.assign(project(ch_init, spaces.CH))


def assign_ics_nsch_split(u0, p0, ch0, ch, spaces, ics):
    """Assign initial conditions from ics to u0[:], p0, c0, mu

    u0 is a vector valued function in spaces.V and p0 is a scalar function in spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """

    # [0] --> u0
    # [1] --> p0
    # [2] --> a0
    # [3] --> b0
    # [4] --> _mu0_AB
    # [5] --> _mu0_AC
    # [6] --> _mu0_BC

    u = as_vector(list(ics[0]))
    p = ics[1]
    u0.assign(project(u, spaces.V))  # , name="u0_init_projection"))
    p0.assign(project(p, spaces.Q))  # , name="p0_init_projection"))

    ch_init = as_vector(
        [ics[2], ics[3], ics[4], ics[5], ics[6]]  # a  # b  # mu_AB  # mu_AC  # mu_BC
    )
    ch0.assign(project(ch_init, spaces.CH))
    ch.assign(project(ch_init, spaces.CH))


def assign_ics_nsch_segregated(u0, p0, ch0, ch, spaces, ics):
    """Assign initial conditions from ics to u0[:], p0, c0, mu

    u0 is a list of scalar functions each in spaces.U and p0 is a scalar function in spaces.Q,
    while ics = (icu, icp); icu = (icu0, icu1, ...).
    """

    # [0] --> u0
    # [1] --> p0
    # [2] --> a0
    # [3] --> b0
    # [4] --> _mu0_AB
    # [5] --> _mu0_AC
    # [6] --> _mu0_BC

    for d in spaces.dims:
        u0[d].assign(project(ics[0][d], spaces.U))  # , name="u0_%d_init_projection"%d))

    p0.assign(project(ics[1], spaces.Q))  # , name="p0_init_projection"))

    ch_init = as_vector(
        [ics[2], ics[3], ics[4], ics[5], ics[6]]  # a  # b  # mu_AB  # mu_AC  # mu_BC
    )
    ch0.assign(project(ch_init, spaces.CH))
    ch.assign(project(ch_init, spaces.CH))
