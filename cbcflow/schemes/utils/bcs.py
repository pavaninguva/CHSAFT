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

from cbcflow.dol import (FacetNormal, DirichletBC, Constant, dot, as_vector, inner,
                         MPI, mpi_comm_world, compile_extension_module, GenericMatrix,
                         GenericVector)
from numpy import array, float_
# --- Boundary condition helper functions for schemes
def _domainargs(problem, D):
    "Helper function to pass domain args if necessary."
    if isinstance(D, int):
        return (problem.facet_domains, D)
    else:
        return (D,)

def make_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcu = [DirichletBC(spaces.V.sub(d), functions[d], *_domainargs(problem, region))
           for functions, region in bcu_raw
           for d in range(len(functions))]
    return bcu

def make_mixed_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcu = [DirichletBC(spaces.W.sub(0).sub(d), functions[d], *_domainargs(problem, region))
           for functions, region in bcu_raw
           for d in range(len(functions))]
    return bcu

def make_segregated_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcu = [[DirichletBC(spaces.U, functions[d], *_domainargs(problem, region))
            for d in range(len(functions))]
           for functions, region in bcu_raw]
    return bcu

def make_pressure_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcp = [DirichletBC(spaces.Q, function, *_domainargs(problem, region))
           for function, region in bcp_raw]
    return bcp

def make_rhs_pressure_bcs(problem, spaces, bcs, v):
    bcu_raw, bcp_raw = bcs[:2]
    ds = problem.ds
    n = FacetNormal(problem.mesh)
    Lbc = -sum(dot(function*n, v)*ds(region) for (function, region) in bcp_raw)
    return Lbc
