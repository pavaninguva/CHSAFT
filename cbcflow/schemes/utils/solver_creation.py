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
from dolfin import (LinearSolver, PETScKrylovSolver, PETScPreconditioner,
                    krylov_solver_methods, krylov_solver_preconditioners,
                    linear_solver_methods)
        
import petsc4py
petsc4py.init()

def create_solver(solver, preconditioner="default"):
    """Create solver from arguments. Should be flexible to handle
    
    - strings specifying the solver and preconditioner types
    - PETScKrylovSolver/PETScPreconditioner objects
    - petsc4py.PETSC.KSP/petsc4py.PETSC.pc objects
    
    or any combination of the above
    """
    # Create solver
    if isinstance(solver, str):
        try:
            linear_solvers = set(dict(linear_solver_methods()).keys())
            krylov_solvers = set(dict(krylov_solver_methods()).keys())
        except:
            linear_solvers = set(linear_solver_methods())
            krylov_solvers = set(krylov_solver_methods())
        direct_solvers = linear_solvers-krylov_solvers

        if solver in direct_solvers:
            s = LinearSolver(solver)
            return s
        elif solver in krylov_solvers:
            s = PETScKrylovSolver(solver)
        else:
            s = PETScKrylovSolver()
            s.ksp().setType(solver)
            if s.ksp().getNormType() == petsc4py.PETSc.KSP.NormType.NONE:
                s.ksp().setNormType(petsc4py.PETSc.KSP.NormType.PRECONDITIONED)
            #raise RuntimeError("Don't know how to handle solver %s" %solver)
    elif isinstance(solver, PETScKrylovSolver):
        s = solver
    elif isinstance(solver, petsc4py.PETSc.KSP):
        s = PETScKrylovSolver(solver)
    else:
        raise ValueError("Unable to create solver from argument of type %s" %type(solver))
    
    assert isinstance(s, PETScKrylovSolver)
    if preconditioner == "default":
        return s
    
    # Create preconditioner
    if preconditioner in [None, "none", "None"]:
        pc = PETScPreconditioner("none")
        pc.set(s)
        return s
    elif isinstance(preconditioner, str):
        if preconditioner in krylov_solver_preconditioners():
            pc = PETScPreconditioner(preconditioner)
            pc.set(s)
            return s
        elif preconditioner in ["additive_schwarz", "bjacobi", "jacobi"]:
            if preconditioner == "additive_schwarz":
                pc_type = "asm"
            else:
                pc_type = preconditioner
            ksp = s.ksp()
            pc = ksp.pc
            pc.setType(pc_type)
            return s
    elif isinstance(preconditioner, PETScPreconditioner):
        pc = preconditioner
        pc.set(s)
        return s
    elif isinstance(preconditioner, petsc4py.PETSc.PC):
        ksp = s.ksp()
        ksp.setPC(preconditioner)
        return s
    else:
        raise ValueError("Unable to create preconditioner from argument of type %s" %type(solver))

    raise RuntimeError("Should not reach this code. The solver/preconditioner (%s/%s) failed to return a valid solver." %(str(solver),str(preconditioner)))