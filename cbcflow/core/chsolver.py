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

from time import time

from cbcpost import ParamDict, Parameterized
from cbcpost.utils import get_memory_usage, cbc_print, Timer, time_to_string
from cbcpost import Restart

from cbcflow.dol import parameters, Mesh, MPI, mpi_comm_world

from cbcflow.fields.converters_nsch import (
    ConcentrationConverter_A,
    ConcentrationConverter_B,
)
from cbcflow.fields.converters_nsch import (
    PotentialConverter_AB,
    PotentialConverter_AC,
    PotentialConverter_BC,
)


class CHSolver(Parameterized):
    """High level Cahn-Hilliard solver. This handles all logic between the cbcflow
    components.

    For full functionality, the user should instantiate this class with a CHProblem
    instance, CHScheme instance and cbcpost.PostProcessor instance.
    """

    def __init__(self, problem, scheme, postprocessor=None, params=None):
        Parameterized.__init__(self, params)

        self.problem = problem
        self.scheme = scheme
        self.postprocessor = postprocessor

        self.concentration_converter_a = ConcentrationConverter_A()
        self.concentration_converter_b = ConcentrationConverter_B()

        self.potential_converter_ab = PotentialConverter_AB()
        self.potential_converter_ac = PotentialConverter_AC()
        self.potential_converter_bc = PotentialConverter_BC()

    @classmethod
    def default_params(cls):
        """Returns the default parameters for a problem.

        Explanation of parameters:
          - restart: bool, turn restart mode on or off
          - restart_time: float, time to search for restart data
          - restart_timestep: int, timestep to search for restart data
          - check_memory_frequency: int, timestep frequency to check memory consumption
          - timer_frequency: int, timestep frequency to print more detailed timing
          - enable_annotation: bool, enable annotation of solve with dolfin-adjoint

        If restart=True, maximum one of restart_time and restart_timestep can be set.
        """
        params = ParamDict(
            restart=False,
            restart_time=-1.0,
            restart_timestep=-1,
            timer_frequency=0,
            check_memory_frequency=0,
            enable_annotation=False,
        )
        return params

    def solve(self):
        """Handles top level logic related to solve.

        Cleans casedir or loads restart data, stores parameters and mesh in
        casedir, calls scheme.solve, and lets postprocessor finalize all
        fields.

        Returns: namespace dict returned from scheme.solve
        """
        # Preserve 'solver.solve()' interface
        for data in self.isolve():
            pass
        return data

    def isolve(self):
        "Experimental iterative version of solve()."
        # Initialize solver
        self._reset()

        # Initialize postprocessor
        if self.postprocessor is not None:
            self._init_postprocessor()

        # Loop over scheme steps
        for data in self.scheme.solve(self.problem, self.timer):
            data.setdefault("d", None)
            self.update(data.ch, data.d, data.t, data.timestep, data.spaces)
            yield data

        # Finalize postprocessor
        if self.postprocessor is not None:
            self.postprocessor.finalize_all()

        # Finalize solver
        self._summarize()

    def _reset(self):
        "Reset timers and memory usage"
        self.timer = Timer(self.params.timer_frequency)
        self.timer._N = -1

        self._initial_time = time()
        self._time = time()
        self._accumulated_time = 0

        if self.params.check_memory_frequency > 0:
            self._initial_memory = get_memory_usage()

    def _init_postprocessor(self):
        self.postprocessor._timer = self.timer

        # Handle restarting or cleaning of fresh casedir
        if self.params.restart:
            Restart(
                self.problem,
                self.postprocessor,
                self.params.restart_time,
                self.params.restart_timestep,
            )
            self.timer.completed("set up restart")
        else:
            # If no restart, remove any existing data coming from cbcflow
            self.postprocessor.clean_casedir()
            self.timer.completed("cleaned casedir")

        # Store parameters
        params = ParamDict(
            solver=self.params,
            problem=self.problem.params,
            scheme=self.scheme.params,
            postprocessor=self.postprocessor.params,
        )

        self.postprocessor.store_params(params)

        # Store mesh
        assert hasattr(self.problem, "mesh") and isinstance(
            self.problem.mesh, Mesh
        ), "Unable to find problem.mesh!"
        self.postprocessor.store_mesh(self.problem.mesh)

        self.timer.completed("stored mesh")

    def _summarize(self):
        "Summarize time spent and memory (if requested)"
        final_time = time()
        msg = "Total time spent in CHSolver: %s" % time_to_string(
            final_time - self._initial_time
        )
        cbc_print(msg)

        if self.params.check_memory_frequency > 0:
            final_memory = get_memory_usage()
            msg = "Memory usage before solve: %s\nMemory usage after solve: %s" % (
                self._initial_memory,
                final_memory,
            )
            cbc_print(msg)

    def _update_timing(self, timestep, t, time_at_top):
        "Update timing and print solver status to screen."
        # Time since last update equals the time for scheme solve
        solve_time = time_at_top - self._time

        # Store time for next update
        self._time = time()

        # Time since time at top of update equals postproc + some update overhead
        pp_time = self._time - time_at_top

        # Accumulate time spent for time left estimation
        if (
            timestep > 1
        ):  # TODO: Call update for initial conditions, and skip timestep 0 instead for better estimates
            # (skip first step where jit compilation dominates)
            self._accumulated_time += solve_time
            self._accumulated_time += pp_time
        if timestep == 2:
            # (count second step twice to compensate for skipping first) (this may be overkill)
            self._accumulated_time += solve_time
            self._accumulated_time += pp_time

        # Report timing of this step
        if t:
            spent = time_to_string(self._accumulated_time)
            remaining = time_to_string(
                self._accumulated_time * (self.problem.params.T - t) / t
            )
        else:
            spent = "--"
            remaining = "--"
        msg = (
            "Timestep %5d finished (t=%.2e, %.1f%%) in %3.1fs (solve: %3.1fs). Time spent: %s Time remaining: %s"
            % (
                timestep,
                t,
                100 * t / self.problem.params.T,
                solve_time + pp_time,
                solve_time,
                spent,
                remaining,
            )
        )
        # TODO: Report to file, with additional info like memory usage, and make reporting configurable
        cbc_print(msg)

    def _update_memory(self, timestep):
        fr = self.params.check_memory_frequency
        if fr > 0 and timestep % fr == 0:
            # TODO: Report to file separately for each process
            cbc_print(
                "Memory usage is: %s" % MPI.sum(mpi_comm_world(), get_memory_usage())
            )

    def update(self, ch, d, t, timestep, spaces):
        """Callback from scheme.solve after each timestep to handle update of
        postprocessor, timings, memory etc."""
        self.timer.completed("completed solve")

        # Do not run this if restarted from this timestep
        if self.params.restart and timestep == self.problem.params.start_timestep:
            return

        # Make a record of when update was called
        time_at_top = time()

        # Disable dolfin-adjoint annotation during the postprocessing
        if self.params.enable_annotation:
            stop_annotating = parameters["adjoint"]["stop_annotating"]
            parameters["adjoint"]["stop_annotating"] = True

        # Run postprocessor
        if self.postprocessor is not None:
            fields = {}

            # Add solution fields
            fields["Concentration_A"] = lambda: self.concentration_converter_a(
                ch, spaces
            )
            fields["Concentration_B"] = lambda: self.concentration_converter_b(
                ch, spaces
            )

            fields["Potential_AB"] = lambda: self.potential_converter_ab(ch, spaces)
            fields["Potential_AC"] = lambda: self.potential_converter_ac(ch, spaces)
            fields["Potential_BC"] = lambda: self.potential_converter_bc(ch, spaces)

            # Trigger postprocessor update
            self.postprocessor.update_all(fields, t, timestep)
        self.timer.completed("postprocessor update")

        # Check memory usage
        self._update_memory(timestep)

        # Update timing data
        self.timer.increment()
        self._update_timing(timestep, t, time_at_top)

        # Enable dolfin-adjoint annotation before getting back to solve
        if self.params.enable_annotation:
            parameters["adjoint"]["stop_annotating"] = stop_annotating
