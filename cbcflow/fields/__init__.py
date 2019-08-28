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
"""
Basic postprocessing fields.

These fields can all be created from the postprocessor from name only.
This is useful when handling dependencies for a postprocessing field: ::

    class DummyField(Field):
        def __init__(self, field_dep):
            self.field_dep = field_dep

        def compute(self, get):
            val = get(field_dep)
            return val/2.0

If a postprocessing field depends only on basic fields to be calculated, the
dependencies will be implicitly added to the postprocessor "on the fly" from
the name alone: ::

    field = DummyField("ABasicField")
    pp = NSPostProcessor()
    pp.add_field(field) # Implicitly adds ABasicField object

For non-basic dependencies, the dependencies have to be explicitly added *before*
the field depending on it: ::

    dependency = ANonBasicField("ABasicField")
    field = DummyField(dependency.name)
    pp.add_field(dependency) # Added before field
    pp.add_field(field) # pp now knows about dependency

"""

# Fields that can be constructed just by name
basic_fields = [
    # The basic solution fields:
    "Velocity",
    "Pressure",
    "Concentration_A",
    "Concentration_B",
    "Concentration_C",
    "Potential_AB",
    "Potential_AC",
    "Potential_BC",
    "Potential_A",
    "Potential_B",
    # The basic problem parameters:
    "FluidDensity",
    "KinematicViscosity",
    "DynamicViscosity",
    # Derived fields:
    "PressureGradient",
    "VelocityCurl",
    "VelocityDivergence",
    "StrainRate",
    "Stress",
    "WSS",
    "StreamFunction",
    "LocalCfl",
    "KineticEnergy",
    "Q",
    "Delta",
    "OSI",
]

for f in basic_fields:
    exec("from cbcflow.fields.%s import %s" % (f, f))
