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
"""A collection of Navier-Stokes schemes."""

# from cbcflow.schemes.official import official_schemes
# for f in official_schemes:
#     exec("from cbcflow.schemes.official import %s" % (f,))
#
# all_schemes = official_schemes
#
# def show_schemes():
#     "Lists which schemes are available."
#     print "Official schemes available:"
#     print "\n".join("    " + f for f in official_schemes)



from cbcflow.schemes.official import official_schemes
for f in official_schemes:
    exec("from cbcflow.schemes.official import %s" % (f,))

from cbcflow.schemes.experimental import experimental_schemes
for f in experimental_schemes:
    exec("from cbcflow.schemes.experimental import %s" % (f,))

all_schemes = official_schemes + experimental_schemes

def show_schemes():
    "Lists which schemes are available."

    print "Official schemes available:"
    print "\n".join("    " + f for f in official_schemes)

    print "Experimental schemes available:"
    print "\n".join("    " + f for f in experimental_schemes)


show_schemes()
