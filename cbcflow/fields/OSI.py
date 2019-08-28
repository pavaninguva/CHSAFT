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

from cbcpost import Field, Magnitude, TimeIntegral
from cbcflow.fields import WSS
from dolfin import Function, project, Constant, conditional

class OSI(Field):

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(
            finalize=True,
            )
        return params

    def add_fields(self):
        params = self.params.copy_recursive()
        params["save"] = False
        params["plot"] = False
        #params["callback"] = False
        #params.pop("finalize")

        fields = []
        #fields.append(WSS(params=params))
        #return fields

        f = TimeIntegral("WSS", params=params, label="OSI")
        fields.append(f)
        fields.append(Magnitude(f, params=params))

        f = Magnitude("WSS", params=params)
        fields.append(f)
        fields.append(TimeIntegral(f, params=params, label="OSI"))

        #f = TimeIntegral("WSS", label="OSI")
        #fields.append(f)
        #fields.append(Magnitude(f))

        #f = Magnitude("WSS")
        #fields.append(f)
        #fields.append(TimeIntegral(f, label="OSI"))

        return fields

    def before_first_compute(self, get):
        tau = get("WSS")
        self.osi = Function(tau.sub(0).function_space().collapse())

    def compute(self, get):
        # Requires the fields Magnitude(TimeIntegral("WSS", label="OSI")) and
        # TimeIntegral(Magnitude("WSS"), label="OSI")
        #self.mag_ta_wss = get("Magnitude_TimeIntegral_WSS_OSI")
        #self.ta_mag_wss = get("TimeIntegral_Magnitude_WSS_OSI")
        self.mag_ta_wss = get("Magnitude_TimeIntegral_WSS-OSI")
        self.ta_mag_wss = get("TimeIntegral_Magnitude_WSS-OSI")

        if self.params.finalize:
            return None
        elif self.mag_ta_wss == None or self.ta_mag_wss == None:
            return None
        else:
            expr = conditional(self.ta_mag_wss < 1e-15,
                               0.0,
                               0.5 * (1.0 - self.mag_ta_wss / self.ta_mag_wss))
            self.osi.assign(project(expr, self.osi.function_space()))
            return self.osi

    def after_last_compute(self, get):
        self.mag_ta_wss = get("Magnitude_TimeIntegral_WSS-OSI")
        self.ta_mag_wss = get("TimeIntegral_Magnitude_WSS-OSI")
        #print self.name, " Calling after_last_compute"

        expr = conditional(self.ta_mag_wss < 1e-15,
                           0.0,
                           0.5 * (1.0 - self.mag_ta_wss / self.ta_mag_wss))
        self.osi.assign(project(expr, self.osi.function_space()))

        return self.osi
