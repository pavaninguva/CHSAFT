from cbcpost import Field, SubFunction, DomainAvg, TimeAverage, Magnitude

class PLC(Field):
    def __init__(self, upstream_planes, downstream_planes, rho, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.upstream_planes = upstream_planes
        self.downstream_planes = downstream_planes
        self.rho = rho

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(use_timeaverage=False,
                      debug=False)
        return params

    def add_fields(self):
        Field.start_recording()
        params = self.params.copy()
        if not self.params.debug:
            params["save"] = False
            params["plot"] = False
        params.pop("debug")
        params.pop("use_timeaverage")
        params.pop("finalize")
            
        T0, T1 = self.params.start_time, self.params.end_time
        assert T0 != Field.default_params().start_time
        assert T1 != Field.default_params().end_time

        fields = []
        dep_fields = []
        if self.params.use_timeaverage:
            f = TimeAverage("Velocity", params=params, label="%2g-%2g" %(T0,T1))
            velocity = f.name
            fields.append(f)
            
            f = TimeAverage("Pressure", params=params, label="%2g-%2g" %(T0,T1))
            pressure = f.name
            fields.append(f)
        else:
            velocity = "Velocity"
            pressure = "Pressure"

        upstream_fields = []
        dynamic_upstream = []
        for i, plane in enumerate(self.upstream_planes):
            fu = DomainAvg(Magnitude(SubFunction(velocity, plane, params=params, label="PLC_upstream_"+str(i)), params=params), params=params)
            fp = DomainAvg(SubFunction(pressure, plane, params=params, label="PLC_upstream_"+str(i)), params=params)

            f = 0.5*self.rho*fu*fu
            f.params.update(params)

            f.name = "PLC_DynamicUpstream_%d" %i
            dynamic_upstream.append(f)

            f += fp
            f.params.update(params)
            f.name = "PLC_Upstream_%d" %i

            upstream_fields.append(f)

        downstream_fields = []
        for i, plane in enumerate(self.downstream_planes):
            fu = DomainAvg(Magnitude(SubFunction(velocity, plane, params=params, label="PLC_downstream_"+str(i)), params=params), params=params)
            fp = DomainAvg(SubFunction(pressure, plane, params=params, label="PLC_downstream_"+str(i)), params=params)

            f = 0.5*self.rho*fu*fu+fp
            f.params.update(params)
            f.name = "PLC_Downstream%d" %i
            downstream_fields.append(f)

        fu = sum(upstream_fields)
        fu *= 1./len(upstream_fields)
        fu.params.update(params)
        fu.name = "PLC_Upstream_Avg"
        
        fd = sum(downstream_fields)
        fd *= 1./len(downstream_fields)
        fd.params.update(params)
        fd.name = "PLC_Downstream_Avg"
        
        fdu = sum(dynamic_upstream)
        fdu *= 1./len(dynamic_upstream)
        fdu.params.update(params)
        fdu.name = "PLC_DynamicUpstream_Avg"
        f = (fu-fd)/fdu
        f.params.update(params)
        
        if not self.params.use_timeaverage:
            f = TimeAverage(f, params=params)
        self.valuename = f.name

        recorded_fields = Field.stop_recording()

        return recorded_fields

    def compute(self, get):
        plc = get(self.valuename)
        return plc
