from cbcpost import MetaField2, Field, SubFunction, ConstantField, Threshold, DomainAvg, TimeAverage, Magnitude, Dot
from cbcpost.utils.slice import create_slice
from dolfin import dx, assemble, Constant

class ICI(Field):
    def __init__(self, neck, pa_planes, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.neck, self.necknormal = neck[0], neck[1]
        self.pa_planes = pa_planes

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(use_timeaverage=False,
                      debug=False)
        return params

    def add_fields(self):
        Field.start_recording()

        params = self.params.copy_recursive()
        if not self.params.debug:
            params["save"] = False
            params["plot"] = False
        params.pop("debug")
        params.pop("use_timeaverage")
        params.pop("finalize")
        T0, T1 = self.params.start_time, self.params.end_time
        assert T0 != Field.default_params().start_time
        assert T1 != Field.default_params().end_time
        if self.params.use_timeaverage:
            u = TimeAverage("Velocity", params=params)
            velocity = u.name
        else:
            velocity = "Velocity"
        

        uneck = SubFunction(velocity, self.neck, params=params, label="neck")

        A = assemble(Constant(1)*dx(domain=self.neck))
        Qin = Dot(ConstantField(self.necknormal), uneck, params=params)

        t = Threshold(Qin, ConstantField(0), dict(threshold_by="above"))
        t.params.update(params)
        t.name = "threshold_neck"

        Ain = A*DomainAvg(t)
        Ain.name = "Ain"
        Ain.params.update(params)
        Qin = A*DomainAvg(Dot(Qin,t, params=params))
        Qin.name = "Qin"
        Qin.params.update(params)
        
        Qpa = 0
        for i, (plane, n) in enumerate(self.pa_planes):
            upa = SubFunction(velocity, plane, params=params, label="pa_%d" %i)
            
            Ai = assemble(Constant(1)*dx(domain=plane))
            Q = Ai*DomainAvg(Dot(ConstantField(n), upa, params=params))
            Q.name = "Qpa%d" %i
            Q.params.update(params)
            Qpa += Magnitude(Q)
            Qpa.params.update(params)
        Qpa.name = "Sum_Qpa"
        f = (Qin/Qpa)/(Ain/A)

        if not self.params.use_timeaverage:
            f = TimeAverage(f.name, params=params)
        
        self.valuename = f.name
        
        self.Qin = Qin.name
        self.Qpa = Qpa.name
        self.Ain = Ain.name
        self.A = A
        

        fields = Field.stop_recording()

        return fields

    def compute(self, get):
        ici = get(self.valuename)
        return ici