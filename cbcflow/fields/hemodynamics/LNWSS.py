from cbcpost import Field, DomainAvg, TimeAverage, Magnitude, MetaField
from dolfin import Function
from numpy import log

class Logarithm(MetaField):   
    def compute(self, get):
        u = get(self.valuename)
        if u == None:
            return

        if isinstance(u, Function):
            if not hasattr(self, "u"):
                self.u = Function(u)
            self.u.vector()[:] = log(u.vector().array())
            m = min(u.vector().array()[u.vector().array()!=0])
            self.u.vector()[u.vector().array()==0] = log(m)
        else:
            self.u = log(u)
        
        return self.u


class LNWSS(Field):
    def __init__(self, aneurysm_domain, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.aneurysm_domain = aneurysm_domain
    
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        #params.replace(finalize=True)
        params.update(use_timeaverage=False,
                      debug=False)
        return params
        
    def add_fields(self):
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

        fields = [Magnitude("WSS", params=params)]

        if self.params.use_timeaverage:
            fields.append(TimeAverage("Magnitude_WSS", params=params, label="%2g-%2g" %(T0,T1)))
        
        fields.append(Logarithm(fields[-1].name, params=params))
        fields.append(DomainAvg(fields[-1].name, cell_domains=self.aneurysm_domain[0], indicator=self.aneurysm_domain[1], params=params, label="aneurysm"))
        
        if not self.params.use_timeaverage:
            fields.append(TimeAverage(fields[-1].name, params=params, label="%2g-%2g" %(T0,T1)))

        self.valuename = fields[-1].name
        for f in fields:
            print f
        return fields

    def compute(self, get):
        u = get(self.valuename)
        return u