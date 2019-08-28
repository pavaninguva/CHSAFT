from cbcpost import Field, Threshold, Magnitude, Norm, DomainAvg, DomainSD, TimeAverage, ConstantField
from dolfin import assemble, Constant, dx

class SCI(Field):
    def __init__(self, aneurysm, near_vessel, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.aneurysm = aneurysm
        self.near_vessel = near_vessel
    
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

        tau = Magnitude("WSS")
        if self.params.use_timeaverage:
            tau = TimeAverage(tau, params=params)

        threshold = DomainAvg(tau, cell_domains=self.near_vessel[0], indicator=self.near_vessel[1], label="nv") + \
                    DomainSD(tau, cell_domains=self.aneurysm[0], indicator=self.aneurysm[1], label="aneurysm")
        threshold.name = "threshold_sci_nv"
        mask = Threshold(tau, threshold, dict(threshold_by="above"))

        mf = self.aneurysm[0]
        idx = self.aneurysm[1]

        Aa = ConstantField(assemble(Constant(1)*dx(idx, domain=mf.mesh(), subdomain_data=mf)))
        Aa.name = "AneurysmArea"
        Fh = Aa*DomainAvg(tau*mask, cell_domains=mf, indicator=idx, params=params)
        Fh.name = "HighShear"
        Ah = Aa*DomainAvg(mask, cell_domains=mf, indicator=idx, params=params)+ConstantField(1e-12)
        Ah.name = "HighShearArea"
        Fa = Aa*DomainAvg(tau, cell_domains=mf, indicator=idx, params=params)
        Fa.name = "TotalShear"

        f = (Fh/Fa)/(Ah/Aa)

        if not self.params.use_timeaverage:
            f = TimeAverage(f, params=params)

        self.valuename = f.name

        fields = Field.stop_recording()
        return fields

    def compute(self, get):
        sci = get(self.valuename)
        return sci