from cbcpost import Field, MetaField2, DomainAvg, DomainSD, TimeAverage, Magnitude, Threshold
from dolfin import Function

class LSA(Field):
    def __init__(self, aneurysm_domain, parent_artery, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.aneurysm_domain = aneurysm_domain
        self.parent_artery = parent_artery
        assert self.params.method in ["Xiang", "Cebral"]
    
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        #params.replace(finalize=True)
        params.update(use_timeaverage=True,
                      debug=False,
                      method="Xiang")
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
        params.pop("method")

        T0, T1 = self.params.start_time, self.params.end_time
        assert T0 != Field.default_params().start_time
        assert T1 != Field.default_params().end_time

        tau = Magnitude("WSS", params=params)
        if self.params.use_timeaverage:
            tau = TimeAverage(tau, params=params, label="%2g-%2g" %(T0,T1))


        if self.params.method == "Xiang":
            threshold = 0.1*DomainAvg(tau, cell_domains=self.parent_artery[0], indicator=self.parent_artery[1])
        elif self.params.method == "Cebral":
            threshold = DomainAvg(tau, cell_domains=self.parent_artery[0], indicator=self.parent_artery[1]) - \
                        DomainSD(tau, cell_domains=self.aneurysm_domain[0], indicator=self.aneurysm_domain[1])
        else:
            raise RuntimeError("Unknown method: "+str(self.params.method))

        f = DomainAvg(Threshold(tau, threshold, params=params.copy_recursive().update(threshold_by="below")),
                                cell_domains=self.aneurysm_domain[0], indicator=self.aneurysm_domain[1], params=params)

        if not self.params.use_timeaverage:
            f = TimeAverage(f, params=params, label="%2g-%2g" %(T0,T1))

        self.valuename = f.name

        fields = Field.stop_recording()
        return fields

    def compute(self, get):
        lsa = get(self.valuename)
        return lsa