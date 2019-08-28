
from cbcpost import Magnitude, TimeAverage, Maximum, Restrict, Field
from cbcpost.utils import create_submesh

class MWSS(Field):
    def __init__(self, aneurysm_domain, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.aneurysm_domain = aneurysm_domain

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(use_timeaverage=True,
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
        
        fields = [Magnitude("WSS")]

        T0, T1 = self.params.start_time, self.params.end_time
        assert T0 != Field.default_params().start_time
        assert T1 != Field.default_params().end_time

        if self.params.use_timeaverage:
            fields.append(TimeAverage("Magnitude_WSS", params=params, label="%2g-%2g" %(T0,T1)))

        submesh = create_submesh(self.aneurysm_domain[0].mesh(), self.aneurysm_domain[0], self.aneurysm_domain[1])

        fields.append(Restrict(fields[-1].name, submesh, label="aneurysm"))
        fields.append(Maximum(fields[-1].name))

        if not self.params.use_timeaverage:
            fields.append(TimeAverage(fields[-1].name, params=params, label="%2g-%2g" %(T0,T1)))

        self.valuename = fields[-1].name
        return fields

    def compute(self, get):
        return get(self.valuename)
