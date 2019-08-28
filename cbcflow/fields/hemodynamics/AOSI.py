from cbcflow.fields.OSI import OSI
from cbcpost import Field, DomainAvg

class AOSI(Field):
    def __init__(self, aneurysm_domain, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.aneurysm_domain = aneurysm_domain

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(debug=False)
        return params

    def add_fields(self):
        params = self.params.copy_recursive()
        if not self.params.debug:
            params["save"] = False
            params["plot"] = False
        params.pop("debug")
        params.pop("finalize")

        fields = [OSI(params=params)]
        fields.append(DomainAvg(fields[-1].name, cell_domains=self.aneurysm_domain[0], indicator=self.aneurysm_domain[1], params=params, label="aneurysm"))
        self.valuename = fields[-1].name

        return fields
    
    def compute(self, get):
        u = get(self.valuename)
        return u
