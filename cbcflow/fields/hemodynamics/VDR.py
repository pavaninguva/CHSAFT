from dolfin import assemble, inner, dx, Constant, grad
from cbcpost import Field, MeshPool, Restrict, TimeAverage
from cbcpost.utils import create_submesh

def epsilon(u):
    return 0.5*(grad(u)+grad(u).T)

class VDR(Field):
    def __init__(self, aneurysm, near_vessel, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        
        self.aneurysm = aneurysm
        self.near_vessel = near_vessel

        self.Vnv = assemble(Constant(1)*dx(near_vessel[1], domain=near_vessel[0].mesh(), subdomain_data=near_vessel[0]))
        self.Va = assemble(Constant(1)*dx(aneurysm[1], domain=aneurysm[0].mesh(), subdomain_data=aneurysm[0]))

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

        u = "Velocity"
        if self.params.use_timeaverage:
            u = TimeAverage(u, params=params).name

        class VD(Field):
            def __init__(self, valuename, nu, subdomain,*args, **kwargs):
                Field.__init__(self, *args, **kwargs)
                self.valuename = valuename
                mf = subdomain[0]
                idx = subdomain[1]
                self.dx = dx(idx, domain=mf.mesh(), subdomain_data=mf)
                self.vol = assemble(Constant(1)*self.dx)

                self.nu = nu

            def compute(self, get):
                u = get(self.valuename)
                if u == None:
                    return None

                return 1.0/self.vol*assemble(2*self.nu*inner(epsilon(u), epsilon(u))*self.dx)

        fa = VD(u, 1.0, self.aneurysm, params=params, label="aneurysm")
        fnv = VD(u, 1.0, self.near_vessel, params=params, label="nv")
        
        f = fa/fnv

        if not self.params.use_timeaverage:
            f = TimeAverage(f, params=params)

        self.valuename = f.name

        recorded_fields = Field.stop_recording()
        return recorded_fields

    def compute(self, get):
        vdr = get(self.valuename)
        return vdr
