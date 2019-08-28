from cbcpost import Field, SpacePool
from dolfin import Function


class Concentration_C(Field):
    def before_first_compute(self, get):
        _a = get("Concentration_A")
        A = _a.function_space()
        spaces = SpacePool(A.mesh())

        if self.params.expr2function == "assemble":
            print("Assemble option")
            V = spaces.get_space(0, 0)
        else:
            print("Other option")
            V = spaces.get_space(2 * A.ufl_element().degree(), 0)

        self._function = Function(V, name=self.name)

    def compute(self, get):
        _a = get("Concentration_A")
        _b = get("Concentration_B")

        # sum of molar fractions is unity, i.e. a+b+c=1
        expr = 1.0 - _a - _b

        return self.expr2function(expr, self._function)
