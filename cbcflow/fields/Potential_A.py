from cbcpost import SolutionField


class Potential_A(SolutionField):
    def __init__(self, params=None, label=None):
        SolutionField.__init__(self, "Potential_A", params, label)
