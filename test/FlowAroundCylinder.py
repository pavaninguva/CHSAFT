#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from cbcpost import PostProcessor

from os import path

files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_0.6k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_2k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_8k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_32k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_129k.xml.gz"),
        ]

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 10.0)

class Cylinder(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt((x[0]-2.0)**2+(x[1]-0.5)**2) < 0.12+DOLFIN_EPS)

class Wall(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], 0.0) or near(x[1], 1.0))

class FlowAroundCylinder(NSProblem):
    "Flow around a cylinder in 2D."

    @classmethod
    def default_params(cls):
        "Default parameters overwriting and adding to NSProblem.default_params()"
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=5.0,
            dt=0.1,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            )
        return params

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Load mesh
        mesh = Mesh(files[self.params.refinement_level])

        # The mesh can also be generated on the fly. This has been
        # commented out because the mesh generator is non-deterministic and thus
        # unsuitable for the test suites.
        """
        refinement_levels=[32,64,128,256,512]
        N = refinement_levels[self.params.refinement_level]
        # Create mesh
        r = Rectangle(0,0, 10, 1)
        c = Circle(2.0, 0.5, 0.12)
        mesh = Mesh(r-c, N)
        """

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        Wall().mark(facet_domains, 0)
        Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    def initial_conditions(self, spaces, controls):
        "Setting the flow at rest as initial conditions"
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        "Setting uniform velocity on inlet, p=0 on outlet and no-slip elsewhere"
        c0 = Constant(0)
        c1 = Constant(1)

        # Create inflow and no-slip boundary conditions for velocity
        inflow = ([c1, c0], 1)
        noslip = ([c0, c0], 0)

        # Create boundary conditions for pressure
        bcp0 = (c0, 2)

        # Collect and return
        bcu = [inflow, noslip]
        bcp = [bcp0]
        return (bcu, bcp)


def main():
    # Create problem and scheme instances
    problem = FlowAroundCylinder({"refinement_level": 2})
    scheme = IPCS()

    # Create postprocessor instance pointing to a case directory
    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    postprocessor = PostProcessor({"casedir": casedir})

    # Creating fields to plot and save
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        StreamFunction(plot_and_save),
        ]

    # Add fields to postprocessor
    postprocessor.add_fields(fields)

    # Create NSSolver instance and solve problem
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()


if __name__ == "__main__":
    main()
