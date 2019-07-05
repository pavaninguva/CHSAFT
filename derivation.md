# Introduction

This document serves to outline the derivation of the binary CH system whilst maintaining the general form of the CH system for ease of integrating a variety of expressions for the homogenous free energy. 

The conversion from the dimensional species transport equation to the dimensionless form and the conversion to the corresponding finite element weak forms will also be covered. 

## Dimensional binary CH 

We work in differences of chemical potentials for convenience. Based on the derivation of Naum1994 and Petr2012,
$$
\mu_{AB} = \frac{\partial g}{\partial a} - \kappa \nabla^{2}a 
$$
where $\mu_{AB}$ is the difference of chemical potentials between species A and B, $g$ is the homogenous free energy expression, $a$ is the mole fraction of species A and $\kappa$ is the gradient energy term. 

The corresponding CH species transport equation can be written as: 
$$
\frac{\partial a}{\partial t} = \nabla \cdot \big(D_{AB}a(1-a)\nabla\mu_{AB} \big)
$$
where $t$ is the dimensional time. 

The mole fraction of species B, $b$ , is inferred by a material balance constraint: 
$$
a + b = 1
$$

## Non-dimensionalization

We introduce the following scalings: 
$$
\tilde{t} = \frac{D_{AB}t}{R_{g}^{2}}
$$
where $\tilde {t}$ is dimensionless time, $N_A$ is the length of the polymer, $R_g$ is the radius of gyration of the polymer. 
$$
\bold{x} = \tilde{\bold{x}}R_{g}
$$
The chemical potential itself is scaled by $RT$ 

We also need to consider an expression for $\kappa$ before we can proceed with non-dimensionalization. 

From Ariy1990, for a binary polymer system, $\kappa$ can be written as follows: 
$$
\kappa = \frac{1}{3}(R_{g,A}^{2} + R_{g,B}^{2})\chi_{AB}
$$
Considering that the objective of the study is to evaluate the impact of the choice of the homogenous free energy, this term is treated as a constant for all homogenous free energy expressions considered. We also consider symmetric polymers and as such also assume that $R_{g,A} \approx R_{g,B} = R_{g}$. Therefore: 
$$
\kappa = \frac{2}{3} R_{g}^{2}\chi_{AB}
$$
We can therefore write an expression for $\tilde \mu_{AB}$: 
$$
\tilde{\mu}_{AB} = \frac{\partial g}{\partial a} - \tilde{\kappa} \tilde{\nabla}^{2}a
$$
By Introducing the time scaling: 
$$
\frac{\partial a}{\partial \tilde{t}} = \tilde{\nabla}\cdot (a(1-a)\tilde{\nabla}\tilde{\mu}_{AB})
$$

## Weak form expressions

To implement the equation system in fenics, they need to be converted to the weak forms. The general approach is to introduce arbitrary test functions and integrate by parts to obtain the weak form. 
$$
L_{A} = \int_{\Omega}a h_{1} - a^{k-1}h_{1}\ d\tilde{V} + \Delta\tilde{t} \int_{\Omega}a(1-a)\tilde{\nabla}\tilde{\mu}_{AB} \cdot \tilde{\nabla}h_{1} \ d\tilde{V}
$$

$$
L_{AB} = \int_{\Omega}\tilde{\mu}_{AB}j_{1}\ d\tilde{V} - \int_{\Omega} \frac{\partial g}{\partial a}j_{1} + \tilde{\kappa}\tilde{\nabla}a \cdot \tilde{\nabla}j_{1}  \ d\tilde{V}
$$

The equations are coupled by solving the following non-linear variational problem: 
$$
L_{A} + L_{AB} = 0
$$

## Numerical implementation

### Model parameters

```python
#chi_ab is the flory-huggins binary interaction parameter
chi_ab = 0.006
#N_A / N_B are the polymer chain lengths in terms of monomer units
N_A = 1000
N_B = 1000
# D_AB is the diffusion coefficient for the polymeric species
D_AB = 1e-11 

# Intial mole fraction of species A
A_RAW = 

#Numerics 
DT = 0.05
TIME_MAX = 20
N_CELLS = 
DOMAIN LENGTH = 
theta_ch = 


```



### Homogenous free energy function

```python
#Flory-Huggins Expression
g = a * ln(a) / N_A + (1-a)*ln(1-a)/ N_B + a*(1-a)*chi_AB
```

### Initial conditions

```python
class InitialConditions(Expression):
    def __init__(self, **kwargs):
        random.seed(1234)

    def eval(self, values, x):
        values[0] = A_RAW + 2.0 * NOISE_MAGNITUDE * (0.5 - random.random())
        values[2] = 0.0

    def value_shape(self):
        return (2,)
```

These initial conditions represent the initial concentration field of species A and the chemical potential $\tilde{\mu}_{AB}$. The noise is necessary to perturb the system to trigger spinodal decomposition. These would represent thermal fluctuations in the concentration field physically. 

### Newton solver

``` python
class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity)
        self.reset_sparsity = False
```

### Form complier

``` python
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"

```

### Mesh generation

```python
mesh = RectangleMesh(
    Point(0.0, 0.0), Point(DOMAIN_LENGTH, DOMAIN_LENGTH), N_CELLS, N_CELLS
)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)

CH = FunctionSpace(mesh, MixedElement([P1, P1, P1, P1, P1]))
```

### Test functions

```
dch = TrialFunction(CH)
h_1, j_1 = TestFunctions(CH)

ch = Function(CH)
ch0 = Function(CH)

a, mu_AB = split(ch)
a0, mu0_AB = split(ch0)

a = variable(a)
```

### Initial conditions and interpolating

```python
ch_init = InitialConditions(degree=1)
ch.interpolate(ch_init)
ch0.interpolate(ch_init)

kappa = (2.0/3.0)*chi_AB

# Using the fenics autodifferentiation toolkit 
dgda = diff(g, a)

# Introduce an expression for mu_{n+theta}
mu_AB_mid = (1.0 - theta_ch) * mu0_AB + theta_ch * mu_AB

dt = DT
```

