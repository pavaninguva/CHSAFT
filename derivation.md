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


