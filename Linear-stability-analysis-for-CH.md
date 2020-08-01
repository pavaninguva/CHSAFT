# Linear stability analysis for CH

## Vanilla CH with quartic potential and simplified mobility

These equations take the form: 
$$
\frac{\partial c}{\partial t} - \nabla \cdot M \nabla \mu = 0
$$

$$
\mu - \frac{\partial g}{\partial c} + \kappa\nabla^{2}c = 0
$$

where: 
$$
M = 1
$$

$$
g = Ac^{2}(1-c)^{2}
$$

$$
k = constant
$$

$A$ is a constant that scales $g$ 

When we expand the equations out fully: 
$$
\frac{\partial c}{\partial t} - \nabla \cdot M \nabla \bigg(\frac{\partial g}{\partial c} - \kappa\nabla^{2}c \bigg) = 0
$$
Applying the distributive property of the gradient operator:
$$
\frac{\partial c}{\partial t} - \nabla \cdot M \bigg( \nabla\frac{\partial g}{\partial c} - \nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
Assuming Mobility is not a constant for now and bringing the mobility in:
$$
\frac{\partial c}{\partial t} - \nabla \cdot \bigg( M\nabla\frac{\partial g}{\partial c} - M\nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
We apply the distributive property of the divergence operator: 
$$
\frac{\partial c}{\partial t} -\bigg( \nabla \cdot M\nabla \frac{\partial g}{\partial c}  \bigg) + \bigg(\nabla \cdot M\nabla(\kappa\nabla^{2}c)\bigg) = 0
$$
Since mobility is not a constant and can be treated as a scalar, we apply the product rule:
$$
\frac{\partial c}{\partial t} -\bigg( M\nabla\cdot \nabla \frac{\partial g}{\partial c} + (\nabla M) \cdot \nabla \frac{\partial g}{\partial c} \bigg) + \bigg(M\nabla \cdot \nabla(\kappa\nabla^{2}c) + (\nabla M) \cdot \nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
Simplifying the above, noting that $\kappa$ is a constant:
$$
\frac{\partial c}{\partial t} -\bigg( M\nabla^{2} \frac{\partial g}{\partial c} + (\nabla M) \cdot \nabla \frac{\partial g}{\partial c} \bigg) + \bigg(M \kappa \nabla^{4} c + (\nabla M) \cdot \nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
We can take the above equation as the general 2 component CH equation for now. 

We introduce the various simplifications made available by the vanilla CH type problems: 
$$
M = Constant 
$$

$$
\frac{\partial g}{\partial c} = A(2c(1-c)^{2} + 2(1-c)(-1)c^{2})
$$

$$
\frac{\partial g}{\partial c} = A(2c(1-2c +c^{2}) - 2(1-c)c^{2})
$$

$$
\frac{\partial g}{\partial c} = A((2c-4c^{2} +2c^{3}) - (2c^{2}-2c^{3}))
$$

$$
\frac{\partial g}{\partial c} = A(2c-6c^{2} +4c^{3})
$$

When $M$ is a constant, all the $\nabla M$ terms are 0
$$
\frac{\partial c}{\partial t} -\bigg( M\nabla^{2} \frac{\partial g}{\partial c}  \bigg) + \bigg(M \kappa \nabla^{4} c  \bigg) = 0
$$
Introducing the expression for the free energy: 
$$
\frac{\partial c}{\partial t} -\bigg( MA\nabla^{2} (2c-6c^{2} +4c^{3})  \bigg) + \bigg(M \kappa \nabla^{4} c  \bigg) = 0
$$

$$
\frac{\partial c}{\partial t} -\bigg( 2MA\nabla^{2}c  - 6MA \nabla^{2} c^{2} + 4MA\nabla^{2} c^{3}  \bigg) + \bigg(M \kappa \nabla^{4} c  \bigg) = 0
$$

To cope with the terms with $c^{2}$ and higher, we employ some vector calculus identities: 

For $\nabla^{2} c^{2}$: 
$$
\nabla ^{2}c^{2} = 2c\nabla^{2}c + 2(\nabla c \cdot \nabla c) = 2c \nabla^{2}c + 2||\nabla c||^{2}
$$


For $\nabla^{2} c^{3}$: 
$$
\nabla^{2}c^{3} = \nabla^{2}(c c^{2}) = c\nabla^{2} c^{2} + c^{2}\nabla^{2} c + 2(\nabla c^{2} \cdot \nabla c )
$$

$$
\nabla^{2}c^{3}  = c\nabla^{2} c^{2} + c^{2}\nabla^{2} c + 4(c\nabla c \cdot \nabla c )
$$

$$
\nabla^{2}c^{3}  = c\nabla^{2} c^{2} + c^{2}\nabla^{2} c + 4c||\nabla c ||^{2}
$$

$$
\nabla^{2}c^{3}  = c (2c\nabla^{2}c + 2 ||\nabla c ||^{2}) + c^{2}\nabla c + 4c||\nabla c ||^{2}
$$

$$
\nabla^{2}c^{3}  = 2c^{2}\nabla^{2}c   + c^{2}\nabla c + 6c||\nabla c ||^{2}
$$

Substituting back into the original equation: 
$$
\frac{\partial c}{\partial t} \\ -   2MA\nabla^{2}c  \\ - 6MA \nabla^{2} c^{2} \\ + 4MA\nabla^{2} c^{3}  \\ + M \kappa \nabla^{4} c = 0
$$



$$
\frac{\partial c}{\partial t} \\ -   2MA\nabla^{2}c  \\ - 6MA (2c \nabla^{2}c + 2||\nabla c||^{2}) \\ + 4MA (2c^{2}\nabla^{2}c   + c^{2}\nabla c + 6c||\nabla c ||^{2} )  \\ + M \kappa \nabla^{4} c \\ = 0
$$
This is the final equation on which we shall perform the linear stability analysis for! We keep $M$, $\kappa$ and $A$ as free constants to understand the impact of these terms. 

We introduce a perturbation: 
$$
c(\textbf{x}, t) = c_{0} + \delta c(\textbf{x}, t)
$$
where: 
$$
\delta c(\textbf{x},t) = \epsilon c(t) e^{-ik \textbf{x}}
$$
where: 
$$
c(t) = e^{\sigma t}
$$
This gives us the ansatz form: 
$$
c(\textbf{x}, t) = c_{0} + \epsilon e^{-ik\textbf{x} + \sigma t}
$$
At this point, we can evaluate the derivatives for substitution: 
$$
\frac{\partial c}{\partial t} = \sigma \epsilon e^{-ik\textbf{x} + \sigma t}
$$

$$
\nabla c = -ik\epsilon e^{-ik\textbf{x} + \sigma t}
$$

$$
\nabla ^{2} c = ik\epsilon e^{-ik\textbf{x} + \sigma t} (-ik) = k^{2}\epsilon e^{-ik\textbf{x} + \sigma t}
$$

$$
\nabla ^{3} c = k^{2}\epsilon e^{-ik\textbf{x} + \sigma t} (-ik) = -ik^{3}\epsilon e^{-ik\textbf{x} + \sigma t}
$$

$$
\nabla^{4} c = -ik^{3}\epsilon e^{-ik\textbf{x} + \sigma t} (-ik) = k^{4}\epsilon e^{-ik\textbf{x} + \sigma t}
$$

Substituting these expressions into the CH equation and discarding terms of $O (\epsilon^{2})$ or higher : 
$$
\sigma \epsilon e^{-ik\textbf{x} + \sigma t}  -2MA(k^{2}\epsilon e^{-ik\textbf{x} + \sigma t}) + M\kappa k^{4}\epsilon e^{-ik\textbf{x} + \sigma t} = 0
$$
Further simplifying and since $M$ is unity:
$$
\sigma = 2A k^{2} - \kappa k^{4}
$$
Interesting that we get two competing effects here where an increasing depth of the free energy potential $A$ results in a more unstable system while larger $k$ values contribute to stability which makes complete sense here cause the creation of an interface penalises demixing, but the free energy of mixing drives demixing. 

### checking consistency of approaches

For sanity's sake, lets test performing a Taylor expansion for the gibbs free energy instead of rejigging the higher order terms. 
$$
\frac{\partial g}{\partial c} = A(2c-6c^{2} +4c^{3})
$$

$$
\frac{\partial^{2} g}{\partial c ^{2}} = A(2-12c + 12c^{2})
$$

Performing the linearization: 
$$
\frac{\partial g}{\partial c} \approx \frac{\partial g}{\partial c} \biggr\vert_{c_{0}} + \frac{\partial ^{2} g }{\partial c^{2}} \biggr \vert_{c_{0}} (c-c_{0})
$$
We end up with the following: 
$$
\frac{\partial g}{\partial c} \approx A(2c_{0} - 6c_{0}^{2} + 4c_{0}^{3}) + A(2 - 12c_{0} + 12 c_{0}^{2})(c-c_{0})
$$

$$
\frac{\partial g}{\partial c} \approx A(6c_{0}^{2 } - 8 c_{0}^{3}) + A(2 - 12c_{0} + 12c_{0}^{2}) c
$$

Throwing these back into the transport equations: 
$$
\frac{\partial c}{\partial t} -\bigg( \nabla^{2} (A(6c_{0}^{2 } - 8 c_{0}^{3}) + A(2 - 12c_{0} + 12c_{0}^{2}) c)  \bigg) + \bigg( \kappa \nabla^{4} c  \bigg) = 0
$$
This simplifies to: 
$$
\frac{\partial c}{\partial t} -\bigg( A(2 - 12c_{0} + 12c_{0}^{2})\nabla^{2} c  \bigg) + \bigg( \kappa \nabla^{4} c  \bigg) = 0
$$


Which upon introducing the ansatz expression: 
$$
\sigma \epsilon e^{-ik\textbf{x} + \sigma t} -A(2 - 12c_{0} + 12c_{0}^{2}) k^{2}\epsilon e^{-ik\textbf{x} + \sigma t} + \kappa k^{4}\epsilon e^{-ik\textbf{x} + \sigma t} = 0
$$
This reduces down to: 
$$
\sigma = A(2 - 12c_{0} + 12c_{0}^{2}) k^{2} - \kappa k^{4}
$$
Which is practically similar to the other way of doing it as outlined above. 

## Simple mobility, but with F-H free energy expression

Now that we got a working framework for the linear stability analysis, we can explore how things are behaving for different CH models. The first one adopts a simple mobility expression, but replaces the quartic polynomial with the Flory-Huggins expression. 

The transport equation: 

Assumptions: 

- Constant $\kappa$ 
- Constant mobility $M$
- $\chi$ is not temporally / spatially variant. 

$$
\frac{\partial c}{\partial t} -\bigg( M\nabla^{2} \frac{\partial g}{\partial c}  \bigg) + \bigg(M \kappa \nabla^{4} c  \bigg) = 0
$$

Lets evaluate the first and second derivatives $g$: 
$$
g(c) = \frac{c}{N_{1}} \ln{c} + \frac{1-c}{N_{2}} \ln(1-c) + \chi c(1-c)
$$

$$
\frac{\partial g}{\partial c} = \bigg(\frac{1}{N_{1}} - \frac{1}{N_{2}} \bigg) + \frac{1}{N_{1}}\ln{c} - \frac{1}{N_2}\ln{(1-c)} + \chi(1-2c)
$$

$$
\frac{\partial g}{\partial c} = \bigg(\frac{1}{N_{1}} - \frac{1}{N_{2}} + \chi \bigg) + \frac{1}{N_{1}}\ln{c} - \frac{1}{N_2}\ln{(1-c)}  -2c \chi
$$

$$
\frac{\partial ^{2} g }{\partial c^{2}} = \frac{1}{N_{1}c} + \frac{1}{N_{2}(1-c)} - 2\chi
$$

We need to carry out a linearization of the free energy function. We perform a Taylor expansion about $c_{0}$ and keep only the first order terms. 
$$
\frac{\partial g}{\partial c} \approx \frac{\partial g}{\partial c} \biggr\vert_{c_{0}} + \frac{\partial ^{2} g }{\partial c^{2}} \biggr \vert_{c_{0}} (c-c_{0})
$$

$$
\frac{\partial g}{\partial c} \approx \bigg(\frac{1}{N_{1}} - \frac{1}{N_{2}} + \chi  + \frac{1}{N_{1}}\ln{c_{0}} - \frac{1}{N_B}\ln{(1-c_{0})} + -2c_{0} \chi \bigg) + \bigg( \frac{1}{N_{1}c_{0}} + \frac{1}{N_{2}(1-c_{0})} - 2\chi \bigg)(c-c_{0})
$$



Throwing this expression back into the transport equation: 
$$
\frac{\partial c}{\partial t} \\ - M\nabla^{2} \bigg( \bigg(\frac{1}{N_{1}} - \frac{1}{N_{2}} + \chi  + \frac{1}{N_{1}}\ln{c_{0}} - \frac{1}{N_B}\ln{(1-c_{0})} + -2c_{0} \chi \bigg) + \bigg( \frac{1}{N_{1}c_{0}} + \frac{1}{N_{2}(1-c_{0})} - 2\chi \bigg)(c-c_{0}) \bigg)  \\ + \bigg(M \kappa \nabla^{4} c  \bigg) = 0
$$
The first term of the free energy function is a constant and is annihilated and expanding out the $(c-c_{0})$ conveniently removes that term. We also acknowledge that $M$ is unity at this point: 
$$
\frac{\partial c}{\partial t} - \nabla^{2}\bigg( \frac{1}{N_{1}c_{0}} + \frac{1}{N_{2}(1-c_{0})} - 2\chi \bigg)c + \bigg( \kappa \nabla^{4} c  \bigg) = 0
$$
We introduce the ansatz form: 
$$
\sigma \epsilon e^{-ik\textbf{x} + \sigma t} - \bigg( \frac{1}{N_{1}c_{0}} + \frac{1}{N_{2}(1-c_{0})} - 2\chi \bigg) k^{2}\epsilon e^{-ik\textbf{x} + \sigma t} + \kappa k^{4}\epsilon e^{-ik\textbf{x} + \sigma t} = 0
$$
Simplifying: 
$$
\sigma = \bigg( \frac{1}{N_{1}c_{0}} + \frac{1}{N_{2}(1-c_{0})} - 2\chi \bigg) k^{2}  - \kappa k^{4}
$$


This reveals really interesting dynamics where increasing $\chi$ actually increases stability. 

## Flory Huggins free energy and real mobility

Now is the time where we introduce the real mobility. 
$$
\frac{\partial c}{\partial t} -\bigg( M\nabla^{2} \frac{\partial g}{\partial c} + (\nabla M) \cdot \nabla \frac{\partial g}{\partial c} \bigg) + \bigg(M \kappa \nabla^{4} c + (\nabla M) \cdot \nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
We need to perform a bunch of linearizations: 

We first write: 
$$
M = c(1-c) = c-c^{2}
$$
Lets throw this into the above equation: 
$$
\frac{\partial c}{\partial t} -\bigg( (c-c^{2})\nabla^{2} \frac{\partial g}{\partial c} + (\nabla (c-c^{2})) \cdot \nabla \frac{\partial g}{\partial c} \bigg) + \bigg((c-c^{2}) \kappa \nabla^{4} c + (\nabla M) \cdot \nabla(\kappa\nabla^{2}c) \bigg) = 0
$$
