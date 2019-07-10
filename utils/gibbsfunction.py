# from dolfin import *
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Define x_a as a variable
# x_a = variable(x_a)
chi_AB= 0.003
x_a = np.arange(0.0, 1.0, 0.001)
N_A = N_B = 1000

#Flory-Huggins Expression
g = ( x_a * np.log(x_a) / N_A ) + ((1.0-x_a)*np.log(1.0-x_a)/ N_B) + x_a*(1.0-x_a)*chi_AB


print(g)

# Using the fenics autodifferentiation toolkit 
# dgdx_a = diff(g,x_a)

fig, ax = plt.subplots()
line1, = ax.plot (x_a, g, label ='Gibbs energy of mixing')
# line2, = ax.plot(x, dgdx_a, label='Derivative')
ax.set(xlabel='mole fraction of A', ylabel='Energy')
fig.savefig('gibbs.png')
plt.show()