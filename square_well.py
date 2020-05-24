# Infinite square well

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
a = 30.0
xvals = np.linspace(0, 30, num = 1000)

A = np.sqrt(2 / a)   # Normalization constant

def psi(n, x):
    wf = A * np.sin((n * pi * x) / a)  # Defining wave function
    return wf

yvals1 = psi(1, xvals)
yvals2 = psi(2, xvals)
yvals3 = psi(3, xvals)

plt.plot(xvals, yvals1, 'b--', label="n=1")
plt.plot(xvals, yvals2, 'ro', label="n=2")
plt.plot(xvals, yvals3, 'k.', label="n=3")
plt.title("Wave function for Infinity Square Well")
plt.xlabel("x values")
plt.ylabel("$\psi$(x)")
plt.grid()

plt.savefig('plot.png')
plt.show()

