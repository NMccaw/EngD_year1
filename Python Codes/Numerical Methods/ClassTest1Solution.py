"""
Class Test Template.

Please do not modify existing lines.
"""

# Average marks for MATH6141 students were 8.4. Penalty
# marks for scripts failing to run (1 mark), failing to output (by
# using a semicolon to suppress the answer 0.5 marks per question), or
# occasionally for excess output (0.5 marks per question) averaged
# 0.25 marks. More detailed feedback as below.
#

from __future__ import division

import numpy
import numpy.linalg
from matplotlib import pyplot
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
import scipy
import scipy.linalg
import scipy.optimize
import scipy.integrate

# If you want to add more import statements, please do it here.

# Answers to questions should follow, in the gaps.

print("\nSolution to question 1 follows.\n")

# Enter and invert a matrix output the inverse
#
# Average mark 1.
#
A = numpy.array([[2, 7], [3, 6]])
print(numpy.linalg.inv(A))

print("\nSolution to question 2 follows.\n")

# Create a linearly spaced vector, 20 points between 1 and 3. Output the 2
# norm of the vector
#
# Average mark 1.
v = numpy.linspace(1, 3, 20)
print(numpy.linalg.norm(v, 2))


print("\nSolution to question 3 follows.\n")
# Using the work in question 2, produce a vector with 20 points containing
# e^{-x} sin(2 pi x) between 1 and 3. Output the sum of the even elements.
#
# Usual error: typos for the exponent of e
#
# Average mark 0.83
y = numpy.exp(-v)*numpy.sin(2 * numpy.pi * v)
print(sum(y[1::2]))


print("\nSolution to question 4 follows.\n")
# Solve the linear system C x = b output x
#
# Average mark 1.

C = numpy.array([[12, 3, -2], [6, 4, 1], [0, 9, 2]])
b = numpy.array([1, 7, 3])
print(numpy.linalg.solve(C,b))


print("\nSolution to question 5 follows.\n")
# Find the eigenvalues of C output the eigenvalue with smallest modulus
#
# Usual error was outputting abs(lambda) instead of lambda.
#
# Average mark 0.63
e, v = numpy.linalg.eig(C)
index = numpy.argmin(numpy.abs(e))
print(e[index])


print("\nSolution to question 6 follows.\n")
# Produce an off-diagonal matrix using diag.
#
# Occasional error: diag of identity matrix, rather than of vector
#
# Average mark 0.92
print(numpy.diag(numpy.array([1, 1]), 1))


print("\nSolution to question 7 follows.\n")
# In Figure 1, plot e^{-x} sin(2 pi x) between 0 and 5. Add axis
# labels and a title.
#
# Usual error was not correctly labelling the axes
#
# Average mark 0.92
x = numpy.arange(0, 5, 0.01)
pyplot.plot(x, numpy.exp(-x) * numpy.sin(2 * numpy.pi * x))
pyplot.xlabel(r"$x$")
pyplot.ylabel(r"$e^{-x} sin(2 \pi x)$")
pyplot.title("Figure for question 7")
pyplot.show()

print("\nSolution to question 8 follows.\n")
# In figure 2, give a surface plot of e^{-(x^2+y^2)} cos(2 pi x^2) sin(4 pi y).
# Plot this on x \in [0, 1], y \in [0, 2] with spacing 0.05 in both
# directions. Use the default viewing angle but a hot colormap and no
# mesh lines
#
# Common errors: not correctly using rstride and cstride,
# leaving the meshgrid lines on
#
# Average mark 0.67
x = numpy.arange(0, 1, 0.05)
y = numpy.arange(0, 2, 0.05)
[X,Y] = numpy.meshgrid(x,y)
Z = numpy.exp(-(X**2 + Y**2)) * numpy.cos(2 * numpy.pi * X**2) * numpy.sin(4 * numpy.pi * Y)
fig = pyplot.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot, linewidth=0.)
pyplot.show()

print("\nSolution to question 9 follows.\n")

# Using inbuilt routines, compute the quadrature of
# e^{-x} sin(x) / (1 + x^2 + 0.1 x^4) between [0,1]. Output the result.
#
# Usual error: Error in the exponent of e again
#
# Average mark 0.79
def f(x):
    return numpy.exp(-x) * numpy.sin(x) / (1 + x**2 + 0.1 * x**4)
I, err = scipy.integrate.quad(f, 0, 1)
print(I)


print("\nSolution to question 10 follows.\n")
# Given that y_n = (a_n, b_n) and y_1 = (0.0, 0.5) iterate the map
#    y_{n+1} = (tan b_n, sin a_n)
# to find y_{10} and y_{100}. Output y_{10}, y_{100}
#
# Usual error: off-by-one errors in the indices
# (e.g. 10 and 100 rather than 9 and 99)
#
# Average mark 0.88
y = numpy.zeros((2, 100))
y[0, 0] = 0.0
y[1, 0] = 0.5
for i in range(1, 100):
    y[0, i] = numpy.tan(y[1, i-1])**3
    y[1, i] = numpy.cos(y[0, i-1])**2

print(y[:, 9])
print(y[:, 99])

print("\nEnd of solutions\n")
