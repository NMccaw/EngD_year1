import pytest

import numpy as np
import matplotlib
from matplotlib import pyplot
from scipy import integrate
from scipy import optimize

def question_1():
    """
    Solution to question 1 goes here
    """
    A=np.matrix([[1,3],[4,5]])
    return np.power(A,3)

def question_2():
    """
    Solution to question 2 goes here
    """
    v= np.linspace(1,2,40)
    return np.dot(v,v)

def question_3():
    """
    Solution to question 3 goes here
    """
    v= np.linspace(1,2,40)
    w=(v**2)*np.cos(np.pi*v)
    
    return np.sum(w[0::2])

def question_4():
    """
    Solution to question 4 goes here
    """
    A=np.array([[8,2,4],[2,-12,6],[4,6,1]])
    b=np.array([4,9,2])
    return np.linalg.solve(A,b)

def question_5():
    """
    Solution to question 5 goes here
    """
    A=np.array([[8,2,4],[2,-12,6],[4,6,1]])
    
    array=A-(2*np.identity(3))
    w,v= np.linalg.eig(array)
    return -w[1],-w[2]
    

def question_6():
    """
    Solution to question 6 goes here
    """
    A=np.array([[8,2,4],[2,-12,6],[4,6,1]])
    
    return np.tril(A,-1)

def question_7():
    """
    Solution to question 7 goes here
    """
    fig=pyplot.figure(1)
    x = numpy.arange(1, 4, 0.01)
    pyplot.plot(x, np.log(x)*np.sin(2*np.pi*x))
    pyplot.xlabel(r"$x$")
    pyplot.ylabel(r"$e^{-x} sin(2 \pi x)$")
    pyplot.title("Figure for question 7")
    pyplot.show()

def question_8():
    """
    Solution to question 8 goes here
    """
    x = np.arange(0, 1, 0.05)
    y = np.arange(0, 2, 0.05)
    [X,Y] = np.meshgrid(x,y)
    Z = np.exp(-(X**2)) * np.cos(2 * np.pi * Y**2)
    fig = pyplot.figure(2)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.hot, linewidth=0.)
    pyplot.show()

def question_9():
    """
    Solution to question 9 goes here
    """
    fun= lambda s:np.cos(s)**3-s**2-s
    s0=1
    return scipy.optimize.newton(fun,s0)


def question_10():
    """
    Solution to question 10 goes here
    """
    y = np.zeros((2, 100))
    y[0, 0] = 0.25
    y[1, 0] = 0.75
    for i in range(1, 100):
        y[0, i] = np.cos(y[1, i-1])
        y[1, i] = np.sin(y[0, i-1])

    return (y[:, 4]),(y[:, 49])


if __name__ == "__main__":
    pytest.main()
    question_7()
    question_8()