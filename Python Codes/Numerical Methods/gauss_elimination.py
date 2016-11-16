# -*- coding: utf-8 -*-
"""
gaussian_elimination
====================

This file is a combination of templates, functions and tests for
implementing a Gaussian Elimination method.
"""

import numpy
import numpy.linalg

def MyBackSubstitution(A, b):
    """
    Solve the upper triangular linear system A x = b.

    Parameters
    ----------

    A : array of float
        real square matrix
    b : vector of float
        real vector

    Returns
    -------

    x : vector of float
        solution

    Notes
    -----

    Simplified method with limited error checking.
    """

    assert(numpy.all(numpy.isreal(b))), "b must be real"
    assert(numpy.all(numpy.isfinite(b))), "b must be finite"
    assert(numpy.ndim(b) == 1), "b must be a vector"
    n = len(b)

    assert(numpy.all(numpy.isreal(A))), "A must be real"
    assert(numpy.all(numpy.isfinite(A))), "A must be finite"
    assert(numpy.ndim(A) == 2), "A must be a matrix"
    assert(A.shape == (n, n)), "A must be a square matrix compatible with b"

    x = numpy.zeros_like(b)

    for i in range(n-1,-1,-1):
        x[i] = b[i] / A[i, i]
        for k in range(i+1,n):
            x[i] -=  A[i, k] * x[k] / A[i, i]

    return x

def MyGaussianElimination(A, b):
    """
    Solve the linear system A x = b using Gaussian Elimination without pivoting.

    Parameters
    ----------

    A : array of float
        real square matrix
    b : vector of float
        real vector

    Returns
    -------

    x : vector of float
        solution

    Notes
    -----

    Simplified method with limited error checking.
    """

    # Error checking here
    assert(numpy.all(numpy.isreal(b))), "b must be real"
    assert(numpy.all(numpy.isfinite(b))), "b must be finite"
    assert(numpy.ndim(b) == 1), "b must be a vector"
    n = len(b)

    assert(numpy.all(numpy.isreal(A))), "A must be real"
    assert(numpy.all(numpy.isfinite(A))), "A must be finite"
    assert(numpy.ndim(A) == 2), "A must be a matrix"
    assert(A.shape == (n, n)), "A must be a square matrix compatible with b"

    # Construct augmented matrix. Slightly tedious.
    aug = numpy.hstack((A, numpy.reshape(b, [len(b), 1])))

    # Put the augmented matrix in triangular form.
    N=len(b)

    for j in range(0,N):
        assert(numpy.abs(aug[j,j]) > 1e-20), "Diagonal element zero!"
        for i in range(j+1,N):
            aug[i,:]=aug[i,:]-((aug[i,j]/aug[j,j]*aug[j,:]))




    # Solve using back substitution.
    x = MyBackSubstitution(aug[:, :-1], aug[:, -1])

    return x


def MyGaussianEliminationWithPivoting(A, b):
    """
    Solve the linear system A x = b using Gaussian Elimination with pivoting.

    Parameters
    ----------

    A : array of float
        real square matrix
    b : vector of float
        real vector

    Returns
    -------

    x : vector of float
        solution

    Notes
    -----

    Simplified method with limited error checking.
    """

    # Error checking here

    # Construct augmented matrix. Slightly tedious.
    aug = numpy.hstack((A, numpy.reshape(b, [len(b), 1])))

    # Put the augmented matrix in triangular form.
    raise NotImplementedError("Code needed here")

    # Solve using back substitution.
    x = MyBackSubstitution(aug[:, :-1], aug[:, -1])

    return x

# What follows are testing functions to validate the code
import pytest

def test_diagonal():
    A = numpy.eye(2)
    b = numpy.array([1.0, 2.0])
    x_my = MyGaussianElimination(A, b)
    check = numpy.allclose(x_my, b)
    assert check

def test_triangular():
    A = numpy.array([[1.0, 2.0], [0.0, 1.0]])
    b = numpy.array([4.0, 1.0])
    x_my = MyGaussianElimination(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_full():
    A = numpy.array([[1.0, 2.0], [3.0, 4.0]])
    b = numpy.array([5.0, 6.0])
    x_my = MyGaussianElimination(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_threebythree():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = numpy.array([4.0, 10.0, 15.0])
    x_my = MyGaussianElimination(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_incompatible():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = numpy.array([4.0, 10.0])
    with pytest.raises(AssertionError):
        MyGaussianElimination(A, b)

def test_input():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = "dog"
    with pytest.raises(AssertionError):
        MyGaussianElimination(A, b)

def test_singular():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 5.0]])
    b = numpy.array([4.0, 10.0])
    with pytest.raises(AssertionError):
        MyGaussianElimination(A, b)

def test_finite():
    A = numpy.array([[1.0, 1.0, 1.0], [0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    b = numpy.array([1.0, 1.0, 2.0])
    with pytest.raises(AssertionError):
        MyGaussianElimination(A, b)

def test_needs_pivoting():
    A = numpy.array([[1.0e-20, 1.0], [1.0, 1.0]])
    b = numpy.array([1.0, 2.0])
    with pytest.raises(AssertionError):
        MyGaussianElimination(A, b)

# Test with pivoting

def test_diagonal_pivoting():
    A = numpy.eye(2)
    b = numpy.array([1.0, 2.0])
    x_my = MyGaussianEliminationWithPivoting(A, b)
    check = numpy.allclose(x_my, b)
    assert check

def test_triangular_pivoting():
    A = numpy.array([[1.0, 2.0], [0.0, 1.0]])
    b = numpy.array([4.0, 1.0])
    x_my = MyGaussianEliminationWithPivoting(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_full_pivoting():
    A = numpy.array([[1.0, 2.0], [3.0, 4.0]])
    b = numpy.array([5.0, 6.0])
    x_my = MyGaussianEliminationWithPivoting(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_threebythree_pivoting():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = numpy.array([4.0, 10.0, 15.0])
    x_my = MyGaussianEliminationWithPivoting(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

def test_incompatible_pivoting():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = numpy.array([4.0, 10.0])
    with pytest.raises(AssertionError):
        MyGaussianEliminationWithPivoting(A, b)

def test_input_pivoting():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 6.0]])
    b = "dog"
    with pytest.raises(AssertionError):
        MyGaussianEliminationWithPivoting(A, b)

def test_singular_pivoting():
    A = numpy.array([[3.0, 0.0, 1.0], [6.0, 2.0, 4.0], [9.0, 2.0, 5.0]])
    b = numpy.array([4.0, 10.0])
    with pytest.raises(AssertionError):
        MyGaussianEliminationWithPivoting(A, b)

def test_finite_pivoting():
    A = numpy.array([[1.0, 1.0, 1.0], [0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    b = numpy.array([1.0, 1.0, 2.0])
    with pytest.raises(AssertionError):
        MyGaussianEliminationWithPivoting(A, b)

def test_needs_pivoting_pivoting():
    A = numpy.array([[1.0e-20, 1.0], [1.0, 1.0]])
    b = numpy.array([1.0, 2.0])
    x_my = MyGaussianEliminationWithPivoting(A, b)
    x_exact = numpy.linalg.solve(A, b)
    check = numpy.allclose(x_my, x_exact)
    assert check

# Run all the tests
pytest.main("-x gauss_elimination.py")