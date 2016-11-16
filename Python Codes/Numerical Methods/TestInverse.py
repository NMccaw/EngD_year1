# This script requires an intro comment
def TestInverse(A, tol):
    """
    Returns the inverse of matrix A. Calulates the condition number and the determinant.
    Determines if the condition number and determinant is too big or small respectively
    given a tolerance
    """
    import numpy as np
    import numpy.linalg

    invtol = 1.0 / tol

    assert(np.linalg.cond((A) < tol)), 'The condition number is too big!'
    assert(np.linalg.det(A < invtol)), 'The determinant of the matrix is too small!'

    InverseA = np.linalg.inv(A)

    return InverseA
