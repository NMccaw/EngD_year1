def TestSequence(N):
    """
    Compute the infinite sum of 2^{-n} starting from n = 0, truncating
    at n = N, returning the value of 2^{-n} and the truncated sum.

    Parameters
    ----------

    N : int
        Positive integer, giving the number of terms in the sum

    Returns
    -------

    limit : float
        The value of 2^{-N}
    sumseq : float
        The value of the truncated sum

    Notes
    -----

    The limiting value should be zero, and the value of the sum should
    converge to 2.
    """

    # Start sum from zero, so give zeroth term
    limit = 1.0
    sumseq = 1.0

    # At each step, increment sum and change summand
    for i in range(1, N+1):
        sumseq += limit
        limit /= 2.0

    return limit, sumseq

print(TestSequence(50))