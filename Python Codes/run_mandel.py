doc = """
Little demo called run_mandel.py, helper for lab 9, FEEG6002.

Hans Fangohr, fangohr@soton.ac.uk, November 2013
Updated for Python 3 by Neil O'Brien, nsob1c12@soton.ac.uk, November 2016

---------------------

Call this script with

    --test [n]

       to compare Python and Cython results for a n x n matrix,
       n defaults to 4

    --demo [n]

       co compute a n x n matrix, and save to file mandel.png and
       display to screen.
       Will use cython version, if it exists.
       n defaults to 400

    --speedup [n]

       Compare Python and Cython speed for n x n matrix.
       Default n is 400.

    --help

       show this text.



Example:

  python run_mandel.py --demo 200

     computes 200 x 200 matrix, saves to file, and displays on screen.

"""


import timeit
import sys

import numpy
import pylab

import pyximport
pyximport.install()

# import module mandel.pyx with cython code
import mandel


def mandel_py(n, itermax=100, xmin=-2, xmax=0.5, ymin=-1.25, ymax=1.25):
    '''
    Mandelbrot fractal computation using Python.

    (n, n) are the output image dimensions
    itermax is the maximum number of iterations to do.
    xmin, xmax, ymin, ymax specify the region of the
    set to compute.

    Returns a list of lists of ints, representing the image matrix with
    n rows and n columns.
    '''

    # create list containing n lists, each containing n zeros
    # (i.e. a matrix, represented as a list of lists)
    its = [ [0] * n for i in range(n)]
    # The data in the matrix are iterations, so 'its' is the plural of
    # IT for ITeration.


    # iterate through all matrix elements
    for ix in range(0, n):
        for iy in range(0, n):
            # compute the position (x, y) corresponding to matrix element
            x = xmin + ix * (xmax - xmin) / float(n)
            y = ymin + iy * (ymax - ymin) / float(n)
            # Need to count iterations
            it = 0
            # c is the complex number with the given
            # x, y coordinates in the complex plane, i.e. c = x + i * y
            # where i = sqrt(-1)
            c = x + y * 1j
            z = 0
            # Here is the actual Mandelbrot criterion: we update z to be
            # z <- z^2 + c until |z| <= 2. We could the number of iterations
            # required. This number of iterations is the data we need to compute
            # (and plot if desired).
            while it < itermax and abs(z) < 2.0:
                z = z ** 2 + c
                it += 1

            #print("ix={}, iy={}, x={}, y={}, c={}, z={}, abs(z)={}, it={}"
            #    .format(ix, iy, x, y, c, z, abs(z), it))

            # Store the result in the matrix
            its[ix][iy] = it

    return its


def test_equality(n):
    I_cy = numpy.array(mandel.mandel_cy(n, 100, -2, .5, -1.25, 1.25))
    I_py = numpy.array(mandel_py(n, 100, -2, .5, -1.25, 1.25))


    diff = I_cy - I_py
    if max(numpy.abs(diff.flat)) > 0:
        print("Something bifferent. Difference is")
        print(diff)
        print("Cython matrix:")
        print(I_cy)
        print("Python matrix:")
        print(I_py)
    else:
        print("Python and Cython results agree for n={}".format(n))

def demo_save_png(n):
    print("Computing {} x {} matrix".format(n, n))

    # try to use cython version
    if hasattr(mandel, "mandel_cy"):
        print("Using Cython version to compute")
        I = mandel.mandel_cy(n, 100, -2, .5, -1.25, 1.25)
    else:
        print("Using Python version to compute")
        I = mandel_py(n, 100, -2, .5, -1.25, 1.25)
    I = numpy.array(I)
    img = pylab.imshow(I.T, origin='lower left')
    print("Writing data to 'mandel.png'")
    img.write_png('mandel.png', noscale=True)
    pylab.show()


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need some command line switch. Use '--help' for help.")
        print("I.e., run 'python run_mandel.py --help'")
        sys.exit(0)
    else:
        mode = sys.argv[1]
        if len(sys.argv) == 3:
            N = int(sys.argv[2])
        else:
            N = None

        # now decide what to do

        # should we run demo, i.e. compute and save as bitmap?
        if mode == '--demo':
            if N is None:
                N = 400
            demo_save_png(N)
        elif mode == '--test':
            if N is None:
                N = 4
            test_equality(N)
        elif mode == '--help':
            print(doc)
        elif mode == '--speedup':
            if N is None:
                N = 400

            print ("Doing timings with {} x {} pixels".format(N, N))
            print ("Python code ...")
            t = timeit.Timer("mandel_py(N)", setup="from run_mandel import mandel_py; N={}".format(N))
            time_python = min(t.repeat(number=1, repeat=5))
            print("Time python: {} s".format(time_python))

            print("Cython code ...")
            t = timeit.Timer("mandel_cy(N)", setup="from mandel import mandel_cy; N={}".format(N))
            time_cython = min(t.repeat(number=1, repeat=5))
            print("Time cython: %f s" % time_cython)
            print("Speed up: %.1f" % (time_python / time_cython))
        else:
            print("Didn't understand '{}'".format(mode))
            print("Run 'python run_mandel.py --help' for list of available options")
        sys.exit(0)
