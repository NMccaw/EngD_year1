def mandel_cy(long n, long itermax=100, long xmin=-2,double xmax=0.5,double ymin=-1.25,double ymax=1.25):
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
    cdef:
        long ix,iy,it, itermax_copy = itermax
        double x,y,z_real,z_imag

    # iterate through all matrix elements
    for ix in range(0, n):
        for iy in range(0, n):
            # compute the position (x, y) corresponding to matrix element
            x = xmin + ix * (xmax - xmin) / (n)
            y = ymin + iy * (ymax - ymin) / (n)
            # Need to count iterations
            it = 0
            # c is the complex number with the given
            # x, y coordinates in the complex plane, i.e. c = x + i * y
            # where i = sqrt(-1)


            z_real = 0
            z_imag = 0


            # Here is the actual Mandelbrot criterion: we update z to be
            # z <- z^2 + c until |z| <= 2. We could the number of iterations
            # required. This number of iterations is the data we need to compute
            # (and plot if desired).

            while it < itermax_copy and (z_real*z_real)+(z_imag*z_imag) < 4.0:




                z_square_real = ((z_real**2)-(z_imag**2))
                z_square_imag = (z_imag*z_real*2)
                z_real = x + z_square_real
                z_imag = y + z_square_imag


                it = it+1

            #print("ix={}, iy={}, x={}, y={}, c={}, z={}, abs(z)={}, it={}"
            #    .format(ix, iy, x, y, c, z, abs(z), it))

            # Store the result in the matrix
            its[ix][iy] = it

    return its
