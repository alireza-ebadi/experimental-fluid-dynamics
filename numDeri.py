import numpy as np

### calculate the first derivative using 3-term Taylor approximation
def numericalDerivative(x, y):
    '''
	calculate the first derivative using 3-term Taylor approximation

	inputs:
    x, y (as arrays)

	outputs:
    dy/dx
	'''
	
	m = len(y)
	dydx = np.zeros(m)

	ii = 0
	dx1 = x[ii + 1] - x[ii]
	dx2 = x[ii + 2] - x[ii]
	f = -1/(dx2 - dx1)*dx1/dx2
	e = -f*dx2**2/dx1**2
	d = -e - f
	dydx[ii] = d*y[ii] + e*y[ii + 1] + f*y[ii + 2]

	for ii in range(1, m - 1):
	    dx1 = x[ii - 1] - x[ii]
	    dx2 = x[ii + 1] - x[ii]
	    f = -1/(dx2 - dx1)*dx1/dx2
	    e = -f*dx2**2/dx1**2
	    d = -e - f
	    dydx[ii] = d*y[ii] + e*y[ii - 1] + f*y[ii + 1]

	ii = m - 1
	dx1 = x[ii - 1] - x[ii]
	dx2 = x[ii - 2] - x[ii]
	f = -1/(dx2 - dx1)*dx1/dx2
	e = -f*dx2**2/dx1**2
	d = -e - f
	dydx[ii] = d*y[ii] + e*y[ii - 1] + f*y[ii - 2]

	return dydx
