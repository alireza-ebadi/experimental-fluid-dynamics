import numpy as np
from numDeri import numericalDerivative as deri

def qw(y, T, vt, Pr, nu, rhocp):

	'''
    calculate wall heat flux from dimensional pramaters using the integral method

    Reference:
    A. Ebadi, F. Mehdi and C.M. White,
    An exact integral method to evaluate wall heat flux in spatially developing two-dimensional wall-bounded flows,
    International Journal of Heat and Mass Transfer, 84 (2015): 856-861

    inputs:
    y (m): wall-normal position
    T (K): mean temperature
    vt (Km/s): = wall-normal Reynolds heat flux
    Pr (-): Prandtl number
    rhocp (kJ/m3.K): density*specific heat

    outputs:
    qw: wall heat flux
    terms: contributing terms
    '''

	n = len(y)
	y =  np.reshape(y, n)
	T =  np.reshape(T, n)
	vt =  np.reshape(vt, n)
	yt = y[-1]

	alpha = nu/Pr
	Tw = T[0]
	TwT = Tw - T
	temp = -deri(y, TwT)
	dTdy = np.reshape(temp, n)

	q = alpha*dTdy - vt
	temp = deri(y, q)
	dqdy = np.reshape(temp, n)

	int = np.zeros((3))
	int[0] = 2*alpha*np.trapz(TwT, y)
	int[1] = 2*np.trapz(vt*(yt - y), y)
	int[2] = np.trapz(dqdy*(yt - y)**2, y)

	terms = int/yt**2
	qw_over_rhocp = np.sum(terms)
	qw = qw_over_rhocp*rhocp

	return qw, terms



def tauw(y, U, uv, nu, rho):

	'''
    calculate wall shear stress from dimensional pramaters using the integral method

    Reference:
    Mehdi, F., Johansson, T.G., White, C.M. and Naughton, J.W.,
    On determining wall shear stress in spatially developing two-dimensional wall-bounded flows,
    Experiments in fluids, 55(1), p.1656 (2014).

	inputs:
    y (m): wall-normal position
    U (m/s): mean streamwise velocity
    uv (m2/s2): wall-normal Reynolds shear stress
    nu (m2/s): kinematic viscosity
    ho (kg/m3): density

    outputs:
    tauw: wall shear stress
    terms: contributing terms
	'''

	n = len(y)
	y =  np.reshape(y, n)
	U =  np.reshape(U, n)
	uv =  np.reshape(uv, n)
	yt = y[-1]

	temp = deri(y, U)
	dUdy = np.reshape(temp, n)

	tau = nu*dUdy - uv;
	temp = deri(y, tau);
	dtaudy = np.reshape(temp, n)

	int = np.zeros((3))
	int[0] = 2*nu*np.trapz(U, y)
	int[1] = 2*np.trapz(-uv*(yt - y), y)
	int[2] = -np.trapz(dtaudy*(yt - y)**2, y)

	terms = int/yt**2
	tauw_over_rho = np.sum(terms)
	tauw = tauw_over_rho*rho

	return tauw, terms
