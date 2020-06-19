from cubevis.utils import CubeVisArgumentError
import numpy
import scipy.optimize
import math

def tau(ratio, isoratio=65, niter=500,
        tol=0.000001, guess=0.04):
    """
    Given a ratio of the integrated intensity of two isotopes
    and the isotopic abundance ratio of those isotopes (the more
    abundant isotope should be in the numerator, this function
    returns the optical depth of the rarer isotope.

    It solves the equation,
      - ratio = (1-exp(-isoratio*S{tau}13))/(1-exp(-S{tau}13))
    for S{tau}13 using the newton raphson method.

    @param ratio: Ratio of the integrated intensity of more abundant
        isotope to less abundant isotope
    @type ratio: float
    @param isoratio: Ratio of abundance of more abundant isotope to
        less abundant isotope
    @type isoratio: float
    @param niter: Number of iterations to try before giving up
    @type niter: int
    @param tol: Tolerance for tau
    @type tol: float
    @return: the optical depth (tau) of rarer isotope. Returns -1 if
         it does not converge
    """
    f = lambda t, R, r: (1-R)+R*math.exp(-t)-math.exp(-r*t)
    fprime = lambda t, R, r: -R*math.exp(-t)+r*math.exp(-r*t)
    if ratio/isoratio > 0.8:
        #guess small tau
        if guess > 0.01:
            guess = 0.01

    best_tau = scipy.optimize.newton(f, guess, fprime, args=(ratio, isoratio),
                                     tol=tol, maxiter=niter)
    if best_tau < 0:
        best_tau = -1
    return best_tau


def tau_simple(ratio, isoratio=65, niter=500,
               tol=0.000001, guess=0.04,
               debug=False):
    """
    Given a ratio of the integrated intensity of two isotopes
    and the isotopic abundance ratio of those isotopes (the more
    abundant isotope should be in the numerator, this function
    returns the optical depth of the rarer isotope.

    It solves the equation,
      - ratio = (1-exp(-isoratio*S{tau}13))/(1-exp(-S{tau}13))
    for S{tau}13 using the newton raphson method.

    @param ratio: Ratio of the integrated intensity of more abundant
        isotope to less abundant isotope
    @type ratio: float
    @param isoratio: Ratio of abundance of more abundant isotope to
        less abundant isotope
    @type isoratio: float
    @param niter: Number of iterations to try before giving up
    @type niter: int
    @param tol: Tolerance for tau
    @type tol: float
    @return: the optical depth (tau) of rarer isotope. Returns -1 if
         it does not converge
    """
    x = guess
    for i in range(niter):
        x1 = x
        a = -isoratio * x
        b = -(isoratio + 1) * x
        c = -2 * x
        d = -x
        try:
            f = ((1-math.exp(a))/(1-math.exp(d))) - ratio
            f1 = ((isoratio * math.exp(a)) - ((isoratio-1)*math.exp(b)) - math.exp(d))/(1.-(2.*math.exp(d))+math.exp(c))
            delta = f/f1                  
            x = x1 - delta
            check = tol * x
            if abs(delta) <= abs(check):
                if debug:
                    print("Convergence in %d iterations" % i)
                if x <= 0.0:
                    if debug:
                        print("tau <= 0")
                    return -1
                if debug:
                    print("normal return")
                return x
        except:
            if debug:
                print("Float point error after %d iterations" % i)
            return -1
    print("Exceeded maxiter %d" % niter)
    return -1


