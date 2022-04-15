"""
RCR.PY: Performs calculations for RCR model

Author:       Brian Krauth
              Department of Economics
              Simon Fraser University
Usage:

       rcr.py [infile outfile logfile]

       where

           infile      An optional argument giving the name
                       of the input file.  Default is IN.TXT.

           outfile     An optional argument giving the name of
                       the output file. Default is OUT.TXT.

           logfile     An optional argument giving the name of
                       the output file. Default is LOG.TXT.

       The RCR program will read in the INFILE, perform
       the calculations, and then write the results to OUTFILE.
       The program may also report information on its status to
       LOGFILE.

"""
# Standard library imports
import sys
# from datetime import datetime

# Third party imports
import numpy as np
from numpy.linalg import inv
# import pandas as pd

# Local application imports
from rcrutil import get_command_arguments, read_data, \
                    write_results, write_to_logfile, warn, die


def estimate_model(moment_vector, lambda_range):
    """Estimate the RCR model.

    Parameters
    ----------
    moment_vector : ndarray of floats
        its elements will be interpreted as the upper triangle of the
        (estimated) second moment matrix E(W'W), where W = [1 X Y Z].
        It is normally constructed by Stata.
    lambda_range : ndarray of floats
        its elements lambda values to consider

    Returns
    -------
    result_matrix : ndarray
        an array of parameter estimates and gradients

    Side effects
    ------------
    None.

    See also
    --------
    To be added.

    Notes
    -----
    To be added.

    Examples
    --------
    To be added.

    """
    write_to_logfile("Estimating model.\n")
    result_matrix = np.full((len(lambda_range) + 3,
                             len(moment_vector) + 1),
                            float('nan'))
    # Check to make sure the moments are consistent
    valid, identified = check_moments(moment_vector)
    # If moments are invalid, just stop there
    if not valid:
        return result_matrix
    # If model is not identified, just stop there
    # TODO: some model elements may still be identified here
    elif not identified:
        return result_matrix
    # We have closed forms for the global parameters lambdastar, thetastar,
    # and lambda(0), so we just estimate them directly.
    result_matrix[0, ] = estimate_parameter(lambdastar_fun, moment_vector)
    result_matrix[1, ] = estimate_parameter(thetastar_fun, moment_vector)
    result_matrix[2, ] = estimate_parameter(lambda0_fun, moment_vector)
    # Here we get to the main estimation problem.  We need to find the range
    # of theta values consistent with the lambda(theta) function falling in
    # lambda_range.  We have a closed form solution for lambda(theta), but
    # finding its inverse is an iterative problem.
    #
    # STEP 1: Estimate THETA_SEGMENTS, which is a global real vector
    #         indicating all critical points (i.e., points where the
    #         derivative is zero or nonexistent) of the function
    #         lambda(theta).  The function is continuous and monotonic
    #         between these points. Note that we don't know a priori how many
    #         critical points there will be, and so we don't know how big
    #         THETA_SEGMENTS will be.
    theta_segments = estimate_theta_segments(moment_vector)
    # STEP 2: For each row of lambda_range (i.e., each pair of lambda values):
    # do i=1,size(lambda_range,1)
    # j is the row in result_matrix corresponding to lambda_range(i,:)
    # j = 2+2*i
    #  Estimate the corresponding theta range, and put it in result_matrix
    result_matrix[3:5, :] = estimate_theta(moment_vector,
                                           lambda_range,
                                           theta_segments)
    return result_matrix


def estimate_theta_segments(moment_vector):
    """Divide real line into segments over which lambda(theta) is monotonic"""
    imax = 30000   # A bigger number produces an FP overflow in fortran
    sm = simplify_moments(moment_vector)
    thetastar = thetastar_fun(moment_vector)
    # THETAMAX is the largest value of theta for which we can calculate both
    # lambda(theta) and lambda(-theta) without generating a floating point
    # exception.
    thetamax = np.sqrt(sys.float_info.max / max(1.0, sm[4], sm[1] - sm[4]))
    # The calculation above seems clever, but it turns out not to always work.
    # So I've put in a hard limit as well
    thetamax = min(1.0e100, thetamax)
    # Create a starting set of theta values at which to calculate lambda
    thetavec = np.sort(np.append(np.linspace(-50.0, 50.0, imax - 2),
                                 (thetamax, -thetamax)))
    if (np.isfinite(thetastar)):
        # Figure out where thetastar lies in thetavec
        i = np.sum(thetavec < thetastar)
        # If i=0 or i=k, then thetastar is finite but outside of
        # [-thetamax,thetamax]. This is unlikely, but we should check.
        if ((i > 0) and (i < imax)):
            # Adjust i to ensure that -thetamax and thetamax are still
            # included in thetavec
            i = min(max(i, 2), imax - 2)
            # Replace the two elements of thetavec that bracket thetastar
            # with two more carefully-chosen numbers.  See BRACKET_THETA_STAR
            # for details
            thetavec[i-1: i+1] = bracket_theta_star(moment_vector)
            # TODO: There is a potential bug here.  The bracket_theta_star
            # function is used to take the two values in thetavec that are
            # closest to thetastar and replace them with values that are
            # guaranteed to give finite and nonzero lambda.  But there's
            # nothing to guarantee that these are still the two values in
            # thetavec that are the closest to thetastar.
            assert thetavec[i-2] < thetavec[i-1]
            assert thetavec[i] < thetavec[i+1]
        else:
            msg = "thetastar (={0}) > thetamax (={1}).".format(thetastar, thetamax)
            warn(msg)
    # Re-sort thetavec
    thetavec = np.sort(thetavec)
    # Calculate lambda for every theta in thetavec
    lambdavec = lambdafast(thetavec, simplify_moments(moment_vector))
    # If a detail_file has been specified, output thetavec and lambdavec to
    # that file
    if (len(detail_file) > 0):
        try:
            with open(detail_file, mode="w") as df:
                df.write("theta, lambda \n")
                for i in range(0, len(thetavec)):
                    df.write("{0}, {1} \n".format(thetavec[i], lambdavec[i]))
        except:
            warn("Unable to write to detail file {0}.".format(detail_file))
    # LOCALMIN = True if the corresponding element of THETAVEC appears to be
    # a local minimum
    localmin = ((lambdavec[1:imax-1] < lambdavec[0:imax-2]) &
                (lambdavec[1:imax-1] < lambdavec[2:imax]))
    # The end points are not local minima
    localmin = np.append(np.insert(localmin, [0], [False]), False)
    # LOCALMAX = True if the corresponding element of THETAVEC appears to be
    # a local maximum
    localmax = ((lambdavec[1:imax-1] > lambdavec[0:imax-2]) &
                (lambdavec[1:imax-1] > lambdavec[2:imax]))
    # The end points are not local max`ima
    localmax = np.append(np.insert(localmax, [0], [False]), False)
    # Figure out where THETASTAR lies in THETAVEC.  We need to do this
    # calculation again because we sorted THETAVEC
    if (np.isfinite(thetastar)):
        i = np.sum(thetavec < thetastar)
        if ((i > 0) and (i < imax)):
            # The two values bracketing THETASTAR are never local optima
            localmin[i-1:i+1] = False
            localmax[i-1:i+1] = False
    # Right now, we only have approximate local optima.  We need to apply
    # an iterative optimization algorithm to improve the precision.
    # do j=1,size(localmin)
    for j in range(1, len(localmin)):
        if localmin[j-1]:
            thetavec[j-1] = brent(thetavec[j-2],
                                  thetavec[j-1],
                                  thetavec[j],
                                  lambdafast,
                                  1.0e-10,
                                  simplify_moments(moment_vector))
        elif localmax[j-1]:
            thetavec[j-1] = brent(thetavec[j-2],
                                  thetavec[j-1],
                                  thetavec[j],
                                  negative_lambdafast,
                                  1.0e-10,
                                  simplify_moments(moment_vector))
    # Now we are ready to create THETA_SEGMENTS.
    if (np.isfinite(thetastar) and (i > 0) and (i < imax)):
        # THETA_SEGMENTS contains the two limits (-Inf,+Inf), the pair of
        # values that bracket thetastar, and any local optima
        theta_segments = np.append(np.concatenate([thetavec[i-1:i+1],
                                                   thetavec[localmin],
                                                   thetavec[localmax]]),
                                   (-thetamax, thetamax))
    else:
        # If thetastar is not finite, then we have two less elements in
        # THETA_SEGMENTS
        theta_segments = np.concatenate([thetavec[i-1:i+1],
                                         thetavec[localmin],
                                         thetavec[localmax]])
    # Sort the result (definitely necessary)
    theta_segments = np.sort(theta_segments)
    return theta_segments


def negative_lambdafast(theta, simplifiedMoments):
    return -lambdafast(theta, simplifiedMoments)


def brent(ax, bx, cx, func, tol, xopt):
    """Optimize by Brent algorithm"""
    itmax = 1000
    cgold = 0.3819660
    zeps = 1.0e-3 * np.finfo(float).eps  # NOT SURE THIS WILL WORK
    a = min(ax, cx)
    b = max(ax, cx)
    v = bx
    w = v
    x = v
    e = 0.0
    fx = func(x, xopt)
    fv = fx
    fw = fx
    # NOTE: I've added the extraneous line below so that code-checking
    # tools do not flag the "e = d" statement below as referencing
    # a nonexistent variable. In practice, this statement will never
    # be reached in the first loop iteration, after which point d will be
    # defined.
    d = e
    for iter in range(1, itmax + 1):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1
        if (abs(x - xm) <= (tol2 - 0.5 * (b - a))):
            brent = x
            break
        if (abs(e) > tol1):
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if (q > 0.0):
                p = -p
            q = abs(q)
            etemp = e
            e = d     # See NOTE above
            if (abs(p) >= abs(0.5 * q * etemp)) or \
               (p <= q * (a - x)) or \
               (p >= q * (b - x)):
                e = (a - x) if (x >= xm) else (b - x)
                d = cgold * e
            else:
                d = p / q
                u = x + d
                if (u - a < tol2) or (b - u < tol2):
                    d = tol1 * np.sign(xm - x)
        else:
            e = (a - x) if (x >= xm) else (b - x)
            d = cgold * e
        u = (x + d) if (abs(d) >= tol1) else (x + tol1 * np.sign(d))
        fu = func(u, xopt)
        if (fu <= fx):
            if (u >= x):
                a = x
            else:
                b = x
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu
        else:
            if (u < x):
                a = u
            else:
                b = u
            if (fu <= fw) or (w == x):
                v = w
                fv = fw
                w = u
                fw = fu
            elif (fu <= fv) or (v == x) or (v == w):
                v = u
                fv = fu
    if (iter == itmax):
        brent = x
        write_to_logfile("brent: exceed maximum iterations.\n")
    return brent


def bracket_theta_star(moment_vector):
    """Find theta valus close to thetastar"""
    # Get the value of THETASTAR.  If we are in this function it should be
    # finite.
    theta_star = thetastar_fun(moment_vector)
    # Get the limit of lambda(theta) as theta approaches THETASTAR,from below
    # and from above. These limits are generally not finite.
    sm = simplify_moments(moment_vector)
    # We may want to use np.inf here
    # NOTE: the np.sign seems extraneous here.
    true_limit = (np.array((1.0, -1.0)) *
                  np.sign(sm[2] - sm[5] * sm[1]/sm[4]) *
                  sys.float_info.max)
    # Pick a default value
    bracket = theta_star + np.array((-1.0, 1.0))*max(abs(theta_star), 1.0)*0.1
    j = 0
    for i in range(1, 101):
        # For the candidate bracket, consider THETASTAR plus or minus some
        # small number epsilon (epsilon gets smaller each iteration)
        candidate = (theta_star +
                     np.array((-1.0, 1.0)) * max(abs(theta_star), 1.0)*0.1**i)
        # To be a good bracket, candidate must satisfy some conditions:
        #    1. The bracket must be wide enough that the system can tell that
        #       CANDIDATE(1) < THETASTAR < CANDIDATE(2)
        #    2. The bracket must be narrow enough that lambda(candidate) is
        #       the same sign as true_limit.
        #    3. The bracket must be wide enough that lambda(candidate) is
        #       finite and nonzero. If candidate is very close to thetastar,
        #       then the calculated lambda(candidate) can be *either* NaN or
        #       zero.  The reason for this is that lambda(candidate) is a
        #       ratio of two things that are going to zero.  Approximation
        #       error will eventually make both the numerator and denominator
        #       indistingushable from zero (NaN), but sometimes the numerator
        #       will reach indistinguishable-from-zero faster (giving zero
        #       for the ratio).
        if ((candidate[0] < theta_star) and (candidate[1] > theta_star)):
            tmp2 = lambdafast(candidate, sm)
            if (np.isfinite(tmp2).all() and
               (tmp2[0]*np.sign(true_limit[0]) > 0.0) and
               (tmp2[1]*np.sign(true_limit[1]) > 0.0)):
                j = i
                bracket = candidate
            else:
                break
    if (j == 0):
        msg = "Unable to find a good bracket for thetastar"
        warn(msg)
    return bracket


def estimate_theta(moment_vector,
                   lambda_range,
                   theta_segments):
    """Estimate theta"""
    ntab = 10
    nmax = 10
    con = 1.4
    con2 = con * con
    big = sys.float_info.max
    safe = 2.0
    h = 1.0e-1
    errmax = 0.0
    estimate_theta = np.zeros((2, len(moment_vector)+1))
    deps = np.zeros(len(moment_vector))
    dmoments = np.zeros(len(moment_vector))
    a = np.zeros((ntab, ntab))
    fac = np.zeros(ntab-1)
    errt = np.zeros(ntab-1)
    # Get lambdastar and thetastar
    lambdastar = lambdastar_fun(moment_vector)
    thetastar = thetastar_fun(moment_vector)
    # Check to make sure that lambdastar is not in lambda_range.  If so,
    # theta is completely unidentified.
    if (lambda_range[0] <= lambdastar) and (lambdastar <= lambda_range[1]):
        estimate_theta[0, 0] = -np.inf
        estimate_theta[1, 0] = np.inf
        estimate_theta[:, 1:] = 0.0
        return estimate_theta
    # IMPORTANT_THETAS is a list of theta values for which lambda(theta) needs
    # to be calculated. We don't know in advance how many important values
    # there will be, so we make IMPORTANT_THETAS way too big, and initialize
    # it to all zeros (this choice is arbitrary).
    # Get simplified moments
    simplified_moments = simplify_moments(moment_vector)
    # k is the number of actual important theta values in IMPORTANT_THETAS
    important_thetas = np.array([])
    k = 1
    # Go piece by piece through theta_segments
    for i in range(1, len(theta_segments)):
        # Get the next pair of thetas.  This represents a range of thetas to
        # check
        current_theta_range = theta_segments[i-1:i+1]
        # Skip ahead to the next pair if thetastar is in the current range
        if ((not np.isfinite(thetastar)) or
           (current_theta_range[0] >= thetastar) or
           (current_theta_range[1] <= thetastar)):
            # Otherwise, calculate the range of lambdas associated with that
            # range of thetas
            current_lambda_range = lambdafast(current_theta_range,
                                              simplified_moments)
            # For each of the values in lambda_range
            for j in range(1, 3):
                # See if that value satisfies lambda(theta)-lambda(j)=0 for
                # some theta in current_theta_range
                if (lambda_range[j-1] > min(current_lambda_range)) and \
                   (lambda_range[j-1] < max(current_lambda_range)):
                    # If so, find theta such that lambda(theta)-lambda(j)=0
                    # and put it inour list of IMPORTANT_THETAS.  Of course,
                    # we can't quite find the exact theta.
                    tmp = zbrent(lambda_minus_lambda,
                                 current_theta_range[0],
                                 current_theta_range[1],
                                 1.0e-200,
                                 np.insert(simplified_moments,
                                           0,
                                           lambda_range[j-1]))
                    important_thetas = np.append(important_thetas, tmp)
                    k = k + 1
    # Add THETA_SEGMENTS to the list of IMPORTANT_THETAS
    important_thetas = np.append(important_thetas, theta_segments)
    # Add the OLS theta to the list of IMPORTANT_THETAS?
    # simplified_moments(3)/simplified_moments(2)
    # Calculate lambda(theta) for every theta in IMPORTANT_THETAS
    lambda_segments = lambdafast(important_thetas, simplified_moments)
    # INRANGE = True if a particular value of theta satisfies the condition
    #     lambda_range(1) <= lambda(theta) <= lambda_range(2)
    # Notice that we have put a little error tolerance in here, since
    # zbrent won't find the exact root.
    # TODO: Make sure the tolerance is big enough for the error in zbrent.
    inrange = ((lambda_segments >= lambda_range[0]-0.001) &
               (lambda_segments <= lambda_range[1]+0.001))
    if (k > 1):
        inrange[0:k-1] = True
    # If the lowest value in IMPORTANT_THETAS is in range, then there is no
    # (finite) lower bound
    if inrange[np.argmin(important_thetas)]:
        estimate_theta[0, 0] = -np.inf
    else:
        # Otherwise the the lower bound for theta is the minimum value in
        # IMPORTANT_THETAS that is in range
        estimate_theta[0, 0] = min(important_thetas[inrange])
    # If the highest value in IMPORTANT_THETAS is in range, then there is no
    # (finite) upper bound
    if inrange[np.argmax(important_thetas)]:
        estimate_theta[1, 0] = np.inf
    else:
        # Otherwise the the upper bound for theta is the maximum value in
        # IMPORTANT_THETAS that is in range
        estimate_theta[1, 0] = max(important_thetas[inrange])
    # Now we find the gradient
    # Take the gradient at both theta_L and theta_H
    for j in range(1, 3):
        theta = estimate_theta[j-1, 0]
        # The gradient can only be calculated if theta is finite!
        if np.isfinite(theta):
            # Gradients are estimated using a simple finite central difference:
            #            df/dx = (f(x+e)-f(x-e))/2e
            # where e is some small step size.  The tricky part is getting the
            # right step size.  The algorithm used here is an adaptation of
            # dfridr in Numerical Recipes.  However, that algorithm needs an
            # input initial step size h.
            #
            # http://www.fizyka.umk.pl/nrbook/c5-7.pdf: "As a function of input
            # h, it is typical for the accuracy to get better as h is made
            # larger, until a sudden point is reached where nonsensical
            # extrapolation produces early return with a large error. You
            # should therefore choose a fairly large value for h, but monitor
            # the returned value err, decreasing h if it is not small. For
            # functions whose characteristic x scale is of order unity, we
            # typically take h to be a few tenths."
            # So we try out starting values (h) until we get one that gives an
            # acceptable estimated error.
            for n in range(1, nmax+1):
                # Our candidate initial stepsize is 0.1, 0.001, ...
                h = 0.1 ** n
                # Initialize errmax
                errmax = 0.0
                # Initialize the finite-difference vector
                deps[:] = 0.0
                # First, we calculate the scalar-as-vector (dlambda / dtheta)
                # hh is the current step size.
                hh = h
                # Calculate an approximate derivative using stepsize hh
                a[0, 0] = (lambdafast(theta + hh, simplified_moments) -
                           lambdafast(theta - hh, simplified_moments)) / \
                          (2.0 * hh)
                # Set the error to very large
                err = big
                # Generate a geometric series
                fac[0:ntab - 1] = geop(con2, con2, ntab - 1)
                # Now we try progressively smaller stepsizes
                for k in range(2, ntab + 1):
                    # The new stepsize hh is the old stepsize divided by 1.4
                    hh = hh / con
                    # Calculate an approximate derivative with the new
                    # stepsize
                    a[0, k-1] = ((lambdafast(theta + hh, simplified_moments) -
                                  lambdafast(theta - hh, simplified_moments)) /
                                 (2.0 * hh))
                    # Then use Neville's method to estimate the error
                    for m in range(2, k + 1):
                        a[m - 1, k - 1] = ((a[m - 2, k - 1] *
                                            fac[m - 2] -
                                            a[m - 2, k - 2]) /
                                           (fac[m - 2] - 1.0))
                    errt[0:k - 2] = np.max((abs(a[1:k, k - 1] -
                                                a[0:k - 1, k - 1]),
                                            abs(a[1:k, k - 1] -
                                                a[0:k - 1, k - 2])))
                    ierrmin = np.argmin(errt[0:k - 1])
                    # If the approximation error is lower than any previous,
                    # use that value
                    if (errt[ierrmin] <= err):
                        err = errt[ierrmin]
                        dfridr = a[1 + ierrmin, k - 1]
                    if abs(a[k - 1, k - 1] - a[k - 2, k - 2]) >= (safe * err):
                        break
                # errmax is the biggest approximation error so far for the
                # current value of h
                errmax = max(errmax, err)
                # Now we have a candidate derivative dlambda/dtheta
                dtheta = dfridr
                # Second, estimate the vector (dlambda / dmoment_vector)
                for i in range(1, len(moment_vector) + 1):
                    hh = h
                    deps[i-1] = hh
                    a[0, 0] = ((lambdafun(moment_vector + deps, theta) -
                                lambdafun(moment_vector - deps, theta)) /
                               (2.0 * hh))
                    err = big
                    fac[0:ntab - 1] = geop(con2, con2, ntab - 1)
                    for k in range(2, ntab + 1):
                        hh = hh / con
                        deps[i-1] = hh
                        a[0, k - 1] = (lambdafun(moment_vector + deps,
                                                 theta) -
                                       lambdafun(moment_vector - deps,
                                                 theta)) / (2.0 * hh)
                        for m in range(2, k + 1):
                            a[m - 1, k - 1] = (a[m - 2, k - 1] * fac[m - 2] -
                                               a[m - 2, k - 2]) / \
                                               (fac[m - 2] - 1.0)
                        errt[0:k - 1] = np.max((abs(a[1:k, k - 1] -
                                                    a[0:k - 1, k - 1]),
                                                abs(a[1:k, k - 1] -
                                                    a[0:k - 1, k - 2])))
                        ierrmin = np.argmin(errt[0:k - 1])
                        if (errt[ierrmin] <= err):
                            err = errt[ierrmin]
                            dfridr = a[1 + ierrmin, k - 1]
                        if abs(a[k - 1, k - 1] - a[k - 2, k - 2]) >= \
                           (safe * err):
                            break
                    # errmax is the biggest approximation error so far for the
                    # current value of h
                    errmax = max(errmax, err)
                    dmoments[i - 1] = dfridr
                    deps[i - 1] = 0.0
                # At this point we have estimates of the derivatives stored in
                # dtheta and dmoments. We also have the maximum approximation
                # error for the current h stored in errmax. If that
                # approximation error is "good enough" we are done and can
                # exit the loop
                if (errmax < 0.01):
                    break
                # Otherwise we will try again with a smaller h
                if (n == nmax):
                    msg1 = "Inaccurate SE for thetaL/H."
                    msg2 = "Try normalizing variables to mean zero."
                    warn(msg1 + " " + msg2)
            # Finally, we apply the implicit function theorem to calculate the
            # gradient that we actually need:
            #   dtheta/dmoments = -(dlambda/dmoments)/(dlambda/dtheta)
            estimate_theta[j-1, 1:] = -dmoments / dtheta
        else:
            # If theta is infinite, then the gradient is zero.
            estimate_theta[j - 1, 1:] = 0.0
    return estimate_theta


def zbrent(func, x1, x2, tol, xopt):
    """Find a root using the Brent algorithm"""
    itmax = 1000
    eps = np.finfo(float).eps   # in fortran was epsilon(x1)
    a = x1
    b = x2
    fa = func(a, xopt)
    fb = func(b, xopt)
    if (((fa > 0.0) and (fb > 0.0)) or ((fa < 0.0) and (fb < 0.0))):
        write_to_logfile("Error in zbrent: Root is not bracketed")
        # call die("root must be bracketed for zbrent") # UPDATE
    c = b
    fc = fb
    for iter in range(1, itmax + 1):
        if (((fb > 0.0) and (fc > 0.0)) or ((fb < 0.0) and (fc < 0.0))):
            c = a
            fc = fa
            d = b - a
            e = d
        if (abs(fc) < abs(fb)):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        # check for convergence
        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        if (abs(xm) <= tol1) or (fb == 0.0):
            zbrent = b
            break
        if (abs(e) >= tol1) and (abs(fa) > abs(fb)):
            s = fb / fa
            if (a == c):
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if (p > 0.0):
                q = -q
            p = abs(p)
            if (2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))):
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d
        a = b
        fa = fb
        b = b + d if abs(d) > tol1 else tol1 * np.sign(xm)
        fb = func(b, xopt)
    if (iter == itmax):
        zbrent = b
        write_to_logfile("zbrent: exceeded maximum iterations")
    return zbrent


def lambda_minus_lambda(theta, simplified_moments_and_lambda):
    """Calculate lamba(thetaval)-lambdaval given thetaval and lambdaval"""
    lambda1 = lambdafast(theta, simplified_moments_and_lambda[1:])
    lambda0 = simplified_moments_and_lambda[0]
    return lambda1 - lambda0


def simplify_moments(moment_vector):
    """Convert moment_vector into the six moments needed for the model"""
    # Get sizes
    m = len(moment_vector)
    k = int(1 + np.floor((np.sqrt(1 + 8 * m) - 1) / 2))
    assert 2*(m + 1) == k ** 2 + k
    mvtmp = np.append(1.0, moment_vector)
    xtmp = np.zeros((k, k))
    # The array XTMP will contain the full cross-product matrix E(WW')
    # where W = [1 X Y Z]
    h = 0
    for i in range(h, k):
        for j in range(i, k):
            xtmp[i, j] = mvtmp[h]
            xtmp[j, i] = mvtmp[h]
            h = h + 1
    # The array XX will contain the symmetric matrix E(XX')
    XX = xtmp[0:(k - 2), 0:(k - 2)]
    # The array XY will contain the vector E(XY)
    XY = xtmp[(k - 2), 0:(k - 2)]
    # The array XZ will contain the vector E(XZ)
    XZ = xtmp[(k - 1), 0:(k - 2)]
    # Now we fill in simplify_moments with the various moments.
    simplify_moments = np.zeros(6)
    # varY
    simplify_moments[0] = (moment_vector[m - 3] -
                           (moment_vector[k - 3]) ** 2)
    # varZ
    simplify_moments[1] = (moment_vector[m - 1] -
                           (moment_vector[k - 2]) ** 2)
    # covYZ
    simplify_moments[2] = (moment_vector[m - 2] -
                           moment_vector[k - 2]*moment_vector[k - 3])
    try:
        # varYhat
        simplify_moments[3] = XY.T @ inv(XX) @ XY - (moment_vector[k - 3]) ** 2
        # varZhat
        simplify_moments[4] = XZ.T @ inv(XX) @ XZ - (moment_vector[k - 2]) ** 2
        # covYZhat
        simplify_moments[5] = ((XY.T @ inv(XX) @ XZ) -
                               (moment_vector[k - 2])*(moment_vector[k - 3]))
    except:
        die("FATAL ERROR: X'X matrix is singular")
    # When there is only one control variable, yhat and zhat are perfectly
    # correlated (positively or negatively) With rounding error, this can lead
    # to a correlation that is > 1 in absolute value.  This can create
    # problems, so we force the correlation to be exactly 1.
    # TODO: This also could happen if there is more than one control variable
    # but only one happens to have a nonzero coefficient.  I don't know how to
    # handle that case.
    if k == 4:
        simplify_moments[5] = (np.sign(simplify_moments[5]) *
                               np.sqrt(simplify_moments[3] *
                                       simplify_moments[4]))
    return simplify_moments


def check_moments(moment_vector):
    """Check to ensure moment_vector is valid"""
    sm = simplify_moments(moment_vector)
    # First make sure that moment_vector describes a valid covariance matrix
    valid = True
    if sm[0] < 0.0:
        valid = False
        warn("Invalid data: var(y) = {0} < 0".format(sm[0]))
    if sm[1] < 0.0:
        valid = False
        warn("Invalid data: var(z) = {0} < 0".format(sm[1]))
    if sm[3] < 0.0:
        valid = False
        warn("Invalid data: var(yhat) = {0} < 0".format(sm[3]))
    if sm[5] < 0.0:
        valid = False
        warn("Invalid data: var(zhat) = {0} < 0".format(sm[5]))
    if np.abs(sm[2]) > np.sqrt(sm[0] * sm[1]):
        valid = False
        covyz = np.abs(sm[2])
        sdyz = np.sqrt(sm[0] * sm[1])
        warn("Invalid data: |cov(y,z)| = {0} > {1} sqrt(var(y)*var(z))".format(covyz, sdyz))
    if np.abs(sm[5]) > np.sqrt(sm[3] * sm[4]):
        valid = False
        covyz = np.abs(sm[5])
        sdyz = np.sqrt(sm[3] * sm[4])
        warn("Invalid data: cov(yhat,zhat) = {0} > {1} sqrt(var(yhat)*var(zhat))".format(covyz,sdyz))
    # Next make sure that the identifying conditions are satisfied.
    # TODO: Maybe these could be addressed with warnings rather than error
    # messages?
    identified = valid
    if sm[0] == 0.0:
        identified = False
        warn("Model not identified: var(y) = 0")
    if sm[1] == 0.0:
        identified = False
        warn("Model not identified: var(z) = 0")
    if sm[3] == 0.0:
        identified = False
        warn("Model not identified: var(yhat) = 0")
    if sm[3] == sm[0]:
        identified = False
        warn("Model not identified: y is an exact linear function of X")
    # TODO: We may also want to check for var(zhat)=0.
    # The model is identified in this case, but we may need to take special
    # steps to get the calculations right.
    return valid, identified


def lambdastar_fun(moment_vector):
    """Calculate lambdastar"""
    sm = simplify_moments(moment_vector)
    # lambdastar is defined as sqrt( var(z)/var(zhat) - 1)
    # The check_moments subroutine should ensure that
    #   var(z) > 0 and that var(z) >= var(zhat) >= 0.
    # This implies that lambdastar >= 0.
    # Special values: If var(zhat) = 0, then lambdastar = +Infinity
    lambdastar = np.sqrt(np.maximum(sm[1] / sm[4], 1.0) - 1.0)
    return lambdastar


def thetastar_fun(moment_vector):
    """Calculate thetastar"""
    sm = simplify_moments(moment_vector)
    # thetastar is defined as
    #   cov(yhat,zhat)/var(zhat)
    # The check_moments subroutine should ensure that
    # var(zhat) >= 0 and that if var(zhat)=0 -> cov(yhat,zhat)=0.
    # Special values: If var(zhat)=0, then thetastar = 0/0 = NaN.
    thetastar = sm[5] / sm[4]
    return thetastar


def lambdafast(theta, simplifiedMoments):
    """Calculate lambda for the given array of thetas"""
    y = simplifiedMoments[0]
    z = simplifiedMoments[1]
    yz = simplifiedMoments[2]
    yhat = simplifiedMoments[3]
    zhat = simplifiedMoments[4]
    yzhat = simplifiedMoments[5]
    # Potential FPE
    lf_num = (yhat -
              2.0 * theta * yzhat +
              theta ** 2 * zhat)
    lf_denom = (y - yhat -
                (2.0) * theta * (yz - yzhat) +
                theta ** 2 * (z - zhat))
    lambdafast = np.asarray(lf_num / lf_denom)
    lambdafast[lambdafast < 0.0] = np.nan
    lambdafast = (yz - yzhat - theta * (z - zhat)) / \
                 (yzhat - theta * zhat) * np.sqrt(lambdafast)
    return lambdafast


def lambdafun(moment_vector, theta):
    """Calculate lambda for the given theta"""
    # Potential FPE
    lambdafun = lambdafast(theta, simplify_moments(moment_vector))
    return lambdafun


def lambda0_fun(moment_vector):
    """"Calculate lambda(theta) for theta = 0"""
    # lambda0 is defined as:
    # (cov(y,z)/cov(yhat,zhat)-1) / sqrt(var(y)/var(yhat)-1)
    # The check_moments subroutine should ensure that
    #  var(y) >= var(yhat) > 0, so the denominator is
    # always positive and finite.
    # Special values: If cov(yhat,zhat)=0, then lambda0 can
    #   be +Infinity, -Infinity, or NaN depending on the sign
    #   of cov(y,z).
    lambda0 = lambdafun(moment_vector, 0.0)
    return lambda0


def geop(first, factor, n):
    """Create a geometric series"""
    geop = np.zeros(n)
    if (n > 0):
        geop[0] = first
    for k in range(1, n-1):
        geop[k] = geop[k - 1] * factor
    return geop


def estimate_parameter(func, moment_vector):
    """Estimate a parameter and its gradient"""
    estimate_parameter = np.zeros(len(moment_vector) + 1)
    estimate_parameter[0] = func(moment_vector)
    nmax = 10
    ntab = 10
    con = 1.4
    con2 = con ** 2
    h = 1.0e-4
    safe = 2.0
    big = 1.0e300
    deps = np.zeros(len(moment_vector))
    errt = np.zeros(ntab - 1)
    fac = np.zeros(ntab - 1)
    a = np.zeros((ntab, ntab))
    if np.isfinite(estimate_parameter[0]):
        for n in range(1, nmax + 1):
            h = 0.1 ** n
            errmax = 0.0
            # We are estimating the gradient, i.e., a vector of derivatives
            # the same size as moment_vector
            for i in range(1, len(moment_vector) + 1):
                # Re-initialize DEPS
                deps[:] = 0.0
                # HH is the step size.  It is chosen by an algorithm borrowed
                # from the dfridr function in Numerical Recipes.  We start
                # with HH set to a predetermined value H.  After that, each
                # successive value of HH is the previous value divided by CON
                # (which is set to 1.4)
                hh = h
                # Set element i of DEPS to HH.
                deps[i - 1] = hh
                # Calculate the first approximation
                a[0, 0] = (func(moment_vector + deps) -
                           func(moment_vector - deps)) / (2.0*hh)
                dfridr = a[0, 0]   # WORKAROUND
                # The error is assumed to be a big number
                err = big
                # Generate a geometric series
                fac = geop(con2, con2, ntab - 1)
                # Try a total of NTAB different step sizes
                for j in range(2, ntab + 1):
                    # Generate the next step size
                    hh = hh / con
                    # Set DEPS based on that step size
                    deps[i - 1] = hh
                    # Calculate the approximate derivative for that step size
                    a[0, j - 1] = (func(moment_vector + deps) -
                                   func(moment_vector - deps)) / (2.0*hh)
                    # Next we estimate the approximation error for the current
                    # step size
                    for k in range(2, j + 1):
                        a[k - 1, j - 1] = (a[k - 2, j - 1] *
                                           fac[k - 2] -
                                           a[k - 2, j - 2]) / \
                                          (fac[k - 2] - 1.0)
                    errt[0:j - 1] = np.maximum(np.abs(a[1:j, j - 1] -
                                                      a[0:j - 1, j - 1]),
                                               np.abs(a[1:j, j - 1] -
                                                      a[0:j - 1, j - 2]))
                    ierrmin = np.argmin(errt[0:j - 1])
                    # If the error is smaller than the lowest previous error,
                    # use that hh
                    if (errt[ierrmin] <= err):
                        err = errt[ierrmin]
                        dfridr = a[1 + ierrmin, j - 1]
                    # If the error is much larger than the lowest previous
                    # error, stop
                    if np.abs(a[j - 1, j - 1] - a[j - 2, j - 2]) >= \
                       (safe * err):
                        break
                errmax = np.maximum(errmax, err)
                estimate_parameter[i] = dfridr
            if (errmax < 0.01):
                break
            if (n == nmax):
                msg1 = "Inaccurate SE for // fname //."
                msg2 = "Try normalizing variables to mean zero."
                warn(msg1 + msg2)
    else:
        estimate_parameter[1:] = 0.0   # change to internal_nan
    return estimate_parameter


#############################################################################
# Begin run code
#############################################################################

# Load in arguments from call to program
infile, outfile, detail_file = get_command_arguments(sys.argv)

# Read in the data from INFILE
# Side effect: allocation/creation of moment_vector, lambda_range, and
# result_matrix
(n_moments, n_lambda, external_big_number, moment_vector,
    lambda_range) = read_data(infile)

# Perform the calculations and put the results in result_matrix
# (side effect: allocation of theta_segments, writing to detail_file)
result_matrix = estimate_model(moment_vector, lambda_range)

# Write out the data to OUTFILE
# Issue: Python breaks lines above a certain length, we don't want to do that
write_results(result_matrix, outfile)

#############################################################################
# End run code
#############################################################################
