# module containing function to add/remove points to quad emit scan
import numpy as np
from scipy.optimize import curve_fit
from pyemittance.optics import get_k1, get_gradient, get_quad_field
# TODO: import m_0
m_0  = 0.000511

def adapt_range(x, y, axis, w=None, energy=0.135,
                cutoff_percent=0.3, num_points=5,
                save_plot=False, show_plots=True,
                verbose=False):
    """Returns new scan quad values AS LIST"""
    x = np.array(x)
    y = np.array(y)
    if w is not None:
        w = np.array(w)

    # Remove points that are too large
    # Do this only if bs more than double from min
    bs_range = y.max() - y.min()

    if bs_range / y.min() > 1:
        cutoff_lim = y.min() + bs_range * cutoff_percent
        idx = np.argwhere(y < cutoff_lim).flatten()

        if len(idx) >= 3:
            x = x[idx]
            y = y[idx]
            if w is not None:
                w = w[idx]

    # Set weights for polyfit (here the weight is sigma)
    if w is not None:
        w = 2 * w * y
    else:
        # Take weight as beamsize
        w = y

    # quad vals are passed in machine units
    x = get_k1(get_gradient(x), energy=energy, m_0=m_0)

    # y-dimensions has opposite polarity
    sign = -1 if axis == "y" else 1
    x = sign * x
    min_x, max_x = x.min(), x.max()

    # we need a poly fit here to find roots, poly_ylim, etc
    fit_coefs, fit_cov = curve_fit(func, x, y ** 2, sigma=w, absolute_sigma=True, method='trf')

    # space over which to fit
    x_fit = np.linspace(min_x, max_x, 100)

    # no more restrictions on quad vals, just staying within
    # the region already scanned (can increase this if need be)
    min_x_range, max_x_range = min([min_x, x[np.argmin(y)] - 2.5]), max([max_x, x[np.argmin(y)] + 2.5])

    c2, c1, c0 = fit_coefs

    if c2 < 0: # same if s11q is negative
        if verbose:
            print("Adjusting concave poly.")
        # go to lower side of concave polynomials
        # (assuming it is closer to the local minimum)
        x_min_concave = x[np.argmin(y)]
        # find the direction of sampling to minimum
        if (x[np.argmin(y)] - x[abs(np.argmin(y) - 2)]) < 0:
            x_max_concave = min_x_range
        else:
            x_max_concave = max_x_range

        x_fine_fit = np.linspace(x_min_concave, x_max_concave, num_points)
        return [sign * get_quad_field(ele) for ele in x_fine_fit]

    # find range within 2-3x the focus size
    # cutoff = 1.2-1.3 for lcls, 2 last MD
    # cutoff = 4 for facet and surrogate
    cutoff = 2
    # from data
    y_lim = y.min()** 2 * cutoff

    # from polyfit
    y_min_poly = np.polyval(fit_coefs, x_fit).min()
    if y_lim < y_min_poly:
        # in this case the roots won't exist
        y_lim = y_min_poly * cutoff

    if y_lim < 0:
        print(f"{axis} axis: min. of poly fit is negative. Setting it to a small val.")
        y_lim = np.mean(y ** 2) / 5

    roots = np.roots((c2, c1, c0 - y_lim))

    # Flag bad fit with complex roots
    if np.iscomplex(roots).any():
        raise ComplexRootError("Cannot adapt quad ranges, complex root encountered.")

    # if roots are outside quad scanning range, set to scan range lim
    if roots.min() < min_x_range:
        roots[np.argmin(roots)] = min_x_range
    if roots.max() > max_x_range:
        roots[np.argmax(roots)] = max_x_range

    if axis=="x":
        x_fine_fit = np.linspace(roots.min(), roots.max(), num_points)
    else:
        # instead of reversing array later
        x_fine_fit = np.linspace(roots.max(), roots.min(), num_points)

    # return the new quad measurement range for this axis (in kG!!)
    return [sign * get_quad_field(ele) for ele in x_fine_fit]


def check_symmetry(x, y, y_err, axis, bs_fn=None, add_meas=False):
    """Check symmetry of quad scan around min of scan
    and find side (left or right) and num of beamsize
    points to add to get a full curve"""

    if len(y) != len(x):
        raise Exception('Array lengths do not match!')

    # get number of points on left and right side
    left_side = np.argmin(y)
    right_side = len(x) - left_side - 1
    stepsize = abs((x[0] - x[-1]) / len(x))

    if left_side == right_side:
        return None

    elif left_side > right_side:
        add_to_side = "right"
        diff = left_side - right_side
        # add points to right_side
        xmin = x[-1] + stepsize
        xmax = xmin + diff * stepsize
        x_add = np.linspace(xmin, xmax, diff)

    elif right_side > left_side:
        add_to_side = "left"
        diff = right_side - left_side
        # add points to left_side
        xmin = x[0] - diff * stepsize
        xmax = x[0] - stepsize
        x_add = np.linspace(xmin, xmax, diff)

    if add_meas:
        return add_measurements(add_to_side, x_add, x, y, y_err, axis, bs_fn)
    else:
        return add_to_side, x_add


def add_measurements(add_to_side, x_add, x, y, y_err, axis, bs_fn=None):
    """Add beamsize measurements on left or right side based on
    symmetry of scan curve.
    x_add are the quad scan values k in units of 1/m^2"""

    if not bs_fn:
        return None
    
    # get new data points
    idx_size = 1 if axis == "y" else 0
    idx_err = 3 if axis == "y" else 2
    new_data = bs_fn(x_add)
    y_add, y_err_add = new_data[idx_size], new_data[idx_err]
        
    # then append to existing dataset
    if add_to_side == "left":
        new_x_list = list(x_add) + list(x)
        new_y_list = list(y_add) + list(y)
        new_y_err_list = list(y_err_add) + list(y_err)
    else:
        new_x_list = list(x) + list(x_add)
        new_y_list = list(y) + list(y_add)
        new_y_err_list = list(y_err) + list(y_err_add)
            
    return new_x_list, new_y_list, new_y_err_list


def find_inflection_pnt(x, y, show_plots=True):
    """Find inflection points in curve and remove
    points outside of convex region around min"""

    y = np.array(y)**2 # since we are fitting sizes**2

    # compute second derivative
    y_d2 = np.gradient(np.gradient(y))

    # find switching points
    infls = np.where(np.diff(np.sign(y_d2)))[0] + 1

    if len(infls)==0:
        # No turning points found
        return None, None

    if len(infls)==1:
        infls = int(infls)
        if infls == np.argmin(y):
            # if the only pnt is the min, don't do anything
            return None, None
        # check if point is on left or right of min
        if infls < np.argmin(y):
            left = infls
            right = None
        elif infls > np.argmin(y):
            left = None
            right = infls
        infls = np.array([infls])


    elif len(infls)>1:
        # cut it down to 2 closest to min
        if ( np.argmin(y) in infls and len(infls) > 2 ) or ( np.argmin(y) not in infls ):
            if np.argmin(y) in infls:
                infls = list(infls)
                infls.remove(np.argmin(y))
                infls = np.array(infls)

            # pick closest point from left side of min
            idx1 = np.argwhere(infls < np.argmin(y)).flatten()
            # pick closest point from right side of min
            idx2 = np.argwhere(infls > np.argmin(y)).flatten()

            if idx1.size == 0:
                left = None
            else:
                left = min(infls[idx1], key=lambda x: abs(x - np.argmin(y)))

            if idx2.size == 0:
                right = None
            else:
                right = min(infls[idx2], key=lambda x: abs(x - np.argmin(y)))
            # make new list of infls pnts
            infls = np.array([left, right])

        elif np.argmin(y) in infls and np.min(infls) == np.argmin(y):
            left = np.min(infls) - 1  # here - not +!
            right = np.max(infls)

        elif np.argmin(y) in infls and np.max(infls) == np.argmin(y):
            left = np.min(infls)
            right = np.max(infls) + 1

        else:
            print("Case not implemented. Keeping data as is.")
            print(f"infls: {infls}, min: {np.argmin(y)}")
            return None, None

    # if we end up with less than 3 points, don't do anything
    # not sure this would ever happen?
    if len(x[left:right]) < 3:
        return None, None

    if show_plots:
        y_new, x_new = y[left:right], x[left:right]

        import matplotlib.pyplot as plt
        from numpy.polynomial import polynomial as P

        x_fit = np.linspace(np.min(x), np.max(x), 50)
        # original polynomial for visuals
        c, stats = P.polyfit(x, y, 2, full=True)
        plt.plot(x_fit, P.polyval(x_fit,c))
        # updated polynomial for visuals
        c, stats = P.polyfit(x_new, y_new, 2, full=True)
        plt.plot(x_fit, P.polyval(x_fit, c), color="blue", linestyle="--")
        # remove nones from infls
        infls  = filter(None, infls)
        # plot the location of each inflection point
        for i, infl in enumerate(infls, 1):
            plt.axvline(x=x[infl], color='gray', label=f'Inflection Point {i}', linestyle="--" )
        plt.scatter(x, y)
        plt.scatter(x_new, y_new, color="blue", label="Use")
        plt.ylim(None, np.max(y)*1.3)
        plt.legend()
        plt.show()
        plt.close()

    return left, right

def func(x, a, b, c):
    """
    Second degree polynomial function of the form
    f(x) = ax^2 + bx + c
    :param x: input variable
    :param a: second deg coeff
    :param b: first deg coeff
    :param c: zeroth deg coeff
    :return: f(x)
    """
    return a*x**2 + b*x + c

class ComplexRootError(Exception):
    """
    Raised when the adapted range emit
    fit results in polynomial with complex root(s)
    """
    pass