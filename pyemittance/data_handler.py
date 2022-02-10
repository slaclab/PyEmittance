# module containing function to add/remove points to quad emit scan
import numpy as np
# on sim
# from beam_io_sim import get_beamsizes
# on lcls
# from beam_io import get_beamsizes, setquad, quad_control

# TODO: implement class

def check_symmetry(x, y):
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

    # return the x-vals to add
    # x_add are in same units as passed to fn
    return add_to_side, x_add

def add_measurements(add_to_side, x_add, x, y, y_err, axis, bs_fn=None):
    """Add beamsize measurements on left or right side based on
    symmetry of scan curve.
    x_add are the quad scan values k in units of 1/m^2"""

    # set indices for get_beamsizes fn data
    # get_beamsizes fn returns [xsize, ysize, xerr, yerr]
    ax_idx_size = 1 if axis == "y" else 0
    ax_idx_err = 3 if axis == "y" else 2
    # y-axis has different field polarity
    sign = -1 if axis == "y" else 1

    y_add, y_err_add = [], []
    # get new beamsizes from machine
    for ele in x_add:
        # this takes B in kG not K
        setquad(sign * get_quad_field(ele))
        # wait for magnet to settle
        time.sleep(3)
        # get new beamsizes at each quad val
        beamsizes = get_beamsizes()
        y_add.append(beamsizes[ax_idx_size])
        y_err_add.append(beamsizes[ax_idx_err])

    # then append to existing dataset
    if add_to_side == "left":
        new_x_list = list(x_add) + list(x)
        new_y_list = list(y_add) + list(y)
        new_y_err_list = list(y_err_add) + list(y_err)
    else:
        new_x_list = list(x) + list(x_add)
        new_y_list = list(y) + list(y_add)
        new_y_err_list = list(y_err) + list(y_err_add)

    # returns values in units of 1/m^2
    return np.array(new_x_list), np.array(new_y_list), np.array(new_y_err_list)

def find_inflection_pnt(x, y, show_plots=True):
    """Find inflection points in curve and remove
    points outside of convex region around min"""

    y = np.array(y)**2 # since we are fitting sizes**2

    # compute second derivative
    y_d2 = np.gradient(np.gradient(y))
    # find switching points
    infls = np.where(np.diff(np.sign(y_d2)))[0]

    if len(infls)==0:
        # No turning points found
        return None, None

    if len(infls)==1:
        infls = int(infls)
        if infls == np.argmin(y):
            # not sure what else to do here
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
        if len(infls)>2:
            #print("old ", infls)
            infls = list(infls)
            # cut it down to 2 closest to min
            if np.argmin(y) in infls:
                infls.remove(np.argmin(y))
            pnt1 = min(infls, key=lambda x: abs(x - np.argmin(y)))
            infls.remove(pnt1)
            pnt2 = min(infls, key=lambda x: abs(x - np.argmin(y)))
            # make new list of infls pnts
            infls = np.array([pnt1, pnt2])
            #print("new ", infls)

        # find min and max limits for slicing
        if np.max(infls) > np.argmin(y) and np.min(infls) > np.argmin(y):
            right = np.min(infls)
            left = None
        elif np.max(infls) < np.argmin(y) and np.min(infls) < np.argmin(y):
            left = np.max(infls)
            right = None
        elif np.min(infls)<np.argmin(y) and np.max(infls)>np.argmin(y):
            left = np.min(infls)
            right = np.max(infls)
        elif np.min(infls)==np.argmin(y):
            left = np.min(infls) + 1
            right = np.max(infls)
        elif np.max(infls)==np.argmin(y):
            left = np.min(infls)
            right = np.max(infls) + 1
        else:
            print("Case not implemented. Keeping data as is.")
            print(f"infls: {infls}, min: {np.argmin(y)}")
            left = None
            right = None

    if show_plots:
        y_new, x_new = y[left:right], x[left:right]

        import matplotlib.pyplot as plt
        from numpy.polynomial import polynomial as P

        fig = plt.figure(figsize=(10,5))

        x_fit = np.linspace(np.min(x), np.max(x), 50)
        # original polynomial for visuals
        c, stats = P.polyfit(x, y, 2, full=True)
        plt.plot(x_fit, P.polyval(x_fit,c))
        # updated polynomial for visuals
        c, stats = P.polyfit(x_new, y_new, 2, full=True)
        plt.plot(x_fit, P.polyval(x_fit, c))
        # plot the location of each inflection point
        for i, infl in enumerate(infls, 1):
            plt.axvline(x=x[infl], color='k', label=f'Inflection Point {i}')
        plt.scatter(x, y)
        plt.scatter(x_new, y_new, color="blue", label="Use")
        plt.ylim(np.min(y)*0.7, np.max(y)*1.3)
        plt.legend()
        plt.show()
        plt.close()

    return left, right
