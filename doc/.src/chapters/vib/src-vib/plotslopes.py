# Taken from
# http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg18418.html

import numpy as np
import matplotlib.pyplot as plt

def slope_marker(origin, slope, size_frac=0.1, pad_frac=0.1, ax=None,
                 invert=False):
    """Plot triangular slope marker labeled with slope.

    Parameters
    ----------
    origin : (x, y)
        tuple of x, y coordinates for the slope
    slope : float or (rise, run)
        the length of the slope triangle
    size_frac : float
        the fraction of the xaxis length used to determine the size of the slope
        marker. Should be less than 1.
    pad_frac : float
        the fraction of the slope marker used to pad text labels. Should be less
        than 1.
    invert : bool
        Normally, the slope marker is below a line for positive slopes and above
        a line for negative slopes; `invert` flips the marker.
    """
    if ax is None:
        ax = plt.gca()

    if np.iterable(slope):
        rise, run = slope
        slope = float(rise) / run
    else:
        rise = run = None

    x0, y0 = origin
    xlim = ax.get_xlim()
    dx_linear = size_frac * (xlim[1] - xlim[0])
    dx_decades = size_frac * (np.log10(xlim[1]) - np.log10(xlim[0]))

    if invert:
        dx_linear = -dx_linear
        dx_decades = -dx_decades

    if ax.get_xscale() == 'log':
        log_size = dx_decades
        dx = _log_distance(x0, log_size)
        x_run = _text_position(x0, log_size/2., scale='log')
        x_rise = _text_position(x0+dx, dx_decades*pad_frac, scale='log')
    else:
        dx = dx_linear
        x_run = _text_position(x0, dx/2.)
        x_rise = _text_position(x0+dx, pad_frac * dx)

    if ax.get_yscale() == 'log':
        log_size = dx_decades * slope
        dy = _log_distance(y0, log_size)
        y_run = _text_position(y0, -dx_decades*slope*pad_frac, scale='log')
        y_rise = _text_position(y0, log_size/2., scale='log')
    else:
        dy = dx_linear * slope
        y_run = _text_position(y0, -(pad_frac * dy))
        y_rise = _text_position(y0, dy/2.)

    x_pad = pad_frac * dx
    y_pad = pad_frac * dy

    va = 'top' if y_pad > 0 else 'bottom'
    ha = 'left' if x_pad > 0 else 'right'
    if rise is not None:
        ax.text(x_run, y_run, str(run), va=va, ha='center')
        ax.text(x_rise, y_rise, str(rise), ha=ha, va='center')
    else:
        ax.text(x_rise, y_rise, str(slope), ha=ha, va='center')

    ax.add_patch(_slope_triangle(origin, dx, dy))


def log_displace(x0, dx_log=None, x1=None, frac=None):
    """Return point displaced by a logarithmic value.

    For example, if you want to move 1 decade away from `x0`, set `dx_log` = 1,
    such that for `x0` = 10, we have `displace(10, 1)` = 100

    Parameters
    ----------
    x0 : float
        reference point
    dx_log : float
        displacement in decades.
    x1 : float
        end point
    frac : float
        fraction of line (on logarithmic scale) between x0 and x1
    """
    if dx_log is not None:
        return 10**(np.log10(x0) + dx_log)
    elif x1 is not None and frac is not None:
        return 10**(np.log10(x0) + frac * np.log10(float(x1)/x0))
    else:
        raise ValueError('Specify `dx_log` or both `x1` and `frac`.')


def _log_distance(x0, dx_decades):
    return log_displace(x0, dx_decades) - x0

def _text_position(x0, dx, scale='linear'):
    if scale == 'linear':
        return x0 + dx
    elif scale == 'log':
        return log_displace(x0, dx)
    else:
        raise ValueError('Unknown value for `scale`: %s' % scale)


def _slope_triangle(origin, dx, dy, ec='none', fc='0.8', **poly_kwargs):
    """Return Polygon representing slope.
          /|
         / | dy
        /__|
         dx
    """
    verts = [np.asarray(origin)]
    verts.append(verts[0] + (dx, 0))
    verts.append(verts[0] + (dx, dy))
    return plt.Polygon(verts, ec=ec, fc=fc, **poly_kwargs)


if __name__ == '__main__':
    plt.plot([0, 2], [0, 1])
    slope_marker((1, 0.4), (1, 2))

    x = np.logspace(0, 2)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.loglog(x, x**0.5)
    slope_marker((10, 2), (1, 2), ax=ax1)

    ax2.loglog(x, x**-0.5)
    slope_marker((10, .4), (-1, 2), ax=ax2)

    ax3.loglog(x, x**0.5)
    slope_marker((10, 4), (1, 2), invert=True, ax=ax3)

    ax4.loglog(x, x**0.5)
    slope_marker((10, 2), 0.5, ax=ax4)

    plt.show()
