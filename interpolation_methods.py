import numpy as np


def bilinear(xc, yc, ic, jc, x, y, z, nx, ny):
    """Bilinear interpolation z(x, y) at (xc, yc)"""

    if xc <= x[0] or xc >= x[nx-1] or yc <= y[0] or yc >= y[ny-1]:
        return z[ic, jc]

    if xc >= x[ic] and yc >= y[jc]:
        x1 = x[ic]
        y1 = y[jc]
        x2 = x[ic+1]
        y2 = y[jc+1]
        z1 = z[ic, jc]
        z2 = z[ic+1, jc]
        z3 = z[ic+1, jc+1]
        z4 = z[ic, jc+1]
    elif xc >= x[ic] and yc < y[jc]:
        x1 = x[ic]
        y1 = y[jc-1]
        x2 = x[ic+1]
        y2 = y[jc]
        z1 = z[ic, jc-1]
        z2 = z[ic+1, jc-1]
        z3 = z[ic+1, jc]
        z4 = z[ic, jc]
    elif xc < x[ic] and yc >= y[jc]:
        x1 = x[ic-1]
        y1 = y[jc]
        x2 = x[ic]
        y2 = y[jc+1]
        z1 = z[ic-1, jc]
        z2 = z[ic, jc]
        z3 = z[ic, jc+1]
        z4 = z[ic-1, jc+1]
    elif xc < x[ic] and yc < y[jc]:
        x1 = x[ic-1]
        y1 = y[jc-1]
        x2 = x[ic]
        y2 = y[jc]
        z1 = z[ic-1, jc-1]
        z2 = z[ic, jc-1]
        z3 = z[ic, jc]
        z4 = z[ic-1, jc]

    t = (xc-x1)/(x2-x1)
    u = (yc-y1)/(y2-y1)

    zc = (1.-t)*(1.-u)*z1+t*(1.-u)*z2+t*u*z3+(1.-t)*u*z4
    return zc


def one_dimensional(xc, yc, ic, jc, x, y, z, nx, ny):
    """One dimensional interpolation z(x, y) at (xc, yc)"""

    from scipy import interpolate

    int_dist = 1

    if (ic-int_dist < 0) or (ic+int_dist+1 > nx-1) or\
       (jc-int_dist < 0) or (jc+int_dist+1 > ny-1):
        return z[ic, jc]

    xi = x[ic-int_dist:ic+int_dist+1]
    yi = y[jc-int_dist:jc+int_dist+1]
    zi = z[ic-int_dist:ic+int_dist+1, jc-int_dist:jc+int_dist+1]

    zy = np.empty_like(yi)
    for j in range(len(yi)):
        fx = interpolate.interp1d(xi, zi[:, j])
        zy[j] = fx(xc)
    fy = interpolate.interp1d(yi, zy)
    return fy(yc)


def bilinear_div(xc, yc, ic, jc, x, y, z, nx, ny):
    """Bilinear interpolation z(x, y) at (xc, yc)"""

    if xc <= x[1] or xc >= x[nx-2] or yc <= y[1] or yc >= y[ny-2]:
        return z[ic, jc]

    if xc >= x[ic] and yc >= y[jc]:
        x1 = x[ic]
        y1 = y[jc]
        x2 = x[ic+1]
        y2 = y[jc+1]
        z1 = z[ic, jc]
        z2 = z[ic+1, jc]
        z3 = z[ic+1, jc+1]
        z4 = z[ic, jc+1]
    elif xc >= x[ic] and yc < y[jc]:
        x1 = x[ic]
        y1 = y[jc-1]
        x2 = x[ic+1]
        y2 = y[jc]
        z1 = z[ic, jc-1]
        z2 = z[ic+1, jc-1]
        z3 = z[ic+1, jc]
        z4 = z[ic, jc]
    elif xc < x[ic] and yc >= y[jc]:
        x1 = x[ic-1]
        y1 = y[jc]
        x2 = x[ic]
        y2 = y[jc+1]
        z1 = z[ic-1, jc]
        z2 = z[ic, jc]
        z3 = z[ic, jc+1]
        z4 = z[ic-1, jc+1]
    elif xc < x[ic] and yc < y[jc]:
        x1 = x[ic-1]
        y1 = y[jc-1]
        x2 = x[ic]
        y2 = y[jc]
        z1 = z[ic-1, jc-1]
        z2 = z[ic, jc-1]
        z3 = z[ic, jc]
        z4 = z[ic-1, jc]

    t = (xc-x1)/(x2-x1)
    u = (yc-y1)/(y2-y1)

    zc = (1.-t)*(1.-u)*z1+t*(1.-u)*z2+t*u*z3+(1.-t)*u*z4
    return zc

