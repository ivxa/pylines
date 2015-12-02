# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np


def resamp_line(all_lines, resamp, a, CGS_units):
    """Downsamples the line to 'resamp' cells.
    All the lines should have the same lenght with a constant
    cell size.
    We are not averaging the data in the same bin because we are
    dealing with relativistic quantities. E.g. to average the
    density we should first average the number of particles, then
    average the volume in the lab frame and then compute again the
    density.
    If needed dx and dy of the new sampling can be computed as the
    projection of dl, i.e. dl*cos theta, where theta con be
    computed from vx and vy.
    After resampling you may notice that there is a gap at the top
    of the lines map, which is different for each line. The
    largest gap corresponds to the largest line (before resampling)
    and its size it's just the bin size. The reason is that we are
    keeping for each bin the first value of the bin, e.g. x[0], the
    opposite would occure (i.e. a gap at the bottom of the map) if
    we would keep the last value of the bin, e.g. x[-1]. The gap
    looks shorter for the other lines because in that case, the gap
    continues outside the grid (with 0 values)"""

    # Length of the largest line and maximum number of points
    sc_max = 0.
    N_max = 0
    for line_values in all_lines:
        xc = zip(*line_values)[0]
        yc = zip(*line_values)[1]
        sc = 0.
        for i in range(len(xc)-1):
            dxc = xc[i+1]-xc[i]
            dyc = yc[i+1]-yc[i]
            sc += np.sqrt(dxc**2.+dyc**2.)
        if sc > sc_max:
            sc_max = sc
        if len(xc) > N_max:
            N_max = len(xc)

    if N_max <= resamp:
        raise RuntimeError('Resampling value too big')

    # Cell size (fixed for all lines)
    dl = sc_max/resamp
    #print dl

    # Resampled line grid, with resamp+1 vertices
    l = [-1.]
    l_aux = 0.
    for k in range(resamp):
        l_aux += dl
        l += [l_aux]
    l[-1] += dl

    # Save line lenght info
    with open('L_dL.dat', 'w') as f:
        if CGS_units == 0:
            f.write(str(l[-1]-dl) + ' ' + str(dl) + '\n')
        else:
            f.write(str((l[-1]-dl)*a) + ' ' + str(dl*a) + '\n')

    all_lines_new = []
    # Regrouping of the physical parameters in the new grid
    for line_values in all_lines:
        xc, yc, ic, jc, densc, epsc, vxc, vyc, divc, tracerc, tc = zip(*line_values)

        # Line grid before resampling
        sc = [0.]
        sc_aux = 0.
        for i in range(len(xc)-1):
            dxc = xc[i+1]-xc[i]
            dyc = yc[i+1]-yc[i]
            sc_aux += np.sqrt(dxc**2.+dyc**2.)
            sc += [sc_aux]

        # Regrouping
        x = []
        y = []
        ii = []
        jj = []
        dens = []
        eps = []
        vx = []
        vy = []
        div = []
        tracer = []
        t = []
        line_values_new = []

        k = 0
        i_first = 0
        while True:
            # Look for the bin fiting the first element
            while l[k] < sc[i_first] <= l[k+1] == False:
                k += 1
            # Look for the other bin elements
            for i in range(i_first, len(xc)):
                if l[k] < sc[i] <= l[k+1]:
                    x += [xc[i]]
                    y += [yc[i]]
                    ii += [ic[i]]
                    jj += [jc[i]]
                    dens += [densc[i]]
                    eps += [epsc[i]]
                    vx += [vxc[i]]
                    vy += [vyc[i]]
                    div += [divc[i]]
                    tracer += [tracerc[i]]
                    t += [tc[i]]
                else:
                    break
            # Once the bin is full, save it as the first element
            x_new = x[0]#np.average(x)
            y_new = y[0]#np.average(y)
            ii_new = ii[0]#np.int(np.floor(np.average(ii)))
            jj_new = jj[0]#np.int(np.floor(np.average(jj)))
            dens_new = dens[0]#np.average(dens)
            eps_new = eps[0]#np.average(eps)
            vx_new = vx[0]#np.average(vx)
            vy_new = vy[0]#np.average(vy)
            div_new = div[0]#np.average(div)
            tracer_new = tracer[0]#np.average(tracer)
            t_new = t[0]#np.average(t)
            line_values_new.append([x_new, y_new, ii_new, jj_new,
                                    dens_new, eps_new,
                                    vx_new, vy_new, div_new,
                                    tracer_new, t_new])

            # Continue with the next bin
            i_first = i
            k += 1
            x = []
            y = []
            ii = []
            jj = []
            dens = []
            eps = []
            vx = []
            vy = []
            div = []
            tracer = []
            t = []

            if k >= len(l)-2 or i_first >= len(sc)-1:
                x += [xc[i]]
                y += [yc[i]]
                ii += [ic[i]]
                jj += [jc[i]]
                dens += [densc[i]]
                eps += [epsc[i]]
                vx += [vxc[i]]
                vy += [vyc[i]]
                div += [divc[i]]
                tracer += [tracerc[i]]
                t += [tc[i]]

                x_new = x[0]#np.average(x)
                y_new = y[0]#np.average(y)
                ii_new = ii[0]#np.int(np.floor(np.average(ii)))
                jj_new = jj[0]#np.int(np.floor(np.average(jj)))
                dens_new = dens[0]#np.average(dens)
                eps_new = eps[0]#np.average(eps)
                vx_new = vx[0]#np.average(vx)
                vy_new = vy[0]#np.average(vy)
                div_new = div[0]#np.average(div)
                tracer_new = tracer[0]#np.average(tracer)
                t_new = t[0]#np.average(t)
                line_values_new.append([x_new, y_new, ii_new, jj_new,
                                        dens_new, eps_new,
                                        vx_new, vy_new, div_new,
                                        tracer_new, t_new])
                break

        non_zero = len(line_values_new)
        for z in range(resamp-non_zero):
            x_new = 0
            y_new = 0
            ii_new = 0
            jj_new = 0
            dens_new = 0
            eps_new = 0
            vx_new = 0
            vy_new = 0
            div_new = 0
            tracer_new = 0
            t_new = 0
            line_values_new.append([x_new, y_new, ii_new, jj_new,
                                    dens_new, eps_new,
                                    vx_new, vy_new, div_new,
                                    tracer_new, t_new])
        all_lines_new.append(line_values_new)
    return all_lines_new
