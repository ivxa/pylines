# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np
import time
import interpolation_methods as im
from resamp_line import resamp_line


def buffer_present_line(xc, yc, ic, jc, densc, epsc, vxc, vyc,
                        divc, tracerc, tc, line_values):
    """Keeps track of all the line parameters"""

    line_values.append([xc, yc, ic, jc, densc, epsc,
                        vxc, vyc, divc, tracerc, tc])
    return


def buffer_all_lines(all_lines, line_values):
    """Keeps track of all the complete lines"""

    all_lines.append(line_values)
    return


def save_all_lines(output_file, all_lines):
    """Save all the lines to a file.
       Two lines for each current line:
       1) length of line
       2) (xc, yc, ic, jc,  densc, epsc, vxc, vyc, divc, tracerc, tc)
       times the length of the line"""

    with open(output_file + '.dat', 'w') as f:
        for line in all_lines:
            line_length = len(line)
            f.write(str(line_length) + '\n')
            for n in range(line_length):
                for var in line[n]:
                    f.write(str(var) + ' ')
            f.write('\n')
    return


def save_one_file_per_line(output_file, all_lines, gammaad, c, fB0, tr0,
                           sf0_array, excluded_lines, input_file):
    """Save one file per line. Each column:
    x, y, z, vx, vy, vz, dens, P, B, -div(v)/3, S, tracer
    We can not give the adiabatic losses because they are for the nonthermal
    particles, and from the RHD simulation we do not have this information.
    Therefore we just give -div(v)/3 instead of dE/dt.
    The line surface is given by the flux conservation from the initial
    surface for each line at the injection point (sf0). The surface
    corresponds to that one perpendicular to the line directon."""

    # Discard the excluded lines from the "sf0_array" and add the excluded
    # line surface to the previous line
    #print np.sum(sf0_array[:, 0])
    aux = 1
    while aux == 1:
        aux = 0
        sf0 = sf0_array[:, 0]
        for ll in range(len(sf0)):
            sf0_ind =  (sf0_array[ll,1], sf0_array[ll,2])
            if sf0_ind in excluded_lines:
                if ll == 0:
                    sf0[1] += sf0[0]
                else:
                    sf0[ll-1] += sf0[ll]
                sf0_array = np.delete(sf0_array, ll, 0)
                aux = 1
                break

    #print np.sum(sf0_array[:, 0]), len(sf0_array), len(all_lines)
    if len(sf0) != len(all_lines):
        raise RuntimeError('Surface indices do not match with line indices')

    # Discard lines with mixing
    aux = 1
    while aux == 1:
        aux = 0
        sf0 = sf0_array[:, 0]
        sf0_ind = sf0_array[:, 1]
        for ll in range(len(sf0)):
            x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*all_lines[ll])
            for tr in tracer:
                if tr > tr0:
                    aux = 1
                    break
            if aux == 1:
                if ll == 0:
                    sf0[1] += sf0[0]
                else:
                    sf0[ll-1] += sf0[ll]
                sf0_array = np.delete(sf0_array, ll, 0)
                del all_lines[ll]
                break
    #print np.sum(sf0_array[:, 0]), len(sf0_array), len(all_lines)
    #print np.sum(sf0_array[:, 0])/(1e12)**2., len(sf0_array), len(all_lines)
    #print 2010612.22971
    if len(sf0) != len(all_lines):
        raise RuntimeError('Surface indices do not match with line indices')

    # Compute secondary variables and create the output files
    for ll, line in enumerate(all_lines):
        #with open('lines/'+output_file + '_' + str(ll+1) + '.dat', 'w') as f:
        with open('lines/lines' + str(ll+1).zfill(3), 'w') as f:
            x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)
            x = np.array(x)
            y = np.array(y)
            z = np.zeros_like(x) # WARNING: modify if 3D
            dens = np.array(dens)
            eps = np.array(eps)
            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.zeros_like(vx) # WARNING: modify if 3D
            div = np.array(div)
            tracer = np.array(tracer)
            time = np.array(time)

            # We will compute secondary variables only for nonzero cells
            nz = np.ones_like(time, dtype=bool)
            for l, t in enumerate(time):
                if l != 0 and t==0:
                    nz[l] = 0
            v = np.zeros_like(time)
            v[nz] = np.sqrt(vx[nz]**2.+vy[nz]**2.)

            gamma = np.zeros_like(time)
            gamma[nz] = 1./np.sqrt(1.-(v[nz]/c)**2.)

            P = np.zeros_like(time)
            P[nz] = (gammaad-1.)*dens[nz]*eps[nz]

            h = np.zeros_like(time)
            h[nz] = 1.+eps[nz]/c**2.+P[nz]/dens[nz]/c**2.
            B = np.zeros_like(time)
            #B0 = np.sqrt(fB0*8.*np.pi*(dens[0]*h[0]*c**2.-dens[0]*c**2.-P[0]))

            # Energy flux equality
            B0 = np.sqrt(fB0*4.*np.pi*(dens[0]*h[0]*c**2.))
            B[nz] = B0*(dens[nz]*v[0]*gamma[0]/dens[0]/v[nz]/gamma[nz])**0.5

            fe = np.zeros_like(time)
            fe[nz] = dens[nz]*gamma[nz]**2.*h[nz]*v[nz]#-dens[nz]*gamma[nz]*v[nz]
            sf = np.zeros_like(time)
            sf[nz] = sf0[ll]*fe[0]/fe[nz] # Energy flux conservation
#            Particle flux conservation
#            sf[nz] = sf0[ll]*(dens[0]*v[0]*gamma[0]/dens[nz]/v[nz]/gamma[nz])

            # Save to file
            line_length = len(time)
            for n in range(line_length):
                f.write(str(x[n]) + ' ')
                f.write(str(y[n]) + ' ')
                f.write(str(z[n]) + ' ')
                f.write(str(vx[n]) + ' ')
                f.write(str(vy[n]) + ' ')
                f.write(str(vz[n]) + ' ')
                f.write(str(dens[n]) + ' ')
                f.write(str(P[n]) + ' ')
                f.write(str(B[n]) + ' ')
                f.write(str(-div[n]/3.) + ' ')
                f.write(str(sf[n]) + ' ')
                f.write(str(tracer[n]))
                f.write('\n')
    return all_lines


def initial_position(xc, yc):
    """Returns the initial position at t=0"""

    return xc, yc


def initial_indices(ic, jc):
    """Returns the initial indices at t=0"""

    return ic, jc


def initial_variables(densc, epsc, vxc, vyc, divc, tracerc):
    """Returns the initial variables at t=0"""

    return densc, epsc, vxc, vyc, divc, tracerc


def update_position(xc, yc, vxc, vyc, tstep):
    """Returns the new position after a time step"""

    xc = xc+vxc*tstep
    yc = yc+vyc*tstep
    return xc, yc


def update_indices(xc, yc, xl, yl, xr, yr, ic, jc, vxc, vyc, nx, ny):
    """Returns the cell indices for a given coordinates"""

    if vxc > 0:
        iend = nx
        di = 1
    else:
        iend = -1  # xrange function does not include the last element
        di = -1
    for i in xrange(ic, iend, di):
        if (xc >= xl[i]) and (xc < xr[i]):
            ic = i
            break

    if vyc > 0:
        jend = ny
        dj = 1
    else:
        jend = -1
        dj = -1
    for j in xrange(jc, jend, dj):
        if (yc >= yl[j]) and (yc < yr[j]):
            jc = j
            break
    return ic, jc


def interpolate(xc, yc, ic, jc, x, y, dens, eps, vx, vy, div, tracer,
                nx, ny, int_method, int_test):
    """Interpolate the physical variables at (xc, yc)"""

    dens_test = 0
    if int_test == 1:
        dens_test = dens[ic, jc]
    elif int_test == 2:
        dens_test = im.bilinear(xc, yc, ic, jc, x, y, dens, nx, ny)
    elif int_test == 3:
        dens_test = im.one_dimensional(xc, yc, ic, jc, x, y, dens, nx, ny)

    if int_method == 0:
        return dens[ic, jc], eps[ic, jc], vx[ic, jc], vy[ic, jc], div[ic, jc], tracer[ic, jc], dens_test
    elif int_method == 1:
        densi = im.bilinear(xc, yc, ic, jc, x, y, dens, nx, ny)
        epsi = im.bilinear(xc, yc, ic, jc, x, y, eps, nx, ny)
        vxi = im.bilinear(xc, yc, ic, jc, x, y, vx, nx, ny)
        vyi = im.bilinear(xc, yc, ic, jc, x, y, vy, nx, ny)
        divi = im.bilinear_div(xc, yc, ic, jc, x, y, div, nx, ny)
        traceri = im.bilinear(xc, yc, ic, jc, x, y, tracer, nx, ny)
        return densi, epsi, vxi, vyi, divi, traceri, dens_test
    elif int_method == 2:
        densi = im.one_dimensional(xc, yc, ic, jc, x, y, dens, nx, ny)
        epsi = im.one_dimensional(xc, yc, ic, jc, x, y, eps, nx, ny)
        vxi = im.one_dimensional(xc, yc, ic, jc, x, y, vx, nx, ny)
        vyi = im.one_dimensional(xc, yc, ic, jc, x, y, vy, nx, ny)
        divi = div[ic, jc]
        traceri = im.one_dimensional(xc, yc, ic, jc, x, y, tracer, nx, ny)
        return densi, epsi, vxi, vyi, divi, traceri, dens_test


def buffer_diff(densc, densc2, int_diff):
    """Keeps the difference between the two interpolation methods"""

    int_diff.append(abs(densc-densc2)/densc*100.)
    return


def tstep_test(ic, jc, ic_aux, jc_aux):
    """Checks if the time step is too large"""

    if abs(ic-ic_aux) > 1 or abs(jc-jc_aux) > 1:
        raise RuntimeError('Step larger than 1 cell. '
                           'Reduce the tstep parameter')


def code_units_to_CGS(all_lines, sf0, c, rho0, a):
    """Converts from code units to the CGS unit system."""

    sf0[:, 0] = sf0[:, 0]*a**2.
    all_lines_new = []
    for line in all_lines:
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)*a
        y = np.array(y)*a
        i = np.array(i)
        j = np.array(j)
        dens = np.array(dens)*rho0
        eps = np.array(eps)*c**2.
        vx = np.array(vx)*c
        vy = np.array(vy)*c
        div = np.array(div)*c/a
        tracer = np.array(tracer)
        time = np.array(time)*a/c

        line_new = []
        for l in range(len(x)):
            line_new.append([x[l], y[l], i[l], j[l], dens[l],
                             eps[l], vx[l], vy[l], div[l], tracer[l],
                             time[l]])
        all_lines_new.append(line_new)
    return all_lines_new, sf0


def compute_lines(x, y, xl, yl, xr, yr, vx, vy, dens, eps, tracer,
                  injec, sf0, tstep,
                  nx, ny, lx, ly, dx, dy, gammaad, div, itemax, resamp,
                  int_method, int_test, CGS_units, c, rho0, a, fB0, tr0,
                  input_file, output_file):
    """Returns the current lines for a given RHD simulation"""

    print '\nComputing the current lines...'

    start_time = time.time()
    all_lines = []
    int_diff = []
    excluded_lines = []

    for i0, j0 in injec:
        line_values = []
        tc = 0.
        xc, yc = initial_position(x[i0], y[j0])
        ic, jc = initial_indices(i0, j0)
        densc, epsc, vxc, vyc, divc, tracerc = initial_variables(dens[ic, jc],
                                                 eps[ic, jc], vx[ic, jc],
                                                 vy[ic, jc], div[ic,jc],
                                                 tracer[ic, jc])
        ic_aux, jc_aux = ic, jc
        buffer_present_line(xc, yc, ic, jc, densc, epsc,
                            vxc, vyc, divc, tracerc, tc, line_values)

        ite = 1
        while True:
            tc += tstep
            xc, yc = update_position(xc, yc, vxc, vyc, tstep)

            if (xc > lx) or (xc < 0) or (yc > ly) or (yc < 0):
                break

            ic, jc = update_indices(xc, yc, xl, yl, xr, yr, ic, jc,
                                        vxc, vyc, nx, ny)

            tstep_test(ic, jc, ic_aux, jc_aux)
            ic_aux, jc_aux = ic, jc

            try:
                (densc, epsc, vxc,
                 vyc, divc, tracerc, densc2) = interpolate(xc, yc, ic, jc, x, y,
                                                           dens, eps, vx, vy, div,
                                                           tracer, nx, ny, int_method,
                                                           int_test)
            except Exception:
                print 'WARNING: Interpolation failed for line starting at [{}, {}]'.format(i0,j0)

            if int_test != 0:
                buffer_diff(densc, densc2, int_diff)

            buffer_present_line(xc, yc, ic, jc, densc, epsc,
                                vxc, vyc, divc, tracerc, tc, line_values)
            ite += 1
            if ite == itemax:
                print 'WARNING: Line starting at [{}, {}] did not converged. Increase itemax parameter?'.format(i0,j0)
                excluded_lines += [(i0,j0)]
                break

        if ite != itemax:
            buffer_all_lines(all_lines, line_values)

    if resamp != 0:
        all_lines = resamp_line(all_lines, resamp, a, CGS_units)

    if CGS_units == 1:
        all_lines, sf0 = code_units_to_CGS(all_lines, sf0, c, rho0, a)

    save_all_lines(output_file, all_lines)

    if CGS_units == 1:
        all_lines = save_one_file_per_line(output_file, all_lines, gammaad, c, fB0,
                                           tr0, sf0, excluded_lines, input_file)
    else:
        print 'WARNING: current lines not saved. CGS unit conversion is off'

    if int_test != 0:
        print "Density difference between the two interpolation methods: " \
            "{:.2f} %".format(np.average(np.array(int_diff)))
    print "Done (elapsed time: {:.0f} seconds) "\
          .format(time.time() - start_time)
    return all_lines
