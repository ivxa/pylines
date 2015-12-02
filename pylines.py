# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import os
import shutil
import numpy as np
import sys
from setup import params
from plots import plots
from compute_lines import compute_lines
import read_dat as rd
import test_rhd as test_rhd


def set_params():
    """Reads the param file and returns the parameters"""

    print '\nSetup.py: '

    input_file = params['input_file']
    print '- Input file: {}'.format(input_file)

    output_file = params['output_file']
    print '- Output filename: {}'.format(output_file + '.dat')

    testing_rhd = params['test_rhd']
    print '- Test RHD simulation: {}'.format(testing_rhd)
    if not (testing_rhd in [0, 1]):
        raise RuntimeError('Wrong input for testing_rhd parameter')

    binary_output_file = params['binary_output_file']
    print '- Binary output file: {}'.format(binary_output_file)
    if not (binary_output_file in [0, 1]):
        raise RuntimeError('Wrong input for binary_output_file parameter')

    tstep = params['tstep']
    print '- Time step fraction: {}'.format(tstep)

    itemax = params['itemax']
    print '- Maximum steps per line: {}'.format(itemax)

    resamp = params['resamp']
    print '- Resampling value: {}'.format(resamp)
    print "WARNING: with resamp = 1: disable surface test: look for two '1==1' in save_one_file_per_line() at compute_lines.py"

    nlines = params['nlines']
    print '- Number of lines: {}'.format(nlines)

    fB0 = params['fB0']
    print '- B0 internal energy fraction: {}'.format(fB0)

    tr0 = params['tr0']
    print '- Tracer cut to avoid mixing: {}'.format(tr0)

    int_method = params['int_method']
    print '- Interpolation method: {}'.format(int_method)
    if not (int_method in [0, 1, 2]):
        raise RuntimeError('Wrong input for int_method parameter')

    int_test = params['int_test']
    print '- Interpolation test: {}'.format(int_test)
    if not (int_test in [0, 1, 2, 3]):
        raise RuntimeError('Wrong input for int_test parameter')

    make_plots = params['make_plots']
    print '- Make plots: {}'.format(make_plots)
    if not (make_plots in [0, 1]):
        raise RuntimeError('Wrong input for make_plots parameter')

    plots_path = params['plots_path']
    print '- Plots path: {}'.format(plots_path)

    plot_maps = params['plot_maps']
    print '- Plot maps: {}'.format(plot_maps)
    if not (plot_maps in [0, 1]):
        raise RuntimeError('Wrong input for plot_maps parameter')

    plot_lines = params['plot_lines']
    print '- Plot lines: {}'.format(plot_lines)
    if not (plot_lines in [0, 1]):
        raise RuntimeError('Wrong input for plot_lines parameter')

    plot_profiles = params['plot_profiles']
    print '- Plot profiles: {}'.format(plot_profiles)
    if not (plot_profiles in [0, 1]):
        raise RuntimeError('Wrong input for plot_profiles parameter')

    return input_file, output_file, testing_rhd, binary_output_file, tstep,\
        itemax, resamp, nlines, fB0, tr0, int_method, int_test, make_plots,\
        plot_maps, plot_lines, plot_profiles, plots_path


def make_folder(f):
    if not os.path.exists(f):
        os.makedirs(f)
    else:
        shutil.rmtree(f)
        os.makedirs(f)


def prepare_data(input_file, nlines, testing_rhd):
    """Read the data file given by input_file. Being:
       (nx, ny): dimensions of the grid
       (x, y): array containing the coordinates of the center of the cell
       (xl, yl): array containing the coordinates of the left of the cell
       (xr, yr): array containing the coordinates of the right of the cell
       (dx, dy): array containing the cell sizes
       (vx, vy): array containing the velocity components
       dens: array containing the density
       eps: array containing the specific interal energy
       injec: array containing the (i, j) indicies of the injector cells
       (lx, ly): grid size
       a: distance unit
       rho0: density unit
       c: light velocity
       gammaad: adiabatic coeficient

       Note: to set a different input file modify this function but
       return the same variables"""

    c = 3e10

    (nx, ny, dens, vx, vy, eps, xl, yl, lx, ly, gammaad,
     tracer) = rd.read(input_file)

    if testing_rhd == 1:
        print '\nTesting the RHD simulation:'
        if input_file == 'PULSAR':
            (i_sw, j_sw) = (xs*nx/lx+2, ly*ny/ly-2)
#            (i_pw, j_pw) = (xp*nx/lx+2, (yp+rin)*ny/ly+2)
            (i_pw, j_pw) = (xp*nx/lx+2, (yp-rin)*ny/ly-2)

            test_rhd.pulsar(i_pw, j_pw, i_sw, j_sw,  a, rho0, xp, yp, xs, ys,
                            dens, eps, vx, vy, xl, yl, gammaad, c, lx, ly,
                            nx, ny, 1,
                            rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins)
            test_rhd.pulsar(i_pw, j_pw, i_sw, j_sw, a, rho0, xp, yp, xs, ys,
                            dens, eps, vx, vy, xl, yl, gammaad, c, lx, ly,
                            nx, ny, 0,
                            rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins)
        elif input_file == 'JET':
            (i_sw, j_sw) = (xs*nx/lx+2, (ys-rins)*ny/ly-2)
            (i_jw, j_jw) = (xp*nx/lx+2, (yp+10)*ny/ly+2)

            test_rhd.jet(i_jw, j_jw, i_sw, j_sw, a, rho0, xp, yp, xs, ys,
                         dens, eps, vx, vy, xl, yl, gammaad, c, lx, ly,
                         nx, ny, 1,
                         rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins)
            test_rhd.jet(i_jw, j_jw, i_sw, j_sw, a, rho0, xp, yp, xs, ys,
                         dens, eps, vx, vy, xl, yl, gammaad, c, lx, ly,
                         nx, ny, 0,
                         rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins)
        print '\nLines not computed\nSTOP\n'
        sys.exit(1)

    x, y, xr, yr, dx, dy, div = rd.build_grid(xl, yl, nx, ny, lx, ly,
                                              vx, vy, c)

    injec, sf0 = rd.injection(nx, ny, x, y, input_file,
                              nlines, xl, xr)

    return nx, ny, x, y, xl, yl, xr, yr, dx, dy, vx, vy,\
        dens, eps, injec, sf0, lx, ly, c, gammaad,\
        div, tracer


def time_step(dx, dy, vx, vy, tstep):
    """Returns the time step"""

    v = np.sqrt(vx**2.+vy**2.)

    if dx.min() < dy.min():
        ds = dx.min()
    else:
        ds = dy.min()
    tstep = tstep*ds/v.max()
    return tstep, ds


def print_info(nx, ny, lx, ly, dx, dy, tstep, c):
    """Print some information on the screen"""

    print '\nSimulation info:'
    print '- Grid dimensions: {:} x {:} cells'.format(nx, ny)
    print '- Grid size: {:.2e} x {:.2e} [cm]'.format(lx, ly)
    print '- Cell size: {:.2e} x {:.2e} [cm]'.format(dx[1], dy[1])
    print '- Time step size: {:.2f} [t0]'.format(tstep)
    print '- Equation of state pres=(gammaad-1)*dens*eps'
    return


def save_binary(output_file, all_lines):
    """Save to a .npy binary file"""
# To read the numpy binary file:
#    all_lines = np.load(output_file)
#    for line in all_lines:
#        x, y, i, j, dens, eps, vx, vy, div, time = zip(*line)

    np.save(output_file + '.npy', all_lines)
    return


def main():
    """Compute the current lines from RHD data"""

    try:
        make_folder('plots/')
        make_folder('lines/')

        (input_file, output_file, testing_rhd, binary_output_file, tstep,
         itemax, resamp, nlines, fB0, tr0, int_method, int_test, make_plots,
         plot_maps, plot_lines, plot_profiles,
         plots_path) = set_params()

        (nx, ny, x, y, xl, yl, xr, yr, dx, dy, vx, vy,
         dens, eps, injec, sf0, lx, ly, c, gammaad,
         div, tracer) = prepare_data(input_file, nlines, testing_rhd)

        tstep, ds = time_step(dx, dy, vx, vy, tstep)

        print_info(nx, ny, lx, ly, dx, dy, tstep, c)

        all_lines = compute_lines(x, y, xl, yl, xr, yr,
                                  vx, vy, dens, eps, tracer,
                                  injec, sf0, tstep,
                                  nx, ny, lx, ly, dx, dy, gammaad, div,
                                  itemax, resamp, int_method, int_test,
                                  c, fB0, tr0,
                                  input_file, output_file)

        if binary_output_file == 1:
            save_binary(output_file, all_lines)

        if make_plots == 1:
            plots(x, y, dens, eps, vx, vy, div, plots_path, all_lines,
                  plot_maps, plot_lines, plot_profiles, lx, ly,
                  c, gammaad, tracer)

        return 0

    except Exception, err:
        sys.stderr.write('\nERROR: %s\n' % str(err))
        return 1


if __name__ == '__main__':
    os.system('clear')
    print '\nRunning PyLines.py'
    os.system('rm -rf lines/current_lines_*')
    main()
