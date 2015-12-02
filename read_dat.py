# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np


def read(input_file):
    """Reads a .dat file as an input"""

    c = 3e10
    Mdot = 3e-7*1.98855e33/31557600
    Rorb = 4e12
    v_w = 2e8
    norm = Mdot/4./np.pi/Rorb**2./v_w

    nx = 150-1 # WARNING: -1 only needed if the input file does not contain the ghost cells
    ny = 300-1 # WARNING: -1 only needed if the input file does not contain the ghost cells

    a = np.loadtxt(input_file, dtype=float, comments='#')
    a = a.reshape(len(a)*3/6, 6)

    vx = a[:, 0].reshape(nx, ny)
    vy = a[:, 1].reshape(nx, ny)

    dens = a[:, 3].reshape(nx, ny)*norm
    pres = a[:, 4].reshape(nx, ny)*norm
    W = 1./np.sqrt(1.-(vx**2.+vy**2.)/c**2.)

    tracer = a[:, 5].reshape(nx, ny)
    tracer = tracer/dens/W*norm
    # print tracer.shape
    # print tracer.max()
    # ii =  np.unravel_index(tracer.argmax(), tracer.shape)
    # print '{:.2e}'.format(ii[0]*4.5e11/(nx+1))
    # print '{:.2e}'.format(ii[1]*9e11/(ny+1))
    # print tracer.min()
    #assert tracer.max() <= 1, 'Tracer bigger than 1' #### asseert posar 1
    assert tracer.max() >= 0, 'Tracer smaller than 0'

    gammaad = 4./3.
    eps = pres/(gammaad-1)/dens

    lx = 4.5e11
    ly = 9e11

    dx = lx/(nx+1) # WARNING: +1 only needed if lx is defined accounting for the ghost cells
    dy = ly/(ny+1) # WARNING: +1 only needed if lx is defined accounting for the ghost cells

    xl = np.zeros(nx)
    for i in range(1, nx):
        xl[i] = xl[i-1]+dx

    yl = np.zeros(ny)
    for i in range(1, ny):
        yl[i] = yl[i-1]+dy

    return nx, ny, dens, vx, vy, eps, xl, yl, lx, ly, gammaad, tracer


def build_grid(xl, yl, nx, ny, lx, ly, vx, vy, c):
    """Computes the cell size and the center of each cell.
    The divergence is computed for each cell and considering
    their neighbors, therefore we can not compute the divergence
    while recovering the line, instead we should compute
    first the divergence for the whole grid and then
    recover it once recovering the line"""

    dx = lx/(nx+1) # WARNING: +1 only needed if lx is defined accounting for the ghost cells
    dy = ly/(ny+1) # WARNING: +1 only needed if lx is defined accounting for the ghost cells

    x = xl+dx/2.
    y = yl+dy/2.
    xr = xl+dx
    yr = yl+dy

    v = np.zeros((nx, ny))
    gamma = np.zeros((nx, ny))
    for i in range(0, nx):
        for j in range(0, ny):
            v[i, j] = np.sqrt(vx[i, j]**2.+vy[i, j]**2.)
            gamma[i, j] = 1./np.sqrt(1.-(v[i, j]/c)**2.)

    div = np.zeros((nx, ny))
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            div[i, j] = ((gamma[i+1, j]*vx[i+1, j]-gamma[i-1, j]*vx[i-1, j])/2./dx +
                         (gamma[i, j+1]*vy[i, j+1]-gamma[i, j-1]*vy[i, j-1])/2./dy)

    dx = np.ones(nx)*dx
    dy = np.ones(ny)*dy
    return x, y, xr, yr, dx, dy, div


def injection(nx, ny, x, y, input_file, nlines, xl, xr):
    """Returns an array with the injection cell indicies"""

    # Locate the pulsar jet cells
    ij = []
    for i in range(nx):
        for j in range(ny):
            if y[j] == y[0]:
                ij += [(i, j)]
    ij = np.array(ij)

    # Select only those at the border (only valid for a jet at the bottom)
    ij_edge = []
    for i in range(ij[:, 0].max()+1):
        jmax = ij[ij[:, 0] == i][:, 1].max()
        ij_edge += [(i, jmax)]

    # Resample the number of lines
    if nlines != -1:
        nlines_init = nlines

        cells_per_line = np.floor(nx/np.float(nlines))

        if cells_per_line % 2 == 0:
            cells_per_line += 1

        while cells_per_line*nlines < nx:
            nlines += 1
        while cells_per_line*nlines > nx:
            nlines -= 1

#        left_index = 0
        center_index = np.floor(cells_per_line/2.)
#        right_index = 2*center_index

        index = 0
        index_aux = 0
        nlines_aux = 0
        iright = -1
        ij_new = []
        sf0 = []
        print '\nNew line resampling: '
        print ' il  ic  ir ir-(il-1) |      xl      x     xr [cm] |   x-xl   xr-x  xr-xl  [cm]'
        while True:
            if index == center_index:
                if nlines_aux != nlines:
                    ileft = iright+1
                    center = ileft+index
                    iright = center+index
                    if (center < cells_per_line*3 or
                        center > cells_per_line*(nlines-3)):
                        print '{:3.0f} {:3.0f} {:3.0f}   {:3.0f}     | '\
                            '{:6.2e} {:6.2e} {:7.2e} | {:6.2e} {:6.2e} {:6.2e}'\
                                .format(ileft, center, iright, iright-ileft+1,
                                        xl[ileft], x[center], xr[iright],
                                        x[center]-xl[ileft], xr[iright]-x[center],
                                        xr[iright]-xl[ileft])
                    elif center > cells_per_line*3 and center < cells_per_line*4:
                        print '...'
                    ij_new += [(center, ij_edge[0][1])]
                    sf0 += [(np.pi*(xr[iright]**2.-xl[ileft]**2.), center, jmax)]
                    nlines_aux += 1
                    index = 0
                    index_aux = 0
                else:
                    break
            else:
                index += 1

        # Add the remaining cells to the last line
        cells_without_line = nx-1-iright
        if cells_without_line > 0.6*cells_per_line:
            ileft = iright+1
            center = ileft+np.int(np.floor(cells_without_line/2.))
            iright = nx-1
            print '{:3.0f} {:3.0f} {:3.0f}   {:3.0f}     | '\
                '{:6.2e} {:6.2e} {:7.2e} | {:6.2e} {:6.2e} {:6.2e} *NEW LINE'\
                .format(ileft, center, iright, iright-ileft+1,
                        xl[ileft], x[center], xr[iright],
                        x[center]-xl[ileft], xr[iright]-x[center],
                        xr[iright]-xl[ileft])
            ij_new += [(center, ij_edge[0][1])]
            sf0 += [(np.pi*(xr[iright]**2.-xl[ileft]**2.), center, jmax)]
        elif cells_without_line > 0:
            # New boundary
            iright = nx-1
            # Recenter the last line
            center += np.int(np.floor(cells_without_line/2.))
            ij_new[-1] = (center, ij_edge[0][1])
            sf0[-1] = (np.pi*(xr[iright]**2.-xl[ileft]**2.), center, jmax)

            print '{:3.0f} {:3.0f} {:3.0f}   {:3.0f}     | '\
                '{:6.2f} {:6.2f} {:7.2f} | {:6.2f} {:6.2f} {:6.2f} *CORREC.'\
                .format(ileft, center, iright, iright-ileft+1,
                        xl[ileft], x[center], xr[iright],
                        x[center]-xl[ileft], xr[iright]-x[center],
                        xr[iright]-xl[ileft])

        cells_without_line_new = nx-1-iright
        print 'Number of cells not included in any line: {}'\
            .format(cells_without_line_new)
        print 'Number of cells added to the last line: {}'\
            .format(cells_without_line)

        ij_edge = ij_new
        sf0 = np.array(sf0)

        # Testing surface
        surf_num = np.sum(sf0[:,0])
        surf_teo = np.pi*(800)**2.-np.pi*(xl[0])**2.  # the first cell starts at x[0] to avoid the singularity
        res = np.abs(surf_num-surf_teo)/surf_teo*100.
        # print surf_num
        # print surf_teo
        # print res

        # Testing
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        nxmax = max(ij_edge, key=lambda x:x[0])[0]
        nymax = max(ij_edge, key=lambda x:x[1])[1]
        ax.set_xlim((0,nxmax+5))
        ax.set_ylim((0,nymax+5))
        if nymax < 10:
            ax.set_ylim((0,30))
        for (ll, (i,j)) in enumerate(ij_edge):
            ax.annotate('{},'.format(ll), xy = (i, j), xytext = (0, 0), textcoords = 'offset points', size=5)
        fig.savefig('plots/injector.eps')

        print 'Number of computed lines: {}'.format(len(ij_edge))
        if len(ij_edge) != nlines_init:
            print 'WARNING: computing {} lines instead of {} lines'\
                .format(len(ij_edge), nlines_init)
            print '(We need an odd number of cells for each line)'

    else:
        raise RuntimeError('Injector not defined')

    return np.array(ij_edge), sf0
