# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np


def read(input_file):
    """Reads a FORTRAN file as an input"""

    offset = 4
    ints = np.array(np.memmap(input_file, dtype='<i4', mode='r', shape=(1, 10),
                              order='F', offset=offset))
    offset += ints.nbytes

    ints = zip(*ints)
    ints = [i[0] for i in ints]
    (nx, ny, nsdim, nstep, igeomx, igeomy, bndmnx,
     bndmxx, bndmny, bndmxy) = ints
    if (nx >= ny):
        mn = nx+4
    else:
        mn = ny+4

    densty = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                                shape=(nx+4, ny+4), order='F', offset=offset))
    offset += densty.nbytes

    velx = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                              shape=(nx+4, ny+4), order='F', offset=offset))
    offset += velx.nbytes

    vely = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                              shape=(nx+4, ny+4), order='F', offset=offset))
    offset += vely.nbytes

    eps = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                             shape=(nx+4, ny+4), order='F', offset=offset))
    offset += eps.nbytes

    tracer = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                                shape=(nx+4, ny+4), order='F', offset=offset))
    offset += tracer.nbytes

    xznl = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                              shape=(mn), order='F', offset=offset))
    offset += xznl.nbytes

    yznl = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                              shape=(mn), order='F', offset=offset))
    offset += yznl.nbytes

    reals = np.array(np.memmap(input_file, dtype='<f8', mode='r',
                               shape=(1, 22), order='F', offset=offset))

    reals = zip(*reals)
    reals = [r[0] for r in reals]
    (gridlx, gridly, gamma, time1,
     rhowp, rhowi, rhod, uwp, uwi, vwp, vwi, vk, rc, rstar,
     xs, ys, xp, yp, rin, rins, rho_0a, R_b) = reals

    densty = densty[0:nx, 0:ny]
    velx = velx[0:nx, 0:ny]
    vely = vely[0:nx, 0:ny]
    eps = eps[0:nx, 0:ny]
    tracer = tracer[0:nx, 0:ny]
    xznl = xznl[0:nx]
    yznl = yznl[0:ny]

    return nx, ny, densty, velx, vely, eps, xznl, yznl, \
        gridlx, gridly, gamma, rho_0a, R_b, xp, yp, rin, \
        rins, tracer, time1,\
        rhowp, uwp, vwp, rhowi, uwi, vwi, xs, ys


def build_grid(xl, yl, nx, ny, lx, ly, vx, vy):
    """Computes the cell size and the center of each cell.
    The divergence is computed for each cell and considering
    their neighbors, therefore we can not compute the divergence
    while recovering the line, instead we should compute
    first the divergence for the whole grid and then
    recover it once recovering the line"""

    dx = np.empty(nx)
    for i in xrange(1, nx):
        dx[i-1] = xl[i]-xl[i-1]
    dx[nx-1] = lx-xl[nx-1]

    dy = np.empty(ny)
    for j in xrange(1, ny):
        dy[j-1] = yl[j]-yl[j-1]
    dy[ny-1] = ly-yl[ny-1]

    x = xl+dx/2.
    y = yl+dy/2.
    xr = xl+dx
    yr = yl+dy

    # WARNING: at this module the physical variables are in CODE UNITS
    v = np.zeros((nx,ny))
    gamma = np.zeros((nx,ny))
    for i in range(0, nx):
        for j in range(0, ny):
            v[i,j] = np.sqrt(vx[i,j]**2.+vy[i,j]**2.)
            gamma[i,j] = 1./np.sqrt(1.-(v[i,j])**2.)

    #print dx[0],dx[1], dx[0]/2.
    div = np.zeros((nx,ny))
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            Hx = dx[i]+dx[i+1]/2.+dx[i-1]/2.
            Hy = dy[j]+dy[j+1]/2.+dy[j-1]/2.
            div[i, j] = (gamma[i+1,j]*vx[i+1,j]-gamma[i-1,j]*vx[i-1,j])/Hx+(gamma[i,j+1]*vy[i,j+1]-gamma[i,j-1]*vy[i,j-1])/Hy
    return x, y, xr, yr, dx, dy, div


def injection(nx, ny, x, y, xp, yp, rin, input_file, nlines, xl, xr):
    """Returns an array with the injection cell indicies"""

    if input_file == 'PULSAR':
        # Locate the pulsar injection cells
        dis = np.empty((nx, ny))
        ij = []
        for i in range(nx):
            for j in range(ny):
                dis[i, j] = np.sqrt((x[i]-xp)**2.+(y[j]-yp)**2.)
                if dis[i, j] <= rin:
                    ij += [(i, j)]
        ij = np.array(ij)

        # Select only those at the border (only valid for an axsymmetric injector)
        ij_edge = []
        for i in range(ij[:, 0].max()+1):
            ij_edge += [(i, ij[ij[:, 0] == i][:, 1].min())]
            if ij[ij[:, 0] == i][:, 1].min() != ij[ij[:, 0] == i][:, 1].max():
                ij_edge += [(i, ij[ij[:, 0] == i][:, 1].max())]

        for j in ij[:, 1]:
            if (ij[ij[:, 1] == j][:, 0].max(), j) not in ij_edge:
                ij_edge += [(ij[ij[:, 1] == j][:, 0].max(), j)]


        # Sort the injector cells clockwise
        ij_edge.sort()
        up = []
        down = []
        for (ll, (i,j)) in enumerate(ij_edge):
            if y[j]-yp >= 0:
                up += [(i,j)]
            else:
                down += [(i,j)]
        down.sort(reverse=True)
        up = sorted(up, key=lambda x:x[1], reverse=True)
        down = sorted(down, key=lambda x:x[1], reverse=True)
        ij_edge_new = up+down
        if len(ij_edge) != len(ij_edge_new):
            raise RuntimeError('Could not sort the injection cells. It is needed to add the excluded line surfaces to the next line')
        ij_edge = ij_edge_new

        # Testing
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        nxmax = max(ij_edge, key=lambda x:x[0])[0]
        nymax = max(ij_edge, key=lambda x:x[1])[1]
        ax.set_xlim((0,nxmax+5))
        ax.set_ylim((0,nymax+5))
        ax.set_title('WARNING:\nThe final line numbering might change if there is\n line binning (joined current lines) and/or\n discarded lines (due to non-convergence and/or tracer cut).\nUse only for testing purposes.\n')
        for (ll, (i,j)) in enumerate(ij_edge):
            ax.annotate('{},'.format(ll+1), xy = (i, j), xytext = (0, 0), textcoords = 'offset points', size=5)
        fig.savefig('plots/injector.eps', bbox_inches="tight")

        # Line surface
        sf0 = []
        for (i,j) in ij_edge:
            sf0 += [(np.pi*(xr[i]**2.-xl[i]**2.), i, j)]
        sf0 = np.array(sf0)

    elif input_file == 'JET':
        # Locate the pulsar injection cells
        ij = []
        for i in range(nx):
            for j in range(ny):
                if y[j] <= 10:
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

            if cells_per_line%2 == 0:
                cells_per_line += 1

            while cells_per_line*nlines < nx:
                nlines += 1
            while cells_per_line*nlines > nx:
                nlines -= 1

            left_index = 0
            center_index = np.floor(cells_per_line/2.)
            right_index = 2*center_index

            index = 0
            index_aux = 0
            nlines_aux = 0
            iright = -1
            ij_new = []
            sf0 = []
            print '\nNew line resampling: '
            print ' il  ic  ir ir-(il-1) |      xl      x     xr |   x-xl   xr-x  xr-xl [a]'
            while True:
                if index == center_index:
                    if nlines_aux != nlines:
                        ileft = iright+1
                        center = ileft+index
                        iright = center+index
                        if (center < cells_per_line*3 or
                            center > cells_per_line*(nlines-3)):
                            print '{:3.0f} {:3.0f} {:3.0f}   {:3.0f}     | '\
                                '{:6.2f} {:6.2f} {:7.2f} | {:6.2f} {:6.2f} {:6.2f}'\
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
                    '{:6.2f} {:6.2f} {:7.2f} | {:6.2f} {:6.2f} {:6.2f} *NEW LINE'\
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
                ax.annotate('{},'.format(ll), xy = (i, j), xytext = (0, 0), textcoords = 'offset points', size=3)
            fig.savefig('plots/injector.eps')

            print 'Number of computed lines: {}'.format(len(ij_edge))
            if len(ij_edge) != nlines_init:
                print 'WARNING: computing {} lines instead of {} lines'\
                    .format(len(ij_edge), nlines_init)
                print '(We need an odd number of cells for each line)'

    else:
        raise RuntimeError('Injector not defined')
    return np.array(ij_edge), sf0
