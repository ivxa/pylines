# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from numpy import ma


def autocrop_img(filename):
    """Call epstools from bash to autocrop image"""
    import subprocess
    import os

    try:
        cwd, img_name = os.path.split(filename)

        bashcmd = 'epstool --copy --bbox %s %s' % (img_name, 'tmp_'+img_name)
        process = subprocess.Popen(bashcmd.split(), stdout=subprocess.PIPE, cwd=cwd)

        process.wait()
        bashcmd2 = 'mv %s %s' % ('tmp_'+img_name, img_name)
        process2 = subprocess.Popen(bashcmd2.split(), stdout=subprocess.PIPE, cwd=cwd)
    except:
        raise RuntimeError('Unable to tight layout. Increase pad_inches?')


def tracer_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines):
    """Plot the tracer map of the RHD simulation"""

    # Set-up
    from setup import params
    input_file = params['input_file']
    numbering = 1
    vel_vectors = 0
    cax_z = [0.855, 0.510, 0.03, 0.390]
    cax_vel = [0.855, 0.100, 0.03, 0.390]
    nmask = 20
    qscale = 0.1
    width0 = 0.005
    counti0 = 11
    countj0 = 16
    #

    print 'Tracer map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~{\rm [cm]}$')
        ax.set_ylabel(r'$z~{\rm [cm]}$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    im = ax.pcolormesh(xv, yv, z, cmap=cm.binary, vmin=0.01, vmax=0.99,
                       rasterized=True)

    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=0.1, pad=0.07)
    cb0 = fig.colorbar(im, cax=cax)
    if CGS_units == 0:
        cb0.set_label(r'${\rm Tracer}$', labelpad=5)
    else:
        cb0.set_label(r'${\rm Tracer}$', labelpad=5)

    if CGS_units == 0:
        minorLocator = MultipleLocator(5)
    else:
        minorLocator = MultipleLocator(0.1e15)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
    class FixedOrderFormatter(ScalarFormatter):
        """Formats axis ticks using scientific notation with a constant order of
        magnitude"""
        def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
            self._order_of_mag = order_of_mag
            ScalarFormatter.__init__(self, useOffset=useOffset,
                                     useMathText=useMathText)
        def _set_orderOfMagnitude(self, range):
            """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
            self.orderOfMagnitude = self._order_of_mag

    numbering_mod = 5
    offset_y = 0.0
    if input_file == 'JET':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(15))
        numbering_mod = 10
        offset_y = 0.01e15
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e15)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
    elif input_file == 'PULSAR':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(12))
        ax.yaxis.set_major_formatter(FixedOrderFormatter(12))
        numbering = 1
        numbering_mod = 10
        offset_y = 0.0
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e12)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)

    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'tracer_map.eps',
                bbox_inches='tight', pad_inches=0.02, dpi=300)
    autocrop_img(plots_path+'tracer_map.eps')

    if vel_vectors == 0:
#        c = ['r','b','c','y','k']
        c = ['0.1','0.25','0.4','0.55','0.7']
        c += 500*c
        for ll, line in enumerate(all_lines):
            x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

            x = np.array(x)
            y = np.array(y)
            nonzero_sel = np.ones_like(time, dtype=bool)
            for l, t in enumerate(time):
                if l != 0 and t==0:
                    nonzero_sel[l] = 0
            if numbering == 1:
                if (ll+1)%numbering_mod == 0 or ll == 0:
                    ax.annotate(str(ll+1), xy=(x[0], y[0]+offset_y), xycoords='data', size=6, color='b')
            ax.plot(x[nonzero_sel], y[nonzero_sel], lw=0.5, color=c[ll])

        fig.savefig(plots_path+'tracer_map_with_lines.eps',
                    bbox_inches='tight', pad_inches=0.07, dpi=300)
        autocrop_img(plots_path+'tracer_map_with_lines.eps')

    plt.close(fig)
    return


def pressure_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines):
    """Plot the pressure map of the RHD simulation"""

    print 'Pressure map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~[cm]$')
        ax.set_ylabel(r'$z~[cm]$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                   cmap=cm.jet, rasterized=True)

    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    cb0 = fig.colorbar(im)
    if CGS_units == 0:
        cb0.set_label(r'${\rm Pressure}~[\rho_0 c^2]$', labelpad=1)
    else:
        cb0.set_label(r'${\rm Pressure}~[erg/cm^3]$', labelpad=1)


    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'pressure_map.eps',
                bbox_inches='tight', pad_inches=0.07, dpi=300)

    for ll, line in enumerate(all_lines):
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)
        y = np.array(y)
        nonzero_sel = np.ones_like(time, dtype=bool)
        for l, t in enumerate(time):
            if l != 0 and t==0:
                nonzero_sel[l] = 0

        if (ll+1)%5 == 0 or ll == 0:
            ax.annotate(str(ll+1), xy=(x[0], y[0]), xycoords='data', size=6)
        ax.plot(x[nonzero_sel], y[nonzero_sel])

    fig.savefig(plots_path+'pressure_map_with_lines.eps',
                bbox_inches='tight', pad_inches=0.02, dpi=300)

    plt.close(fig)
    return


def beta_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines):
    """Plot the beta map of the RHD simulation"""

    print 'Beta map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~[cm]$')
        ax.set_ylabel(r'$z~[cm]$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                   cmap=cm.jet, rasterized=True)

    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    cb0 = fig.colorbar(im)
    if CGS_units == 0:
        cb0.set_label(r'${\rm \beta}$', labelpad=1)
    else:
        cb0.set_label(r'${\rm \beta}$', labelpad=1)


    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'beta_map.eps',
                bbox_inches='tight', pad_inches=0.07, dpi=300)

    for ll, line in enumerate(all_lines):
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)
        y = np.array(y)
        nonzero_sel = np.ones_like(time, dtype=bool)
        for l, t in enumerate(time):
            if l != 0 and t==0:
                nonzero_sel[l] = 0

        if (ll+1)%5 == 0 or ll == 0:
            ax.annotate(str(ll+1), xy=(x[0], y[0]), xycoords='data', size=6)
        ax.plot(x[nonzero_sel], y[nonzero_sel])

    fig.savefig(plots_path+'beta_map_with_lines.eps',
                bbox_inches='tight', pad_inches=0.07, dpi=300)

    plt.close(fig)
    return


def eps_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines):
    """Plot the internal energy map of the RHD simulation"""

    print 'Specific internal energy map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~[cm]$')
        ax.set_ylabel(r'$z~[cm]$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                   cmap=cm.jet, rasterized=True)

    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    cb0 = fig.colorbar(im)
    if CGS_units == 0:
        cb0.set_label(r'${\rm Specific~internal~energy}~[c^2]$', labelpad=1)
    else:
        cb0.set_label(r'${\rm Specific~internal~energy}~[erg/g]$', labelpad=1)


    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'eps_map.eps',
                bbox_inches='tight', pad_inches=0.07, dpi=300)

    for ll, line in enumerate(all_lines):
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)
        y = np.array(y)
        nonzero_sel = np.ones_like(time, dtype=bool)
        for l, t in enumerate(time):
            if l != 0 and t==0:
                nonzero_sel[l] = 0

        if (ll+1)%5 == 0 or ll == 0:
            ax.annotate(str(ll+1), xy=(x[0], y[0]), xycoords='data', size=6)
        ax.plot(x[nonzero_sel], y[nonzero_sel])

    fig.savefig(plots_path+'eps_map_with_lines.eps',
                bbox_inches='tight', pad_inches=0.07, dpi=300)

    plt.close(fig)
    return


def density_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines):
    """Plot the density map of the RHD simulation"""

    # Set-up
    from setup import params
    input_file = params['input_file']
    numbering = 1
    vel_vectors = 0
    cax_z = [0.855, 0.510, 0.03, 0.390]
    cax_vel = [0.855, 0.100, 0.03, 0.390]
    nmask = 20
    qscale = 0.1
    width0 = 0.005
    counti0 = 11
    countj0 = 16
    #

    print 'Density map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~{\rm [cm]}$')
        ax.set_ylabel(r'$z~{\rm [cm]}$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                   cmap=cm.jet, rasterized=True)

    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    if vel_vectors == 0:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=0.1, pad=0.07)
        cb0 = plt.colorbar(im, cax=cax)
        if CGS_units == 0:
            cb0.set_label(r'$\rho~{\rm[\rho_0]}$', labelpad=5)
        else:
            cb0.set_label(r'$\rho~{\rm[g/cm^3]}$', labelpad=5)
    else:
        from numpy import ma
        M = np.ones(vx.shape, dtype='bool')
        counti = counti0
        for i in range(len(vx[:,0])):
            if counti == nmask:
                countj = countj0
                for j in range(len(vx[0,:])):
                    if countj == nmask:
                        M[i,j] = False
                        countj = 0
                    countj += 1
                counti = 0
            counti += 1
        vxc = ma.masked_array(vx, mask=M)
        vyc = ma.masked_array(vy, mask=M)
        color_scale = np.sqrt(vxc**2.+vyc**2.)/3e10
        vxc_norm = vxc/3e10/color_scale # equal length for all vectors
        vyc_norm = vyc/3e10/color_scale
        Q = ax.quiver(xv,yv,vxc_norm,vyc_norm,color_scale,angles='xy',
                      scale_units='dots', scale=qscale)#, width=width0)

        cax = fig.add_axes(cax_z)
        cb0 = fig.colorbar(im,cax = cax)
        cax2 = fig.add_axes(cax_vel)
        cb = fig.colorbar(Q,cax=cax2,ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
        cb.set_label(r'$\beta$',labelpad=9.5)
        if CGS_units == 0:
            cb0.set_label(r'$\rho~{\rm[\rho_0]}$', labelpad=6)
        else:
            cb0.set_label(r'$\rho~{\rm[g/cm^3]}$', labelpad=6)

    if CGS_units == 0:
        minorLocator = MultipleLocator(5)
    else:
        minorLocator = MultipleLocator(0.1e15)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
    class FixedOrderFormatter(ScalarFormatter):
        """Formats axis ticks using scientific notation with a constant order of
        magnitude"""
        def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
            self._order_of_mag = order_of_mag
            ScalarFormatter.__init__(self, useOffset=useOffset,
                                     useMathText=useMathText)
        def _set_orderOfMagnitude(self, range):
            """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
            self.orderOfMagnitude = self._order_of_mag

    numbering_mod = 5
    offset_y = 0.0
    if input_file == 'JET':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(15))
        numbering_mod = 10
        offset_y = 0.01e15
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e15)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
    elif input_file == 'PULSAR':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(12))
        ax.yaxis.set_major_formatter(FixedOrderFormatter(12))
        numbering = 1
        numbering_mod = 10
        offset_y = 0.0
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e12)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)

    if x.max() >= 1.5e15:
        from matplotlib.ticker import AutoMinorLocator
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'density_map.eps',
                bbox_inches='tight', pad_inches=0.02, dpi=300)
    autocrop_img(plots_path+'density_map.eps')

    if vel_vectors == 0:
#        c = ['r','b','c','y','k']
        c = ['0.1','0.25','0.4','0.55','0.7']
        c += 500*c
        for ll, line in enumerate(all_lines):
            x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

            x = np.array(x)
            y = np.array(y)
            nonzero_sel = np.ones_like(time, dtype=bool)
            for l, t in enumerate(time):
                if l != 0 and t==0:
                    nonzero_sel[l] = 0
            if numbering == 1:
                if (ll+1)%numbering_mod == 0 or ll == 0:
                    ax.annotate(str(ll+1), xy=(x[0], y[0]+offset_y), xycoords='data', size=6, color='w')
            ax.plot(x[nonzero_sel], y[nonzero_sel], lw=0.5, color=c[ll])

        fig.savefig(plots_path+'density_map_with_lines.eps',
                    bbox_inches='tight', pad_inches=0.07, dpi=300)
        autocrop_img(plots_path+'density_map_with_lines.eps')

    plt.close(fig)
    return


def doppler_plot(x, y, z, vx, vy, plots_path, CGS_units, all_lines, fsuff):
    """Plot the Doppler Boosting map of the RHD simulation"""

    # Set-up
    from setup import params
    input_file = params['input_file']
    numbering = 1
    vel_vectors = 0
    cax_z = [0.855, 0.510, 0.03, 0.390]
    cax_vel = [0.855, 0.100, 0.03, 0.390]
    nmask = 20
    qscale = 0.1
    width0 = 0.005
    counti0 = 11
    countj0 = 16
    #

    print 'Doppler Boosting map'

    np.seterr(divide='ignore', invalid='ignore')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
    else:
        ax.set_xlabel(r'$r~{\rm [cm]}$')
        ax.set_ylabel(r'$z~{\rm [cm]}$')
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')

    from setup import params
    try:
        nick_name = params['nick_name']
    except:
        nick_name = []
        pass

    if params['input_file'] == 'JET':
        if x.max() >= 1.5e15:
            if nick_name == 'brac':
                if fsuff == '_gamma_times_beta_':
                    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=1e-1, vmax=2.5e0),
                                   cmap=cm.jet, rasterized=True)
                elif fsuff == '_factor_':
                    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=1e0, vmax=4e2),
                                   cmap=cm.jet, rasterized=True)
            else:
                if fsuff == '_gamma_times_beta_':
                    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=7e-2, vmax=3e0),
                                   cmap=cm.jet, rasterized=True)
                elif fsuff == '_factor_':
                    im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=2e0, vmax=5e2),
                                   cmap=cm.jet, rasterized=True)
        else:
            if fsuff == '_gamma_times_beta_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=1e-1, vmax=4e0),
                               cmap=cm.jet, rasterized=True)
            elif fsuff == '_factor_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=0.95e0, vmax=1.5e3),#5e3
                               cmap=cm.jet, rasterized=True)
    else:
        if nick_name == 'steady':
            if fsuff == '_gamma_times_beta_':
                #im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
            elif fsuff == '_factor_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
        elif nick_name == 'clump1':
            if fsuff == '_gamma_times_beta_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
            elif fsuff == '_factor_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
        elif nick_name == 'clump5':
            if fsuff == '_gamma_times_beta_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
            elif fsuff == '_factor_':
                im = ax.pcolor(xv, yv, z, norm=LogNorm(vmin=z[z!=0].min(), vmax=z.max()),
                               cmap=cm.jet, rasterized=True)
    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))

    if vel_vectors == 0:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=0.1, pad=0.07)
        cb0 = plt.colorbar(im, cax=cax)
        if fsuff == '_gamma_times_beta_':
            cb0.set_label(r'$\Gamma~\beta$', labelpad=5)
        elif fsuff == '_factor_':
            cb0.set_label(r'$\delta^4$', labelpad=5)
    else:
        from numpy import ma
        M = np.ones(vx.shape, dtype='bool')
        counti = counti0
        for i in range(len(vx[:,0])):
            if counti == nmask:
                countj = countj0
                for j in range(len(vx[0,:])):
                    if countj == nmask:
                        M[i,j] = False
                        countj = 0
                    countj += 1
                counti = 0
            counti += 1
        vxc = ma.masked_array(vx, mask=M)
        vyc = ma.masked_array(vy, mask=M)
        color_scale = np.sqrt(vxc**2.+vyc**2.)/3e10
        vxc_norm = vxc/3e10/color_scale # equal length for all vectors
        vyc_norm = vyc/3e10/color_scale
        Q = ax.quiver(xv,yv,vxc_norm,vyc_norm,color_scale,angles='xy',
                      scale_units='dots', scale=qscale)#, width=width0)

        cax = fig.add_axes(cax_z)
        cb0 = fig.colorbar(im,cax = cax)
        cax2 = fig.add_axes(cax_vel)
        cb = fig.colorbar(Q,cax=cax2,ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
        cb.set_label(r'$\beta$',labelpad=9.5)
        if CGS_units == 0:
            cb0.set_label(r'$\rho~{\rm[\rho_0]}$', labelpad=6)
        else:
            cb0.set_label(r'$\rho~{\rm[g/cm^3]}$', labelpad=6)

    if CGS_units == 0:
        minorLocator = MultipleLocator(5)
    else:
        minorLocator = MultipleLocator(0.1e15)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
    class FixedOrderFormatter(ScalarFormatter):
        """Formats axis ticks using scientific notation with a constant order of
        magnitude"""
        def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
            self._order_of_mag = order_of_mag
            ScalarFormatter.__init__(self, useOffset=useOffset,
                                     useMathText=useMathText)
        def _set_orderOfMagnitude(self, range):
            """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
            self.orderOfMagnitude = self._order_of_mag

    numbering_mod = 5
    offset_y = 0.0
    if input_file == 'JET':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(15))
        numbering_mod = 10
        offset_y = 0.01e15
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e15)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
    elif input_file == 'PULSAR':
        ax.xaxis.set_major_formatter(FixedOrderFormatter(12))
        ax.yaxis.set_major_formatter(FixedOrderFormatter(12))
        numbering = 1
        numbering_mod = 10
        offset_y = 0.0
        if CGS_units == 0:
            minorLocator = MultipleLocator(5)
        else:
            minorLocator = MultipleLocator(0.1e12)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_locator(minorLocator)


    if x.max() >= 1.5e15:
        from matplotlib.ticker import AutoMinorLocator
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(plots_path+'doppler_'+fsuff+'map.eps',
                bbox_inches='tight', pad_inches=0.02, dpi=300)
    autocrop_img(plots_path+'doppler_'+fsuff+'map.eps')

    if vel_vectors == 0:
#        c = ['r','b','c','y','k']
        c = ['0.1','0.25','0.4','0.55','0.7']
        c += 500*c
        for ll, line in enumerate(all_lines):
            x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

            x = np.array(x)
            y = np.array(y)
            nonzero_sel = np.ones_like(time, dtype=bool)
            for l, t in enumerate(time):
                if l != 0 and t==0:
                    nonzero_sel[l] = 0
            if numbering == 1:
                if (ll+1)%numbering_mod == 0 or ll == 0:
                    ax.annotate(str(ll+1), xy=(x[0], y[0]+offset_y), xycoords='data', size=6, color='w')
            ax.plot(x[nonzero_sel], y[nonzero_sel], lw=0.5, color=c[ll])

        fig.savefig(plots_path+'doppler_'+fsuff+'map_with_lines.eps',
                    bbox_inches='tight', pad_inches=0.07, dpi=300)
        autocrop_img(plots_path+'doppler_'+fsuff+'map_with_lines.eps')

    plt.close(fig)
    return


def lines_plot(plots_path, all_lines, CGS_units, lx, ly, a):
    """Plot the computed current lines"""

    print '\nPlotting lines...'
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    if CGS_units == 0:
        ax.set_xlabel(r'$r~[a]$')
        ax.set_ylabel(r'$z~[a]$')
        ax.set_xlim((0, lx))
        ax.set_ylim((0, ly))
    else:
        ax.set_xlabel(r'$r~[cm]$')
        ax.set_ylabel(r'$z~[cm]$')
        ax.set_xlim((0, lx*a))
        ax.set_ylim((0, ly*a))
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect='equal')
    ax2.set_xlabel(r'$i$')
    ax2.set_ylabel(r'$j$')
    for ll, line in enumerate(all_lines):
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)
        y = np.array(y)
        i = np.array(i)
        j = np.array(j)
        dens = np.array(dens)
        eps = np.array(eps)
        vx = np.array(vx)
        vy = np.array(vy)
        time = np.array(time)
        nonzero_sel = np.ones_like(time, dtype=bool)
        for l, t in enumerate(time):
            if l != 0 and t==0:
                nonzero_sel[l] = 0

        ax.plot(x[nonzero_sel], y[nonzero_sel])
        ax2.plot(i[nonzero_sel], j[nonzero_sel])
        if (ll+1)%5 == 0 or ll == 0:
            ax.annotate(str(ll+1), xy=(x[0], y[0]), xycoords='data', size=6)
            ax2.annotate(str(ll+1), xy=(i[0], j[0]), xycoords='data', size=6)

    fig.savefig(plots_path+'lines_xy.eps', bbox_inches='tight',
                pad_inches=0.07)
    plt.close(fig)

    ax2.set_xlim(left=-1)
    fig2.savefig(plots_path+'lines_ij.eps', bbox_inches='tight',
                 pad_inches=0.07)
    plt.close(fig2)

    print "Done"
    return


def profile_plot(plots_path, all_lines, CGS_units):
    """Plot the phyisical quantities a long the line"""

    print '\nPlotting profiles along the lines...'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if CGS_units == 0:
        ax.set_xlabel(r'$t~[t_0]$')
        ax.set_ylabel(r'${\rm Density}~[\rho_0]$')
    else:
        ax.set_xlabel(r'$t~[s]$')
        ax.set_ylabel(r'${\rm Density}~[g/cm^3]$')
    ax.set_yscale('log')

    for line in all_lines:
        x, y, i, j, dens, eps, vx, vy, div, tracer, time = zip(*line)

        x = np.array(x)
        y = np.array(y)
        i = np.array(i)
        j = np.array(j)
        dens = np.array(dens)
        eps = np.array(eps)
        vx = np.array(vx)
        vy = np.array(vy)
        time = np.array(time)
        nonzero_sel = np.ones_like(time, dtype=bool)
        for l, t in enumerate(time):
            if t==0:
                nonzero_sel[l] = 0

        ax.plot(time, dens)
#        ax.plot(time, div)
#    ax.set_ylim(1e-22,1e-3)
    fig.savefig(plots_path+'density_profile.eps',
                bbox_inches='tight', pad_inches=0.07)
    plt.close(fig)

    print "Done"
    return


def plots(x, y, dens, eps, vx, vy, div, plots_path, all_lines,
          plot_maps, plot_lines, plot_profiles, CGS_units, lx, ly,
          a, rho0, c, gammaad, tracer):
    """Perform the selected plots"""

    # GENERAL PLOT PARAMETERS
    fig_width_pt = 0.5*512.1496             # From Latex \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    if x.max() > y.max():
        fig_size = [fig_width, fig_height]
    else:
        fig_size = [fig_height, fig_width]
    params = {'backend': 'ps',
              'font.family': 'serif',
              'axes.labelsize': 9,
              'axes.linewidth': 0.5,
              'ytick.major.width': 0.5,
              'ytick.minor.width': 0.5,
              'font.size': 9,
              'legend.fontsize': 9,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'text.latex.preamble': [r'\usepackage{txfonts}'],
              'ps.usedistiller': 'xpdf',
              'figure.figsize': fig_size}
    plt.rcdefaults()
    plt.rcParams.update(params)

    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

    beta = np.sqrt(vx**2.+vy**2.)
    P = (gammaad-1.)*dens*eps
    if CGS_units == 1:
        x = np.array(x)*a
        y = np.array(y)*a
        dens = np.array(dens)*rho0
        eps = np.array(eps)*c**2.
        P = (gammaad-1.)*dens*eps
        beta = np.sqrt(vx**2.+vy**2.)
        vx = np.array(vx)*c
        vy = np.array(vy)*c
        div = np.array(div)*c/a
        #time = np.array(time)*a/c

    if plot_maps == 1:
        print '\nPlotting maps...'
###        pressure_plot(x, y, P, vx, vy, plots_path, CGS_units, all_lines)
###        beta_plot(x, y, beta, vx, vy, plots_path, CGS_units, all_lines)
###        eps_plot(x, y, eps, vx, vy, plots_path, CGS_units, all_lines)

        tracer_plot(x, y, tracer, vx, vy, plots_path, CGS_units, all_lines)
        density_plot(x, y, dens, vx, vy, plots_path, CGS_units, all_lines)

        b = np.sqrt(vx**2.+vy**2.)/c
        bx = vx/c
        by = vy/c
        g = 1./np.sqrt(1.-b**2.)
        D4 = (1./g/(1.-by))**4.
        doppler_plot(x, y, D4, vx, vy, plots_path, CGS_units, all_lines, '_factor_')
        doppler_plot(x, y, g*b, vx, vy, plots_path, CGS_units, all_lines, '_gamma_times_beta_')



    if plot_lines == 1:
        lines_plot(plots_path, all_lines, CGS_units, lx, ly, a)

    if plot_profiles == 1:
        profile_plot(plots_path, all_lines, CGS_units)
    return
