# Author: Xavier Paredes-Fortuny (xparedesfortuny@gmail.com)
# License: MIT, see LICENSE.md

import numpy as np


def pulsar(i, j, i2, j2, R_b, rho_0a, xp, yp, xs, ys, densty, eps,
           velx, vely, xznl, yznl, gamma, c, gridlx, gridly, nx, ny, theo,
           rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins):
    """Computes some physical quantities for testing purposes"""

    t_b = R_b/c

    if theo == 1:
        print '\nTHEORETICAL EXPECTED VALUES:'
        rwp = 1.
        rwi = 1.

        # Convert to CGS units
        rhowp = rhowp*rho_0a
        uwp   = uwp*c**2
        vwp   = vwp*c
        rwp   = rwp*R_b
        rhowi = rhowi*rho_0a
        uwi   = uwi*c**2.
        vwi   = vwi*c
        rwi   = rwi*R_b
        yp    = yp*R_b
        ys    = ys*R_b
        gridlx = gridlx*R_b
        gridly = gridly*R_b
    else:
        print '\nSIMULATED VALUES:'
        # Convert to CGS units
        rhowp = densty[i,j]*rho_0a
        uwp = eps[i,j]*c**2
        vwp = np.sqrt(velx[i,j]**2.+vely[i,j]**2.)*c
        rwp = np.sqrt((xznl[i]-xp)**2.+(yznl[j]-yp)**2.)*R_b
        rhowi = densty[i2,j2]*rho_0a
        uwi = eps[i2,j2]*c**2.
        vwi = np.sqrt(velx[i2,j2]**2.+vely[i2,j2]**2.)*c
        rwi = np.sqrt((xznl[i2]-xs)**2.+(yznl[j2]-ys)**2.)*R_b
        yp = yp*R_b
        ys = ys*R_b
        gridlx = gridlx*R_b
        gridly = gridly*R_b

    # Compute physical quantities
    dorb = (ys-yp)

    S_s = 4*np.pi*rwi**2.
    Mdot = rhowi*vwi*S_s

    pres_p = (gamma-1.0)*rhowp*uwp
    W_p = 1./np.sqrt(1.-(vwp/c)**2.)

    h_p = 1+uwp/c**2.+pres_p/rhowp/c**2.
    S_p = 4.*np.pi*rwp**2.

    eta = S_p*(vwp**2.*W_p**2.*rhowp*h_p+pres_p)/Mdot/vwi

    Lsd = S_p*(vwp*c**2.*W_p**2.*rhowp*h_p-W_p*rhowp*vwp*c**2.)

    Rp = eta**0.5/(1.+eta**0.5)*dorb

    print 'Adiabatic Coeficient = {:.2f}'.format(gamma)
    print 'lr = {:.1e} [cm]'.format(gridlx)
    print 'lz = {:.1e} [cm]'.format(gridly)
    print 'nr = {} [cells]'.format(nx)
    print 'nz = {} [cells]'.format(ny)

    print '\nAt d_pulsar = {:.2e} [cm] and d_star = {:.2e} [cm]:'.format(rwp, rwi)
    print 'Density of the pulsar = {:.2e} [g/cm^3]'.format(rhowp)
    print 'Density of the star/obstacle = {:.2e} [g/cm^3]'.format(rhowi)
    print 'Internal energy of the pulsar = {:.2e} [erg/g]'.format(uwp)
    print 'Internal energy of the star/obstacle = {:.2e} [erg/g]'.format(uwi)
    print 'Velocity of the pulsar wind = {:.4f} [c]'.format(vwp/c)
    print 'Velocity of the stellar wind = {:.2e} [cm/s]'.format(vwi)

    print '\nLsd = {:.1e} [erg/s]'.format(Lsd)
    print 'Mdot (Mo/yr) = {:.3e}'.format(Mdot/1.98855e33*31557600)
    print 'Mdot (g/s) = {:.3e}'.format(Mdot)
    print 'Mo/yr (g/s) = {:.3e}'.format(1.98855e33/31557600)
    print 'eta = {:.3f}'.format(eta)

    print "\nContact discontinuity measured with respect the pulsar position:"
    print 'Rp (cm) = {:.2e}'.format(Rp)
    print 'Rp (R_b) = {:.2f}'.format(Rp/R_b)

    print "\nPulsar position with respect the coordinate origin:"
    print 'yp (cm) = {:.2e}'.format(yp)
    print 'yp (R_b) = {:.2f}'.format(yp/R_b)

    print "\nStar position with respect the coordinate origin:"
    print 'ys (cm) = {:.2e}'.format(ys)
    print 'ys (R_b) = {:.2f} '.format(ys/R_b)

    print "\nDistance between the pulsar and the star:"
    print 'ys-yp (cm) = {:.2e}'.format(ys-yp)
    print 'ys-yp (R_b) = {:.2f}'.format((ys-yp)/R_b)



def jet(i, j, i2, j2, R_b, rho_0a, xp, yp, xs, ys, densty, eps,
        velx, vely, xznl, yznl, gamma, c, gridlx, gridly, nx, ny, theo,
        rhowp, uwp, vwp, rhowi, uwi, vwi, rin, rins):
    """Computes some physical quantities for testing purposes"""

    t_b = R_b/c

    if theo == 1:
        print '\n------------------------\n\nTHEORETICAL EXPECTED VALUES:'
        # Convert to CGS units
        rhowp = rhowp*rho_0a
        uwp   = uwp*c**2
        vwp   = vwp*c
        rwp   = gridlx*R_b
        rhowi = rhowi*rho_0a
        uwi   = uwi*c**2.
        vwi   = vwi*c
        rwi   = rins*R_b
        yp    = yp*R_b
        ys    = ys*R_b
        gridlx = gridlx*R_b
        gridly = gridly*R_b
        yy = 0.
    else:
        print '\n------------------------\n\nSIMULATED VALUES:'
        # Convert to CGS units
        rhowp = densty[i,j]*rho_0a
        uwp = eps[i,j]*c**2
        vwp = np.sqrt(velx[i,j]**2.+vely[i,j]**2.)*c
        rwp = gridlx*R_b
        rhowi = densty[i2,j2]*rho_0a
        uwi = eps[i2,j2]*c**2.
        vwi = np.sqrt(velx[i2,j2]**2.+vely[i2,j2]**2.)*c
        rwi = np.sqrt((xznl[i2]-xs)**2.+(yznl[j2]-ys)**2.)*R_b
        yp = yp*R_b
        ys = ys*R_b
        gridlx = gridlx*R_b
        gridly = gridly*R_b
        yy = yznl[j]*R_b

    # Compute physical quantities
    dorb  = (ys-yp)

    S_s    = 4*np.pi*rwi**2.
    Mdot   = rhowi*vwi*S_s

    pres_p = (gamma-1.0)*rhowp*uwp
    W_p    = 1./np.sqrt(1.-(vwp/c)**2.)
    h_p    = 1+uwp/c**2.+pres_p/rhowp/c**2.
    S_p    = np.pi*rwp**2.

    eta = S_p*(vwp**2.*W_p**2.*rhowp*h_p+pres_p)/Mdot/vwi

    Lsd = S_p*(rhowp*W_p**2.*c**2.*h_p*vwp-rhowp*W_p*vwp*c**2.)
    eta = Lsd/Mdot/c/vwi

    Rp = ys-np.sqrt(1/4./eta)*gridlx

    print 'Adiabatic Coeficient = {:.2f}'.format(gamma)
    print 'lr = {:.1e} [cm]'.format(gridlx)
    print 'lz = {:.1e} [cm]'.format(gridly)
    print 'nr = {} [cells]'.format(nx)
    print 'nz = {} [cells]'.format(ny)

    print '\nAt d_jet = {:.2e} [cm] and d_star/obstacle = {:.2e} [cm]'.format(yy, rwi)
    print 'Density of the jet = {:.2e} [g/cm^3]'.format(rhowp)
    print 'Density of the star/obstacle = {:.2e} [g/cm^3]'.format(rhowi)
    print 'Internal energy of the jet = {:.2e} [erg/g]'.format(uwp)
    print 'Internal energy of the star/obstacle = {:.2e} [erg/g]'.format(uwi)
    print 'Velocity of the jet wind = {:.2e} [c]'.format(vwp/c)
    print 'Velocity of the stellar wind = {:.2e} [km/s]'.format(vwi*1e-5)


    print '\nLsd = {:.1e} [erg/s]'.format(Lsd)
    print 'Mdot (Mo/yr) = {:.3e}'.format(Mdot/1.98855e33*31557600)
    print 'Mdot (g/s) = {:.3e}'.format(Mdot)
    print 'Mo/yr (g/s) = {:.3e}'.format(1.98855e33/31557600)
    print 'eta = {:.3f}'.format(eta)

    print "\nContact discontinuity measured with respect the jet position:"
    print 'Rp (cm) = {:.2e}'.format(Rp)
    print 'Rp (R_b) = {:.2f}'.format(Rp/R_b)

    print "\nJet position with respect the coordinate origin:"
    print 'yp (cm) = {:.2e}'.format(yp)
    print 'yp (R_b) = {:.2f}'.format(yp/R_b)

    print "\nStar/obstacle position with respect the coordinate origin:"
    print 'ys (cm) = {:.2e}'.format(ys)
    print 'ys (R_b) = {:.2f} '.format(ys/R_b)

    print "\nDistance between the jet and the star:"
    print 'ys-yp (cm) = {:.2e}'.format(ys-yp)
    print 'ys-yp (R_b) = {:.2f}'.format((ys-yp)/R_b)
