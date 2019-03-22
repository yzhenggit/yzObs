def calc_r200(mstar = 0., mhalo = 1e12):

    """
    Assuming an isothermal sphere. Based on Joo's code.
    """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np
    import sys

    if mstar != 0:
        mhalo = mstar_2_mhalo(mstar)

    # print(mhalo)

    mhalo = mhalo*u.Msun
    delta = 200 ## rho/rho_crit=200

    # calculate r200 with respect to critical density
    rho_c = Planck15.critical_density0
    r200_c = ((3*mhalo/(4.*np.pi*delta*rho_c))**(1./3.)).to(u.kpc)
    print("If use critical density, delta_c=rho/rho_c=200")
    print("r200 = %.2f kpc\n"%(r200_c.value))

    # calculate r200 with respect to matter density
    rho_m = Planck15.critical_density0 * Planck15.Om0
    r200_m = ((3*mhalo/(4.*np.pi*delta*rho_m))**(1./3.)).to(u.kpc)
    print("If use critical MATTER density, delta_c=rho/(rho_c*Omega_m)=200")
    print("r200 = %.2f kpc\n"%(r200_m.value))

    return r200_c, r200_m, mhalo.value

def mstar_2_mhalo(mstar):
    """Moster+2010's halo abundance matching. Based on Joo's code. """

    import astropy.constants as const
    from astropy.cosmology import Planck15
    import astropy.units as u
    import numpy as np

    M1 = 10**11.884 # Msun, from Table 4 or 6 in Moster+2010
    mM0 = 0.0282
    beta = 1.06
    gamma = 0.556
    Mhalo = 10**np.arange(7,13,0.001)

    # from moster equation 2
    Mstar = Mhalo * (2.*mM0*((Mhalo/M1)**(-beta)+(Mhalo/M1)**(gamma))**(-1.))

    # interpolate this to mstar
    from scipy import interpolate
    func = interpolate.interp1d(Mstar, Mhalo)
    mhalo = func(mstar)
    print("logMh = %.2f (%.1e) for logMstar = %.2f (%.1e)\n"%(np.log10(mhalo), mhalo,
                                                            np.log10(mstar), mstar))
    return mhalo
