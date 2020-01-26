def search_mast_qso_cosfuv(gal_name, gal_ra, gal_dec, gal_dist_kpc=300,
                           within_radius_kpc=100, within_radius_deg=0.):
    """
    Archival search for COS/FUV data within certain radius of (gal_ra, gal_dec)

    gal_ra: ra of host galaxy, in unit of degree
    gal_dec: dec of host galaxy, in unit of degree
    gal_dist_kpc: distance of host galaxy, in unit of kpc, this is only useful
                  when doing the within_radius_kpc thing
    within_radius_kpc: search sightlines within this radius, in unit of kpc
    within_radius_deg: search sightlines within this radius, in unit of degree
                  Note that if within_radius_deg is not 0, then this code
                  will always go with within_radius_deg choice. Otherwise,
                  use within_radius_kpc and gal_dist_Mpc options.
    """

    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery import mast
    from yzObs.kpc2deg import kpc2deg
    from yzObs.deg2kpc import deg2kpc

    gal_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg)
    print("*"*80)
    if within_radius_deg == 0.:
        sepdeg = kpc2deg(within_radius_kpc, gal_dist_kpc)
        if sepdeg > 30.:
            print("This sep=%.1f deg is too big, gonna just search within 5 deg instead."%(sepdeg))
            sepdeg = 5
            within_radius_kpc = deg2kpc(sepdeg, gal_dist_kpc)
        print("Searching MAST within %.1f kpc of %s (RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f)"%(within_radius_kpc,
                                                                                      gal_name, gal_ra, gal_dec,
                                                                                      gal_coord.galactic.l.degree,
                                                                                      gal_coord.galactic.b.degree))
    else:
        sepdeg = within_radius_deg
        print("Searching MAST within %.1f deg of %s (RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f)"%(within_radius_deg,
                                                                                      gal_name, gal_ra, gal_dec,
                                                                                      gal_coord.galactic.l.degree,
                                                                                      gal_coord.galactic.b.degree))

    mast_table = mast.Observations.query_criteria(coordinates=gal_coord,
                                                  radius=sepdeg*u.degree,
                                                  instrument_name='COS/FUV')
    print(len(mast_table))
    if len(mast_table) == 0:
        print("Find nothing.")
    else:
        # Check how many targets
        arx_qsos = np.unique(mast_table['target_name'])
        arx_ras = np.unique(mast_table['s_ra'])
        arx_decs = np.unique(mast_table['s_dec'])
        for i in range(len(arx_qsos)):
            print("-"*15 + " MAST Archival Search "+"-"*15)
            ind = np.where(mast_table['target_name'] == arx_qsos[i])[0][0]
            arx_ra = mast_table['s_ra'][ind]
            arx_dec = mast_table['s_dec'][ind]
            if within_radius_deg == 0.:
                arx_coord = SkyCoord(ra=arx_ra*u.deg, dec=arx_dec*u.deg, distance=gal_dist_kpc*u.kpc)
                # same as gal_coord, but with distance now
                tar_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg, distance=gal_dist_kpc*u.kpc)
                arx_dist = tar_coord.separation_3d(arx_coord)

                print("Found: "+arx_qsos[i])
                print('impact=%.1f kpc, ra=%.4f, dec=%.4f, l=%.4f, b=%.4f'%(arx_dist.kpc,
                                                                        arx_ra, arx_dec,
                                                                        arx_coord.galactic.l.degree,
                                                                        arx_coord.galactic.b.degree))
            else:
                arx_coord = SkyCoord(ra=arx_ra*u.deg, dec=arx_dec*u.deg)
                tar_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg)
                arx_dist = tar_coord.separation(arx_coord)

                print("Found: "+arx_qsos[i])
                print('deg=%.1f deg, ra=%.4f, dec=%.4f, l=%.4f, b=%.4f'%(arx_dist.degree,
                                                                        arx_ra, arx_dec,
                                                                        arx_coord.galactic.l.degree,
                                                                        arx_coord.galactic.b.degree))
            print("Observation Details: ")
            print('%10s %11s %10s %06s'%("Filter", "ExpTime", "PI", "ID"))
            for j in range(len(mast_table)):
                if mast_table['target_name'][j] == arx_qsos[i]:
                    print('%10s %10ds %10s %06s'%(mast_table['filters'][j],
                                                  mast_table['t_exptime'][j],
                                                  mast_table['proposal_pi'][j].split(',')[0],
                                                  mast_table['proposal_id'][j]))
            print("\n")
