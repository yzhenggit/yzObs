def search_mast_cos_stis(gal_name, gal_ra, gal_dec, gal_dist_kpc,
                         instrument_name='COS/FUV',
                         search_r_kpc=100,
                         search_r_deg=0.):
    """
    Archival search for COS/FUV data within certain radius of (gal_ra, gal_dec)

    gal_ra: ra of host galaxy, in unit of degree
    gal_dec: dec of host galaxy, in unit of degree
    gal_dist_kpc: distance of host galaxy, in unit of kpc, this is only useful
                  when doing the search_r_kpc thing
    search_r_kpc: search sightlines within this radius, in unit of kpc
    search_r_deg: search sightlines within this radius, in unit of degree
                  Note that if search_r_deg is not 0, then this code
                  will always go with search_r_deg choice. Otherwise,
                  use search_r_kpc and gal_dist_Mpc options.

    History:
    Sometime in 2017 probably... YZ wrote this.
    01/26/2020, change search_mast_qso_cosfuv.py to search_mast_cos_stis.py
               because now this can handle STIS data. YZ.
    """

    import sys
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery import mast
    from yzObs.kpc2deg import kpc2deg
    from yzObs.deg2kpc import deg2kpc

    ### first check if the input instrument_name is ok:
    allow_list = ['COS/FUV', 'STIS/FUV-MAMA', 'STIS/NUV-MAMA', 'STIS/CCD']
    print("ok, input instrument_name = %s"%(instrument_name))
    if instrument_name not in allow_list:
        print("Sorry, I can only take instrument_name = COS/FUV, STIS/FUV-MAMA, STIS/NUV-MAMA, STIS/CCD.")
        sys.exit()

    gal_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg)
    gal_l = gal_coord.galactic.l.degree
    gal_b = gal_coord.galactic.b.degree
    print("*"*100)
    if search_r_deg == 0.:
        search_r_deg = kpc2deg(search_r_kpc, gal_dist_kpc)
        if search_r_deg > 30.:
            print("%.1f deg is too big, gonna just search within 5 deg instead."%(search_r_deg))
            search_r_deg = 5
            search_r_kpc = deg2kpc(search_r_deg, gal_dist_kpc)
    else:
        search_r_kpc = deg2kpc(search_r_deg, gal_dist_kpc)

    print("Searching MAST within %.2f deg/%.1f kpc of %s "%(search_r_deg, search_r_kpc, gal_name))
    print("%s at RA=%.4f, DEC=%.4f, l=%.4f, b=%.4f"%(gal_name, gal_ra, gal_dec, gal_l, gal_b))

    ##### ok, now search! #####
    mast_table = mast.Observations.query_criteria(coordinates=gal_coord,
                                                  radius=search_r_deg*u.degree,
                                                  instrument_name=instrument_name)
    if len(mast_table) == 0:
        print("Sorry, Find nothing.")
    else:
        # Check how many targets
        arx_objs = np.unique(mast_table['target_name'])
        # arx_ras = np.unique(mast_table['s_ra'])
        # arx_decs = np.unique(mast_table['s_dec'])
        for i in range(len(arx_objs)):
            print("-"*15 + " MAST/%s Search "%(instrument_name)+"-"*15)
            ind = np.where(mast_table['target_name'] == arx_objs[i])[0][0]
            arx_ra = mast_table['s_ra'][ind]
            arx_dec = mast_table['s_dec'][ind]
            arx_coord = SkyCoord(ra=arx_ra*u.deg, dec=arx_dec*u.deg, distance=gal_dist_kpc*u.kpc)

            # same as gal_coord, but with distance now
            tar_coord = SkyCoord(ra=gal_ra*u.deg, dec=gal_dec*u.deg, distance=gal_dist_kpc*u.kpc)
            arx_dist = tar_coord.separation_3d(arx_coord)

            print("Found: "+arx_objs[i])
            print('impact=%.1f kpc, ra=%.4f, dec=%.4f, l=%.4f, b=%.4f'%(arx_dist.kpc,
                                                                        arx_ra, arx_dec,
                                                                        arx_coord.galactic.l.degree,
                                                                        arx_coord.galactic.b.degree))
            print("Observation Details: ")
            print('%10s %11s %10s %06s'%("Filter", "ExpTime", "PI", "ID"))
            for j in range(len(mast_table)):
                if mast_table['target_name'][j] == arx_objs[i]:
                    print('%10s %10ds %10s %06s'%(mast_table['filters'][j],
                                                  mast_table['t_exptime'][j],
                                                  mast_table['proposal_pi'][j].split(',')[0],
                                                  mast_table['proposal_id'][j]))
            print("\n")
