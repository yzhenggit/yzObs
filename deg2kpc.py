def deg2kpc(theta_deg, distance_kpc):
    """ 
    Calculation the angular size of certain length (impact) at distance (dist)
    
    dist and impact both in kpc 
    return value in degree 
    """
    # impact and dist both in kpc 
    impact_kpc = theta_deg * (distance_kpc/206265) * 3600
    print('%.4f deg at %.1f kpc is'%(theta_deg, distance_kpc))
    print('%.2f kpc'%(impact_kpc))
    return impact_kpc
