def kpc2deg(impact_kpc, distance_kpc):
    """ 
    Calculation the angular size of certain length (impact) at distance (dist)
    
    dist and impact both in kpc 
    return value in degree 
    """
    # impact and dist both in kpc 
    sepdeg = impact_kpc/(distance_kpc/206265)/3600
    return sepdeg 
