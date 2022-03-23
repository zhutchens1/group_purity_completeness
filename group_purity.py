import numpy as np

def get_metrics(groupid, haloid, galproperty):
    """
    For a group catalog constructed from a mock galaxy catalog,
    compute galaxy-wise purity and completeness metrics using
    true halo IDs for comparison.

    Parameters
    ---------------------------
    groupid : iterable
        Group ID numbers after applying group-finding algorithm, length = # galaxies.
    haloid : iterable
        Halo ID numbers extracted from mock catalog halos, length = # galaxies = len(groupid).
    galproperty : iterable
        Group property by which to determine the central galaxy in the group. If all values are
        >-15 and <-27, then galproperty is assumed to be a magnitude, and the central will be the
        brightest galaxy. If all values are >0, this value is assumed to be mass, and the central
        will be selected by the maximum.

    Returns
    ---------------------------
    Suppose we map groups to halos, and define
        - N_g as the number of galaxies in the group
        - N_h as the number of galaxies in the corresponding true halo
        - N_s as the number of galaxies in the group that are correctly classified as
            members of the corresponding halo
        - N_i as the number of interlopers in the group; galaxies classified to the group
            but that do not belong to the true halo.
    
    purity : np.array
        At index `i`, purity of the group to which galaxy `i` belongs (duplicated for every
        group member). Purity is defined as the percentage of the number of galaxies in the group
        that are correctly identified as part of the halo, N_s/N_g. Because N_g = N_s + N_i, we
        the contamination fraction is given by 1 - purity = (N_g - N_s)/N_g = (N_i)/N_g. 
    completeness : np.array
        At index `i`, completeness of the group to which galaxy `i` belongs (duplicated for every
        group member). Completeness is defined as the percentage of galaxies in the halo that are 
        correctly identified as part of the group, N_s/N_h. 
    """
    groupid = np.array(groupid)
    haloid = np.array(haloid)
    return None
