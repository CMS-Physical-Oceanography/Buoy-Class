def coriolis (lat):
    # this function inputs the lattitude in degrees
    # and calculates the coriolis parameter
    phi = np.deg2rad(lat)
    Omega_e =  2*np.pi/86400 # rad/s
    return 2*Omega_e*np.sin(phi)