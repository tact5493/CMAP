import numpy as np

def get_coordinate_params():
    """座標パラメータを返す"""
    ir_min, ir_max = 0, 100
    dr = 0.02
    r_min, r_max = 0.0, 2.0
    r = np.arange(r_min, r_max + dr, dr)

    iz_min, iz_max = 0, 200
    dz = 0.02
    z_min, z_max = -2.0, 2.0
    z = np.arange(z_min, z_max + dz, dz)

    dl = 0.02
    return (ir_min, ir_max, dr, r_min, r_max, r,
            iz_min, iz_max, dz, z_min, z_max, z, dl)