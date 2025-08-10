import numpy as np
from src.utils.mfield_sub import cal_vecp_2

def read_or_generate_A0(mfile_bin_path, ir_max, iz_max, r, z, dr, dz):
    """A_0配列を生成またはバイナリから読み込み"""
    if not mfile_bin_path.exists():
        A_0 = np.zeros((ir_max+1, iz_max+1, ir_max+1, iz_max+1))
        for iz2 in range(iz_max+1):
            for ir2 in range(ir_max+1):
                r_mid = r[ir2] + 0.5*dr
                z_mid = z[iz2] + 0.5*dz        
                for iz in range(iz_max+1):
                    for ir in range(ir_max+1):
                        A_0[ir, iz, ir2, iz2] = cal_vecp_2(1, np.array([r_mid]), np.array([z_mid]), np.array([1]), r[ir], z[iz])[3]
        A_0.tofile(mfile_bin_path)
    A_0 = np.fromfile(mfile_bin_path, dtype=np.float64).reshape(ir_max+1, iz_max+1, ir_max+1, iz_max+1)
    return A_0