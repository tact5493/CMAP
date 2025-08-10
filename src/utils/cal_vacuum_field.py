import numpy as np
from src.utils.mfield_sub import cal_vecp2

def calc_vacuum_field(ir_max, iz_max, r, z, k, r_c, z_c, I_pf_c, elect_posi, cal_vecp_2):
            """
            真空場（A_phi, Bz, mv_field, mv_field_l）を計算する
            """
            A_phi = np.zeros((ir_max+1, iz_max+1))
            mv_field, mv_field_l = np.zeros_like(A_phi), np.zeros_like(A_phi)
            Bz = np.zeros_like(A_phi)
            gmin_0, gmax_0 = 1e29, -1e29

            for iz in range(iz_max+1):
                for ir in range(ir_max+1):
                    _, sbz, _, ssbz = cal_vecp_2(k, r_c, z_c, I_pf_c, r[ir], z[iz])
                    A_phi[ir, iz] = ssbz
                    Bz[ir, iz] = sbz

                    if elect_posi(r[ir], z[iz]):
                        mv_field_l[ir, iz] = 1e29
                    else:
                        mv_field_l[ir, iz] = ssbz * 2*np.pi*r[ir]

                    if gmin_0 > mv_field[ir, iz]:
                        gmin_0 = mv_field[ir, iz]
                    if gmax_0 < mv_field[ir, iz]:
                        gmax_0 = mv_field[ir, iz]
            return A_phi, Bz, mv_field, mv_field_l, gmin_0, gmax_0
        