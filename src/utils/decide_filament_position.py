import numpy as np
from src.utils.mfield_sub import cal_vecp2

def decide_filament_indices(
    t_ana, t_decay, t_inj_0, iz_puc_fit, z, z_re, bz_re, z_min, dz, nint,
    r_min, dr
):
    """
    フィラメント電流のZ, R位置インデックスと方向を決定する関数

    Returns
    -------
    ikr0 : int
        R方向のインデックス
    ikz0 : int
        Z方向のインデックス
    ikz_dir : int
        Z方向のずらし方向
    ikr : int
        R方向のインデックス（初期値）
    ikz : int
        Z方向のインデックス（初期値）
    """
    # Z position
    if t_decay <= t_ana < t_inj_0:
        ikz0 = iz_puc_fit
        ikz_dir = 1
    elif t_ana >= t_inj_0:
        ibz_re_min = np.argmin(bz_re)
        z_re_min = z_re[ibz_re_min]
        iz_puc_min = np.abs(z - z_re_min).argmin()
        ikz0 = iz_puc_min
        ikz_dir = 1
    else:
        kz0 = -1.5
        ikz0 = nint((kz0 - z_min)/dz)
        ikz_dir = 1

    # R position
    if t_ana >= t_decay:
        t1, r1 = 18.9, 0.6
        t2, r2 = 19.7, 0.3
        a = (r2 - r1) / (t2 - t1)
        b = r1 - a * t1
        kr0 = a * t_ana + b
    else:
        kr0 = 0.3

    ikr0 = nint((kr0 - r_min)/dr)
    ik = 0
    ikr, ikz = ikr0, ikz0

    return ikr0, ikz0, ikz_dir, ikr, ikz