import numpy as np
from src.utils.mfield_sub import fitting_bz, cal_r2

def process_pickup_coil_data(t_ana, t_puc, bz_puc, z_puc, z, r):
    """
    Pick up coilデータの整形・フィッティング・インデックス取得・r2計算までをまとめて行う関数
    """
    # t_anaに最も近いインデックスを取得
    it = np.where(np.isclose(t_puc*1e3, t_ana, atol=5e-4))[0][0]

    # Bzデータを該当時刻で抽出
    bz_puc_it = bz_puc[it, :]

    # フィッティング
    z_puc_fit, bz_puc_fit, weights, fit_params = fitting_bz(z_puc, bz_puc_it)

    # Pick up coil位置インデックス
    iz_puc = np.array([np.where(np.isclose(z, zp, atol=1e-2))[0][0] for zp in z_puc])
    iz_puc_fit = np.array([np.where(np.isclose(z, zpf, atol=1e-2))[0][0] for zpf in z_puc_fit])
    ir_puc = np.where(r == 0.22)[0][0]

    # Bz分布の再構成
    z_re = np.sort(z_puc)[::-1]
    bz_re = np.concatenate([bz_puc_it[1:6], [bz_puc_it[0]], bz_puc_it[6:]])
    z_re = z_re[1:]
    bz_re = bz_re[1:]

    z_re = np.delete(z_re, [3, 5])
    bz_re = np.delete(bz_re, [3, 5])

    # r2計算
    r2_lin, r2_quad = cal_r2(z_re, bz_re)

    return bz_puc_it, z_puc_fit, bz_puc_fit, weights, fit_params, iz_puc, iz_puc_fit, ir_puc, z_re, bz_re, r2_lin, r2_quad