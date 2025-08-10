import numpy as np

def setup_coil_and_injection(t_ana, t_g, G_array, t_inj_0, t_ip, ip, nint):
    """
    コイル電流・インジェクション領域・トロイダル電流のセットアップをまとめて行う関数

    Returns
    -------
    G_round : np.ndarray
        PFコイル電流（丸め済み）
    I_inA : np.ndarray
        インジェクション領域のインデックス配列
    I_inj : np.ndarray
        インジェクション電流配列
    I_tor_def : float
        トロイダル電流（丸め済み, [A]）
    """
    # PFコイル電流
    it_g = np.where(np.isclose(t_g*1e3, t_ana, atol=1e-4))[0][0]
    G = np.array([G_array[it_g,0], np.sum(G_array[it_g,1:4]), G_array[it_g,4]])
    G[G < -100] = 0
    G_round = np.round(G, -2)
    if t_ana >= t_inj_0 and G_round[1] == 0:
        G_round[1] = -100

    # インジェクション領域
    I_inA = np.zeros(21)
    I_inA[2:4] = 1      # Fortran 3, 4 -> Python 2, 3
    I_inA[4:8] = 2      # Fortran 5, 6, 7, 8 -> Python 4, 5, 6, 7
    I_inA[8:10] = 3     # Fortran 9, 10 -> Python 8, 9

    I_inj_total = np.sum(G)

    # インジェクション電流
    I_inj = np.zeros(21)
    I_inj[1:4] = G_round[0:3]

    # トロイダル電流
    it_ip = np.where(np.isclose(t_ip*1e3, t_ana, atol=5e-4))[0][0]
    I_tor_def = -np.round(ip[it_ip]*1e3, -2)   # [A]

    return G_round, I_inA, I_inj, I_tor_def, I_inj_total