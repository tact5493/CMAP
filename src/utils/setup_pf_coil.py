import numpy as np

def setup_pf_coil(PF_coil):
    """PFコイル電流配列をセットアップ"""
    I_pf_c = np.zeros(55)
    I_pf_c[[0, 1]] = PF_coil[[0, 1]] * 41           # PF3-1, 2
    I_pf_c[[2, 5]] = PF_coil[[2, 5]] * 36           # PF2, 6
    I_pf_c[[3, 4]] = PF_coil[[3, 4]] * 12           # PF1, 7
    I_pf_c[[6, 7]] = PF_coil[[6, 7]] * 41           # PF5-1, 2
    I_pf_c[[8, 9, 10]] = PF_coil[[8, 9, 10]] * 36   # PF4-1, 2, 3
    I_pf_c[8] /= 13
    I_pf_c[9] /= 17
    I_pf_c[10] /= 13
    k = 10
    for i in range(8, 11):
        for j in range(1, 9):
            k += 1
            I_pf_c[k] = I_pf_c[i]
            k += 1
            I_pf_c[k] = I_pf_c[i]
            if (i == 9 or i == 11) and j == 6:
                break
    return I_pf_c
