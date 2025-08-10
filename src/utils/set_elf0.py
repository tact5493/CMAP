import numpy as np

# なんとなくで合わせています．今後はここを要検討
def set_elf0(t_ana, t_decay, t_inj_0):
    # Decay phase
    if t_decay <= t_ana < t_inj_0:
        return np.linspace(1, 0.3, 8)
    elif t_ana >= t_inj_0:
        # 3点を通る曲線を2次関数でFitting
        t_points = np.array([19.2, 19.4, 19.6])
        elf0_points = np.array([1, 3, 6.5])
        coeff = np.polyfit(t_points, elf0_points, 2)
        a, b, c = coeff
        y = a*t_ana**2 + b*t_ana + c
        return np.array([1, 2, 3, 4, 5, 4, 3, 2, 1])/5 * y
    # Ramp up phase
    else:
        return np.ones(8)
        #elf0 = np.array([1, 2, 3, 4, 5, 4, 3, 2, 1])/5 * y
        #elf0 = 1.   # 容器壁間を繋ぐ電流 = 1
