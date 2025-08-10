import cmap.src.scripts.get_data as g

def load_experiment_data(count):
    s = g.get_CHI_Data(count, True)
    t_ip, ip = s.get_ip() # Plasma current
    t_inj, inj = s.get_inj() # Injection current
    t_g, G_array = s.get_G() # wall-Injection current
    t_puc, bz_puc = s.get_bz() # Bz (pick up coil)
    return s, t_ip, ip, t_inj, inj, t_g, G_array, t_puc, bz_puc