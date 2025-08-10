from src.utils.mfield_sub import copy_file
from src.scripts.plot import plot_field, plot_field2, plot_psi, plot_z_Bz, plot_z_Bz2

def save_best_iteration_files(I_close_store, t_ana, path, m_in, plot_field2, plot_psi, plot_z_Bz2, z_puc, bz_puc_it, sq_error):
    i_best = int(I_close_store[-1,0])
    copy_file(f"{path}/{t_ana:.3f}/M_field_{i_best:03}.csv", f"{path}/data/M_field_{t_ana:.3f}.csv")
    copy_file(f"{path}/{t_ana:.3f}/J_field_{i_best:03}.csv", f"{path}/data/J_field_{t_ana:.3f}.csv")
    copy_file(f"{path}/{t_ana:.3f}/MV_field_flux_{i_best:03}.csv", f"{path}/data/MV_field_flux_{t_ana:.3f}.csv")
    plot_field2(m_in, t_ana, path)
    plot_psi(m_in, t_ana, path)
    plot_z_Bz2(path, t_ana, z_puc, bz_puc_it, min(sq_error[2:]))