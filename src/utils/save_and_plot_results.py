import numpy as np
from src.utils.mfield_sub import minimum_position, nint
from src.scripts.plot import plot_field, plot_field2, plot_psi, plot_z_Bz, plot_z_Bz2


def save_and_plot_results(i, t_ana, path, time_path, mf_output, jf_array, mv_flux, Bz_cal_puc, Bz_vac_puc, z_puc, bz_puc_it, sq_error, I_tor, ele_lim, dr, dz, r, A_phi, plot_field, plot_z_Bz):
    np.savetxt(f"{time_path}/M_field_{i:03}.csv", mf_output.T, delimiter = ",", fmt = "%12.4e")
    np.savetxt(f"{time_path}/J_field_{i:03}.csv", jf_array, delimiter = ",", fmt = "%12.4e")
    np.savetxt(f"{time_path}/MV_field_flux_{i:03}.csv", mv_flux, fmt = "%.3e", delimiter = ",")
    plot_z_Bz(Bz_cal_puc, Bz_vac_puc, time_path, i, z_puc, bz_puc_it, sq_error[i])
    z_Bz_data = np.array([Bz_cal_puc[:,0], Bz_cal_puc[:,1], Bz_vac_puc[:,1]])
    np.savetxt(f"{path}/data/z_Bz_{t_ana:.3f}.csv", z_Bz_data.T, fmt = "%.3e", delimiter = ",")
    ir_in, iz_in = minimum_position(A_phi)
    ir_in, iz_in = ir_in + 1, iz_in + 1
    m_in = A_phi[ir_in,iz_in] * 2*np.pi*r[ir_in]
    ir_out = nint((ele_lim[0] - r[0])/dr) + 1
    iz_out = nint((ele_lim[1] - r[0])/dz)
    m_out = A_phi[ir_out,iz_out] * 2*np.pi*ele_lim[0]
    plot_field(m_in, m_out, I_tor[-1], path, time_path, i)