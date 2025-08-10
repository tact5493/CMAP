import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import cmap.src.scripts.get_data as g
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"

# -- plot --
def plot_field(m_in, m_out, I_tor, path, time_path, num):
    """
    Plot current vector and mag. field on background pic.

    Parameter
    ---
    m_in  (np.float) : Inner electode line
    m_out (np.float) : Last closed flux surface line
    path  (str)      : path of background pic., M_field, J_field and electrode
    num   (int)      : M_field and J_field file name
    """
    
    # Qvessel.pngを背景画像として読み込む
    qvessel_img = Image.open(f"{DATA_DIR}/Qvessel.png")

    # get data from file
    M_field = np.loadtxt(f"{time_path}/M_field_{num:03}.csv", delimiter = ",")
    J_field = np.loadtxt(f"{time_path}/J_field_{num:03}.csv", delimiter = ",")
    electrode = np.loadtxt(f"{path}/electrode.csv", delimiter = ",")

    # plot setting
    fig, ax = plt.subplots(figsize = (5, 10))
    plt.rcParams["font.family"] = "Arial"

    ax.tick_params(labelsize = 14)
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", which = "major", direction = "in")
    ax.tick_params(axis = "y", which = "major", direction = "in")

    # install Qvessel.png
    ax.imshow(qvessel_img, extent = [0, 2, -2, 2], aspect = 'auto')

    # plot vector (J-field)
    J_f_non_zero = J_field[(J_field[:, 2] != 0) | (J_field[:, 3] != 0)]
    ax.quiver(J_f_non_zero[:, 0], J_f_non_zero[:, 1], J_f_non_zero[:, 2], J_f_non_zero[:, 3], 
              color = 'red', scale=1, scale_units='xy', width=0.005)

    # plot mag. contour
    Z = M_field.reshape(201, 101)
    x = np.linspace(0, 201, Z.shape[1])
    y = np.linspace(-201, 201, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    ax.contour(X*1e-2, Y*1e-2, Z,
                        levels = sorted([m_in, m_out]), 
                        colors = "blue", alpha = 0.5, linewidths = 2)

    # plot electrode
    ax.plot(electrode[:, 0], electrode[:, 1], color = 'steelblue', linewidth = 5)

    # plot Pick up coil position
    z_puc = np.array([0 if i == 0 else 687e-3 - ((i-1) * 150e-3) for i in range(12)])
    r_puc = np.full_like(z_puc, 0.215)
    ax.scatter(r_puc, z_puc, marker = ",", color = "black")

    # axis setting
    ax.set_title(f"No. = {num:03}, " + r"$I_\mathrm{tor \text{-} close} = $" + f"{I_tor*1e-3:.1f} [kA]",
                 fontsize = 16)
    ax.set_xlabel(r"$\mathrm{R \, [m]}$", fontsize = 15)
    ax.set_ylabel(r"$\mathrm{Z \, [m]}$", fontsize = 15)

    # save
    plt.tight_layout()
    fig.savefig(f"{time_path}/img/field_cont_{num:03}.png",
                dpi = 300,
                bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()



def plot_field2(m_in, t_ana, path):
    """
    Plot current vector and mag. field on background pic.

    Parameter
    ---
    m_in  (np.float) : Inner electode line
    m_out (np.float) : Last closed flux surface line
    path  (str)      : path of background pic., M_field, J_field and electrode
    num   (int)      : M_field and J_field file name
    """
    
    # Qvessel.pngを背景画像として読み込む
    qvessel_img = Image.open(f"{DATA_DIR}/Qvessel.png")

    # get data from file
    M_field = np.loadtxt(f"{path}/data/M_field_{t_ana:.3f}.csv", delimiter = ",")
    J_field = np.loadtxt(f"{path}/data/J_field_{t_ana:.3f}.csv", delimiter = ",")
    electrode = np.loadtxt(f"{path}/electrode.csv", delimiter = ",")

    # plot setting
    fig, ax = plt.subplots(figsize = (5, 10))
    plt.rcParams["font.family"] = "Arial"

    ax.tick_params(labelsize = 14)
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", which = "major", direction = "in")
    ax.tick_params(axis = "y", which = "major", direction = "in")

    # install Qvessel.png
    ax.imshow(qvessel_img, extent = [0, 2, -2, 2], aspect = 'auto')

    # plot vector (J-field)
    J_f_non_zero = J_field[(J_field[:, 2] != 0) | (J_field[:, 3] != 0)]
    ax.quiver(J_f_non_zero[:, 0], J_f_non_zero[:, 1], J_f_non_zero[:, 2], J_f_non_zero[:, 3], 
              color = 'red', scale=1, scale_units='xy', width=0.005)

    # plot mag. contour
    Z = M_field.reshape(201, 101)
    x = np.linspace(0, 201, Z.shape[1])
    y = np.linspace(-201, 201, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    ax.contour(X*1e-2, Y*1e-2, Z,
                levels = np.linspace(np.min(Z)/8, m_in, 25),
                colors = "blue", alpha = 0.5, linewidths = 2)

    # plot electrode
    ax.plot(electrode[:, 0], electrode[:, 1], color = 'steelblue', linewidth = 5)

    # plot Pick up coil position
    z_puc = np.array([0 if i == 0 else 687e-3 - ((i-1) * 150e-3) for i in range(12)])
    r_puc = np.full_like(z_puc, 0.215)
    ax.scatter(r_puc, z_puc, marker = ",", color = "black")

    # axis setting
    ax.set_title(f"t = {t_ana:.3f} ms", fontsize = 16)
    ax.set_xlabel(r"$\mathrm{R \, [m]}$", fontsize = 15)
    ax.set_ylabel(r"$\mathrm{Z \, [m]}$", fontsize = 15)

    # save
    plt.tight_layout()
    fig.savefig(f"{path}/img/field_cont_{t_ana:.3f}.png",
                dpi = 300,
                bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()



def plot_z_Bz(mv_flux, vac_flux, path, num, z_puc, bz_puc, sq_error):
    """
    Plot Z-Bz of flux loop in center stack position

    Parameters
    ---
    mv_flux (np.array) : (118, 5)
        118 -> 0:24 : flux loop on center stack
        5 -> 1      : z position [m]
            0 : -1.150 [m],   1 : -1.050 [m]
            ...
            22 : 1.050 [m],  23 : -1.150 [m]
        5 -> 4      : Bz [T]
    path (str)         : Save folder
    num  (int)         : itteration number
    """
    # get data
    z = mv_flux[:,0]
    bz = mv_flux[:,1]
    bz_vac = vac_flux[:,1]

    # plot figure
    fig, ax = plt.subplots(1, 1, figsize = (6, 8))    

    fig.subplots_adjust(wspace = 0.2, hspace = 0.16)
    plt.rcParams["font.family"] = "Arial"

    ax.tick_params(labelsize = 14)
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", which = "major", direction = "in")
    ax.tick_params(axis = "y", which = "major", direction = "in")

    ax.scatter((bz-bz_vac)*1e3, z, lw = 2, marker = "o", 
            label = "Result", c = "steelblue")
    ax.scatter(bz_puc*1e3, z_puc, lw = 2, marker = "x", 
               label = "Pick up coil", c = "orangered")

    ax.legend(fontsize = 15,
                framealpha = 0.0,
                facecolor = "white",
                markerscale = 2,
                handlelength = 1)

    ax.set_ylim(np.min(z)*1.1, np.max(z)*1.1)

    ax.set_title(f"No. {num}, sq_error = {sq_error:.1e}", fontsize = 16)
    ax.set_xlabel(r"$B_\mathrm{z} \, \mathrm{[mT]}$", fontsize = 15)
    ax.set_ylabel(r"$\mathrm{Z} \, \mathrm{[m]}$", fontsize = 15)

    fig.savefig(f"{path}/img/Bz-Z_{num}.png",
                dpi = 300, bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()



def plot_z_Bz2(path, t_ana, z_puc, bz_puc, sq_error):
    """
    Plot Z-Bz of flux loop in center stack position

    Parameters
    ---
    mv_flux (np.array) : (118, 5)
        118 -> 0:24 : flux loop on center stack
        5 -> 1      : z position [m]
            0 : -1.150 [m],   1 : -1.050 [m]
            ...
            22 : 1.050 [m],  23 : -1.150 [m]
        5 -> 4      : Bz [T]
    path (str)         : Save folder
    num  (int)         : itteration number
    """
    # get data
    data = np.loadtxt(f"{path}/data/z_Bz_{t_ana:.3f}.csv", delimiter = ",")

    z = data[:,0]
    bz = data[:,1]
    bz_vac = data[:,2]

    # plot figure
    fig, ax = plt.subplots(1, 1, figsize = (6, 8))    

    fig.subplots_adjust(wspace = 0.2, hspace = 0.16)
    plt.rcParams["font.family"] = "Arial"

    ax.tick_params(labelsize = 14)
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", which = "major", direction = "in")
    ax.tick_params(axis = "y", which = "major", direction = "in")

    ax.scatter((bz-bz_vac)*1e3, z, lw = 2, marker = "o", 
            label = "Result", c = "steelblue")
    ax.scatter(bz_puc*1e3, z_puc, lw = 2, marker = "x", 
               label = "Pick up coil", c = "orangered")

    ax.legend(fontsize = 15,
                framealpha = 0.0,
                facecolor = "white",
                markerscale = 2,
                handlelength = 1)

    ax.set_ylim(np.min(z)*1.1, np.max(z)*1.1)

    ax.set_title(f"t = {t_ana:.3f}, sq_error = {sq_error:.1e}", fontsize = 16)
    ax.set_xlabel(r"$B_\mathrm{z} \, \mathrm{[mT]}$", fontsize = 15)
    ax.set_ylabel(r"$\mathrm{Z} \, \mathrm{[m]}$", fontsize = 15)

    fig.savefig(f"{path}/img/Bz-Z_{t_ana:.3f}.png",
                dpi = 300, bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()



def plot_t_Bz(t_puc, bz_puc, t_ip, ip, it, count, path):
    """
    Plot time - Ip, Iinj, Bz

    Parameters
    ---
    t_puc  (np.array)       : [s]
    Bz_puc (np.array)       : [T]
        [0]  -> Z =    0 mm,     [1]  -> Z = +687 mm
        [2]  -> Z = +537 mm,     [3]  -> Z = +387 mm
        [4]  -> Z = +237 mm,     [5]  -> Z = + 87 mm
        [6]  -> Z = - 63 mm,     [7]  -> Z = -213 mm
        [8]  -> Z = -363 mm,     [9]  -> Z = -513 mm
        [10] -> Z = -663 mm,     [11] -> Z = -813 mm
        ref. TAKEDA, master thesis. P18
    t      (np.array)      : [s] 1MS/s
    Ip     (np.array)      : [kA]
    it     (int)           : Analysis time index 
    count  (int)           : Shot number
    path   (str)           : Save folder path
    """

    # Get data
    s = g.get_CHI_Data(count, True)
    #Inj
    t_inj, inj = s.get_inj()


    fig, ax = plt.subplots(2,1,figsize = (6,4),sharex = True)
    fig.subplots_adjust(wspace = 0.2, hspace = 0.16)
    plt.rcParams["font.family"] = "Arial"

    for a in ax:
        a.tick_params(labelsize = 14)
        a.xaxis.set_ticks_position("both")
        a.tick_params(axis = "x", which = "major", direction = "in")
        a.axvline(t_puc[it]*1e3, ls = ":", c = "black")
        a.set_xlim(18.65, 20)

    ax[0].plot(t_ip*1e3, -ip, label = r"$I_{\mathrm{p}}$", lw = 2)
    ax[0].plot(t_inj*1e3, -inj, label = r"$I_{\mathrm{inj}}$", lw = 2)
    ax[0].axhline(0, c = "black", alpha = 0.5)
    ax[0].invert_yaxis()

    [ax[1].plot(t_puc*1e3, bz_puc[:,i]*1e3) for i in range(12)]
    
    ax[-1].set_xlabel(r"$\mathrm{Time \, [ms]}$", fontsize = 15)
    ax[0].set_ylabel(r"$I \, \mathrm{[kA]}$", fontsize = 15)
    ax[1].set_ylabel(r"$B_\mathrm{z} \, \mathrm{[mT]}$", fontsize = 15)
    ax[0].set_title(f"#{count}, t_ana = {t_puc[it]*1e3:.3f} ms",
                    fontsize = 16)

    fig.savefig(f"{path}/img/{count}_Bz.png",
                dpi = 300, bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()



def plot_psi(m_in, t_ana, path):
    # Qvessel.pngを背景画像として読み込む
    qvessel_img = Image.open(f"{DATA_DIR}/Qvessel.png")

    # get data from file
    M_field = np.loadtxt(f"{path}/data/M_field_{t_ana:.3f}.csv", delimiter = ",")
    electrode = np.loadtxt(f"{DATA_DIR}/electrode.csv", delimiter = ",")

    # plot setting
    fig, ax = plt.subplots(figsize = (5, 10))
    plt.rcParams["font.family"] = "Arial"

    ax.tick_params(labelsize = 14)
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", which = "major", direction = "in")
    ax.tick_params(axis = "y", which = "major", direction = "in")

    # install Qvessel.png
    ax.imshow(qvessel_img, extent = [0, 2, -2, 2], aspect = 'auto')

    # plot mag. contour
    Z = M_field.reshape(201, 101)
    x = np.linspace(0, 201, Z.shape[1])
    y = np.linspace(-201, 201, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    ax.contourf(X*1e-2, Y*1e-2, Z,
                cmap = "rainbow", alpha = 0.5)

    #. plot mag. line
    ax.contour(X*1e-2, Y*1e-2, Z,
                levels = np.linspace(np.min(Z)/8, m_in, 25),
                colors = "blue", alpha = 0.5, linewidths = 1.5)

    # plot electrode
    ax.plot(electrode[:, 0], electrode[:, 1], color = 'black', linewidth = 5)

    # plot Pick up coil position
    z_puc = np.array([0 if i == 0 else 687e-3 - ((i-1) * 150e-3) for i in range(12)])
    r_puc = np.full_like(z_puc, 0.215)
    ax.scatter(r_puc, z_puc, marker = ",", color = "black")

    # axis setting
    ax.set_title(f"t = {t_ana:.3f} ms", fontsize = 16)
    ax.set_xlabel(r"$\mathrm{R \, [m]}$", fontsize = 15)
    ax.set_ylabel(r"$\mathrm{Z \, [m]}$", fontsize = 15)

    # save
    plt.tight_layout()
    fig.savefig(f"{path}/img/psi_cont_{t_ana:.3f}.png",
                dpi = 300,
                bbox_inches = "tight", pad_inches = 0.1)
    plt.clf()
    plt.close()
