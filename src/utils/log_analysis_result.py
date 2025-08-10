def log_analysis_result(
    filepath, count, t_ana, I_pf_c, I_tor_def, G_round, elf0, t_decay, z_puc_fit, weights
):
    """
    解析結果をログファイルに出力する関数
    """
    with open(filepath, "w") as f:
        f.write(f"#{count}, t_ana = {t_ana:.3f} ms\n")
        f.write("PF3-1, 3-2, 2, 1, 7, 6, 5-2, 5-1, 4-1, 4-2, 4-3\n")
        for i, val in enumerate(I_pf_c[:11]):
            if i == len(I_pf_c[:11]) - 1:
                f.write(f"{val}\n")
            else:
                f.write(f"{val}, ")
        f.write(f"I_tor_def [A], {I_tor_def}\n")
        f.write(f"G [A], {G_round[0]}, {G_round[1]}, {G_round[2]}\n")
        f.write(f"elf0 ,{elf0}\n")
        f.write(f"t_decay [ms],{t_decay:.2f}\n")
        f.write("z [m], ")
        for zf in z_puc_fit:
            f.write(f"{zf:.2f}, ")
        f.write("\nbz weight, ")
        for w in weights:
            f.write(f"{w}, ")