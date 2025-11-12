import os
import subprocess
import itertools
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np

def run_usalign(pdb1, pdb2, usalign_path="USalign"):
    """
    调用 USalign / TMalign 对两条 PDB 主链进行比对，
    返回 (pdb1, pdb2, RMSD, TM_1, TM_2)
    """
    try:
        cmd = [usalign_path, pdb1, pdb2]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out = result.stdout

        rmsd, tm1, tm2 = None, None, None

        for line in out.splitlines():
            if "RMSD=" in line:
                try:
                    rmsd = float(line.split("RMSD=")[1].split(",")[0].strip())
                except:
                    pass
            if "TM-score=" in line and "normalized by length of Structure_1" in line:
                try:
                    tm1 = float(line.split("TM-score=")[1].split()[0])
                except:
                    pass
            if "TM-score=" in line and "normalized by length of Structure_2" in line:
                try:
                    tm2 = float(line.split("TM-score=")[1].split()[0])
                except:
                    pass

        return pdb1, pdb2, rmsd, tm1, tm2

    except Exception as e:
        print(f"Error running USalign/TMalign on {pdb1} vs {pdb2}: {e}")
        return pdb1, pdb2, None, None, None


def batch_tmalign(pdb_folder, output_dir="results", usalign_path="USalign", num_workers=4, save_results=False):
    """
    批量运行 USalign，对 data 文件夹下所有 PDB 文件两两比对
    输出 RMSD、TM_1、TM_2 三个结果矩阵（TM-score 完整方阵）
    """
    output_dir = os.path.join(output_dir, "tmalign_results")
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
    pdb_names = [os.path.splitext(os.path.basename(f))[0] for f in pdb_files]

    if len(pdb_files) < 2:
        raise ValueError("Error: PDB 文件数量不足两条，无法比对。")

    tasks = list(itertools.combinations(pdb_files, 2))
    results = []

    print(f"共检测到 {len(pdb_files)} 个结构，将进行 {len(tasks)} 次两两比对...")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_usalign, a, b, usalign_path) for a, b in tasks]
        for f in tqdm(as_completed(futures), total=len(futures), desc="Running USalign"):
            results.append(f.result())

    df = pd.DataFrame(results, columns=["PDB1", "PDB2", "RMSD", "TM_1", "TM_2"])
    df["PDB1"] = df["PDB1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df["PDB2"] = df["PDB2"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

    # 初始化矩阵
    n = len(pdb_names)
    rmsd_matrix = pd.DataFrame(np.zeros((n, n)), index=pdb_names, columns=pdb_names)
    tm_matrix = pd.DataFrame(np.eye(n), index=pdb_names, columns=pdb_names)  # 对角线设为1

    # 填充矩阵
    for _, row in df.iterrows():
        i, j = row["PDB1"], row["PDB2"]
        rmsd_matrix.loc[i, j] = rmsd_matrix.loc[j, i] = row["RMSD"]
        tm_matrix.loc[i, j] = row["TM_1"]  # 上三角（A→B）
        tm_matrix.loc[j, i] = row["TM_2"]  # 下三角（B→A）

    # =============== 保存结果 ===============
    if save_results:
        df_path = os.path.join(output_dir, "TMalign_results.xlsx")
        matrix_path_rmsd = os.path.join(output_dir, "RMSD_matrix.xlsx")
        matrix_path_tm = os.path.join(output_dir, "TMscore_matrix.xlsx")

        df.to_excel(df_path, index=False)
        rmsd_matrix.to_excel(matrix_path_rmsd)
        tm_matrix.to_excel(matrix_path_tm)

        print(f"\n比对结果保存完成：")
        print(f"- Pairwise 数据表：{df_path}")
        print(f"- RMSD 矩阵：{matrix_path_rmsd}")
        print(f"- TM-score 矩阵（对称填充完成）：{matrix_path_tm}")
    else:
        print("未保存结果（save_results=False）")

    return df, rmsd_matrix, tm_matrix
