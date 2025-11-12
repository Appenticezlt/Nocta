import os
import time
import pandas as pd
import numpy as np
from pyrosetta import *
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm




# ======================== 单个 PDB 的能量计算 ========================
def compute_energy_for_pdb(pdb_path, relax=True):
    """
    计算单个结构的总能量和结合能
    如果有多条链：假定链1是蛋白、链2是配体
    """

    import pyrosetta
    pyrosetta.init("-ignore_unrecognized_res true -mute all -multithreading:total_threads 4")

    scorefxn = pyrosetta.rosetta.core.scoring.get_score_function()
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    start_time = time.time()
    print(f"启动计算: {pdb_name} | PID={os.getpid()}", flush=True)
    try:
        pose = pose_from_file(pdb_path)

        # ---------- Relax ----------
        if relax:
            from pyrosetta.rosetta.protocols.relax import FastRelax
            relaxer = FastRelax()
            relaxer.set_scorefxn(scorefxn)
            relaxer.apply(pose)

        # ---------- 能量计算 ----------
        total_energy = scorefxn(pose)
        num_chains = pose.num_chains()

        if num_chains >= 2:
            protein_pose = pose.split_by_chain(1)
            ligand_pose = pose.split_by_chain(2)

            E_complex = total_energy
            E_protein = scorefxn(protein_pose)
            E_ligand = scorefxn(ligand_pose)
            binding_energy = E_complex - (E_protein + E_ligand)
        else:
            E_complex = total_energy
            E_protein = total_energy
            E_ligand = 0.0
            binding_energy = 0.0

        elapsed = time.time() - start_time
        print(f" [{pdb_name}] 完成计算 | Relax={relax} | Total={total_energy:.3f} | ΔG={binding_energy:.3f} | 耗时={elapsed:.2f}s",flush=True)

        return {
            "PDB": pdb_name,
            "E_complex": E_complex,
            "E_protein": E_protein,
            "E_ligand": E_ligand,
            "Binding_Energy": binding_energy,
            "Total_Energy": total_energy,
            "Has_Ligand": num_chains >= 2,
            "Relaxed": relax,
            "Time_sec": round(elapsed, 2)
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f" [{pdb_name}] 计算出错: {e} (耗时 {elapsed:.2f}s)")
        return {"PDB": pdb_name, "Error": str(e), "Time_sec": round(elapsed, 2)}


# ======================== 批量能量计算 ========================
def batch_energy(pdb_folder, output_dir="results", num_workers=4, relax=True, save_results=False):

    
    """
    批量计算 PDB 文件的能量信息（支持多进程）
    """
    output_dir = os.path.join(output_dir, "energy_results")
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]

    if not pdb_files:
        raise ValueError("未找到任何 PDB 文件！")

    print(f"共检测到 {len(pdb_files)} 个结构，将使用 {num_workers} 个进程进行能量分析...\n")

    start_all = time.time()
    results = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(compute_energy_for_pdb, pdb_path, relax): pdb_path for pdb_path in pdb_files}
        for f in tqdm(as_completed(futures), total=len(futures), desc="Computing energies"):
            results.append(f.result())

    df = pd.DataFrame(results)
    elapsed_all = time.time() - start_all
    print(f"\n所有计算完成，总耗时 {elapsed_all / 60:.2f} 分钟")

    # =============== 保存结果 ===============
    if save_results:
        output_path = os.path.join(output_dir, "Energy_results.xlsx")
        df.to_excel(output_path, index=False)
        np.save(output_path.replace(".xlsx", ".npy"), df.to_numpy())
        print(f"能量结果已保存到: {output_path}")
    else:
        print("未保存结果（save_results=False）")

    return df
