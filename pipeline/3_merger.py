import pymongo
import time
import sys
import os
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors
from tqdm import tqdm

# --- 1. 环境配置与 SA打分器加载 ---
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# 尝试加载 RDKit 的 SAscore (可合成性评价)
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
try:
    import sascorer
except ImportError:
    logging.error("❌ 无法加载 sascorer，请确保 RDKit 贡献库路径正确。")
    sascorer = None

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MolNexusRefiner:
    def __init__(self):
        self.client = pymongo.MongoClient(config.MONGO_URI)
        self.db_raw = self.client["drug_raw_data"]
        self.db_refined = self.client["drug_refined_data"] # 新的精选数据库

    def _check_lipinski(self, mol):
        """类药五原则检查"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        return violations <= 1

    # ==========================================
    # 🧪 1. 处理 PubChem -> refined_pubchem
    # ==========================================
    def refine_pubchem(self):
        col_raw = self.db_raw["raw_pubchem"]
        col_ref = self.db_refined["refined_pubchem"]
        col_ref.drop() # 重置精选表
        col_ref.create_index("inchikey", unique=True)

        raw_data = list(col_raw.find())
        logging.info(f"🧪 正在精炼 PubChem (结构与类药性), 共 {len(raw_data)} 条...")

        for doc in tqdm(raw_data, desc="Refining PubChem"):
            data = doc.get("data", {})
            smiles = data.get("smiles") or data.get("isomeric_smiles")
            if not smiles: continue
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol: continue

            # 计算核心评价指标
            qed_val = QED.qed(mol)
            sa_score = sascorer.calculateScore(mol) if sascorer else None
            
            refined_doc = {
                "name": doc["query_name"],
                "inchikey": data.get("inchikey"),
                "smiles_canonical": Chem.MolToSmiles(mol, canonical=True),
                "molecular_props": {
                    "mw": Descriptors.MolWt(mol),
                    "logp": Descriptors.MolLogP(mol),
                    "tpsa": Descriptors.TPSA(mol),
                    "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol)
                },
                "quality_scores": {
                    "qed": qed_val,
                    "sa_score": sa_score,
                    "is_drug_like": (qed_val > 0.5 and self._check_lipinski(mol))
                },
                "meta": {"source": "PubChem", "timestamp": time.time()}
            }
            try:
                col_ref.insert_one(refined_doc)
            except: pass # 忽略重复项

    # ==========================================
    # 🧬 2. 处理 ChEMBL -> refined_chembl
    # ==========================================
    def refine_chembl(self):
        col_raw = self.db_raw["raw_chembl"]
        col_ref = self.db_refined["refined_chembl"]
        col_ref.drop()

        raw_data = list(col_raw.find())
        logging.info(f"🧬 正在精炼 ChEMBL (实验活性与临床), 共 {len(raw_data)} 条...")

        for doc in tqdm(raw_data, desc="Refining ChEMBL"):
            data = doc.get("data", {})
            mol_rec = data.get("molecule_full_record", {})
            activities = data.get("all_activities", [])
            
            # 只提取高质量活性数据 (IC50/Ki 且有 pChEMBL 值)
            valid_acts = []
            for act in activities:
                p_val = act.get("pchembl_value")
                if p_val and act.get("standard_type") in ["IC50", "Ki"]:
                    valid_acts.append({
                        "target_id": act.get("target_chembl_id"),
                        "pchembl_value": float(p_val),
                        "type": act.get("standard_type"),
                        "relation": act.get("standard_relation")
                    })
            
            if not valid_acts: continue

            refined_doc = {
                "name": doc["query_name"],
                "chembl_id": mol_rec.get("molecule_chembl_id"),
                "max_phase": mol_rec.get("max_phase"),
                "is_approved": mol_rec.get("max_phase") == 4,
                "activities": valid_acts,
                "activity_count": len(valid_acts),
                "avg_pchembl": sum([a['pchembl_value'] for a in valid_acts]) / len(valid_acts)
            }
            col_ref.insert_one(refined_doc)

    # ==========================================
    # 💎 3. 处理 PDB -> refined_pdb
    # ==========================================
    def refine_pdb(self):
        col_raw = self.db_raw["raw_pdb"]
        col_ref = self.db_refined["refined_pdb"]
        col_ref.drop()

        raw_data = list(col_raw.find())
        logging.info(f"💎 正在精炼 PDB (蛋白质质量校验), 共 {len(raw_data)} 条...")

        for doc in tqdm(raw_data, desc="Refining PDB"):
            data = doc.get("data", {})
            # 路径: rcsb_entry_info.resolution_combined
            res = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
            # 检查配体 (是否有钥匙孔)
            ligands = data.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])
            
            refined_doc = {
                "pdb_id": doc["query_id"],
                "title": data.get("struct", {}).get("title"),
                "resolution": res,
                "method": data.get("exptl", [{}])[0].get("method"),
                "has_ligands": len(ligands) > 0,
                "ligand_count": len(ligands),
                "is_high_quality": (res is not None and res <= 2.5)
            }
            col_ref.insert_one(refined_doc)

    def run_all(self):
        print("\n" + "="*40)
        print("🚀 MolNexus 数据精炼流水线启动")
        print("="*40)
        self.refine_pubchem()
        self.refine_chembl()
        self.refine_pdb()
        print("\n🎉 精炼任务完成！数据已存入数据库: drug_refined_data")

if __name__ == "__main__":
    refiner = MolNexusRefiner()
    refiner.run_all()