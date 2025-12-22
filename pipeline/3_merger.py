import pymongo
import time
import sys
import os
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from tqdm import tqdm

# --- 1. 环境配置 ---
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 2. 数据库连接 ---
try:
    client = pymongo.MongoClient(config.MONGO_URI)
    
    db_raw = client[config.DB_RAW]
    col_pubchem = db_raw["raw_pubchem"]
    col_chembl = db_raw["raw_chembl"]
    col_pdb = db_raw["raw_pdb"]
    
    db_core = client["drug_core"]
    col_merged = db_core["merged_rough_data"]
    
    # 【修复步骤 1】先删除旧索引，防止之前的错误残留
    try:
        col_merged.drop_index("identity.std_inchi_key_1")
        logging.info("🧹 已清理旧索引")
    except:
        pass # 如果索引本来就不存在，忽略报错
    
    # 【修复步骤 2】重新创建稀疏唯一索引
    # sparse=True 的意思是：只有当文档里包含这个字段时，才检查唯一性。不包含就不检查。
    col_merged.create_index("identity.std_inchi_key", unique=True, sparse=True)
    col_merged.create_index("identity.primary_name")
    
    logging.info(f"🔌 数据库连接成功: {config.MONGO_URI}")

except Exception as e:
    logging.error(f"❌ 数据库连接失败: {e}")
    sys.exit(1)


def generate_3d_conformer(mol):
    """生成 3D 构象"""
    try:
        mol_h = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG(randomSeed=42))
        if res == 0:
            AllChem.MMFFOptimizeMolecule(mol_h)
            return Chem.MolToMolBlock(mol_h), "Computed_RDKit"
    except:
        pass
    return None, "Failed"

def process_small_molecules():
    """处理小分子"""
    raw_mols = list(col_pubchem.find({}))
    # 过滤掉没有 SMILES 的数据再计数，让进度条更准
    valid_mols = [m for m in raw_mols if m.get("data", {}).get("isomeric_smiles")]
    
    logging.info(f"🧪 开始处理小分子，有效记录 {len(valid_mols)} 条...")

    success_count = 0
    for p_doc in tqdm(valid_mols, desc="Processing Molecules", unit="mol"):
        name = p_doc.get("query_name")
        p_data = p_doc.get("data", {})
        input_smiles = p_data.get("isomeric_smiles")

        mol = Chem.MolFromSmiles(input_smiles)
        if not mol: continue
        
        std_inchi_key = Chem.MolToInchiKey(mol)
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        mol_block, str_source = generate_3d_conformer(mol)

        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "qed": QED.qed(mol)
        }

        # 融合 ChEMBL
        c_doc = col_chembl.find_one({"query_name": name})
        bio_targets = []
        chembl_id = None
        
        if c_doc and c_doc.get("data"):
            c_data = c_doc.get("data", {})
            if "molecule_info" in c_data:
                chembl_id = c_data["molecule_info"].get("molecule_chembl_id")
            if "bio_activities" in c_data:
                for act in c_data["bio_activities"]:
                    if act.get("standard_value") and act.get("standard_type") in ["IC50", "Ki", "Kd", "EC50"]:
                        bio_targets.append({
                            "target_id": act.get("target_chembl_id"),
                            "type": act.get("standard_type"),
                            "value": float(act.get("standard_value")),
                            "unit": act.get("standard_units")
                        })

        merged_doc = {
            "type": "Small_Molecule",
            "identity": {
                "primary_name": name,
                "std_inchi_key": std_inchi_key, # 小分子有这个字段
                "ids": {
                    "pubchem_cid": p_data.get("cid"),
                    "chembl_id": chembl_id
                }
            },
            "structure": {
                "smiles_clean": canonical_smiles,
                "mol_block_3d": mol_block,
                "3d_source": str_source
            },
            "properties": props,
            "bio_activity": {
                "targets": bio_targets,
                "count": len(bio_targets)
            },
            "meta": {
                "source_pipeline": "MolNexus_V1",
                "last_updated": time.time()
            }
        }

        try:
            col_merged.update_one(
                {"identity.std_inchi_key": std_inchi_key}, 
                {"$set": merged_doc}, 
                upsert=True
            )
            success_count += 1
        except Exception as e:
            pass

    print(f"\n✅ 小分子处理完成，入库 {success_count} 条。")


def process_proteins():
    """处理蛋白质"""
    raw_pdbs = list(col_pdb.find({}))
    logging.info(f"\n🧬 开始处理蛋白质，共 {len(raw_pdbs)} 个记录...")
    
    success_count = 0
    for doc in tqdm(raw_pdbs, desc="Processing Proteins ", unit="pdb"):
        pdb_id = doc.get("query_id")
        data = doc.get("data", {})
        
        info = data.get("rcsb_entry_info", {})
        struct = data.get("struct", {})
        
        merged_doc = {
            "type": "Protein_Target",
            "identity": {
                "primary_name": pdb_id,
                # 【修复步骤 3】这里彻底删掉了 "std_inchi_key": None
                # 因为只有字段不存在，sparse 索引才会忽略它
                "ids": {
                    "pdb_id": pdb_id,
                }
            },
            "structure": {
                "title": struct.get("title"),
                "method": data.get("exptl", [{}])[0].get("method"),
                "resolution": info.get("resolution_combined", [None])[0]
            },
            "properties": {
                "polymer_count": info.get("polymer_entity_count_protein"),
                "atom_count": info.get("atom_count_lebedev")
            },
            "meta": {
                "source_pipeline": "MolNexus_V1",
                "last_updated": time.time()
            }
        }
        
        try:
            col_merged.update_one(
                {"identity.primary_name": pdb_id},
                {"$set": merged_doc},
                upsert=True
            )
            success_count += 1
        except Exception as e:
            logging.error(f"Error saving {pdb_id}: {e}")
        
    print(f"\n✅ 蛋白质处理完成，入库 {success_count} 条。")

if __name__ == "__main__":
    print(f"{'='*30} MolNexus Merger 启动 {'='*30}")
    process_small_molecules()
    process_proteins()
    print(f"\n🎉 全部清洗与融合完毕！")

