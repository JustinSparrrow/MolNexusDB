import pymongo
import sys
import os

# 确保能找到配置
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from pipeline import config

def get_client():
    return pymongo.MongoClient(config.MONGO_URI)

def clear_refined_db():
    """清空精炼层数据库 (Refined Layer)"""
    client = get_client()
    db = client["drug_refined_data"]
    
    # 获取该库下所有的精炼表名
    collections = ["refined_pubchem", "refined_chembl", "refined_pdb"]
    
    for col_name in collections:
        db[col_name].drop()
        print(f"🗑️ [Refined] 集合已清空: {col_name}")
    print("✅ 精炼层数据库已完全重置。")

def show_refined_stats():
    """精炼层数据质量统计"""
    client = get_client()
    db = client["drug_refined_data"]
    
    print("\n🥈 --- MolNexus 精炼层 (Refined) 状态统计 ---")
    
    # 1. PubChem 统计
    col_pc = db["refined_pubchem"]
    pc_total = col_pc.count_documents({})
    pc_drug_like = col_pc.count_documents({"quality_scores.is_drug_like": True})
    print(f"🧪 PubChem (结构): {pc_total} 条 (类药性合格: {pc_drug_like})")
    
    # 2. ChEMBL 统计
    col_cb = db["refined_chembl"]
    cb_total = col_cb.count_documents({})
    cb_approved = col_cb.count_documents({"max_phase": 4})
    print(f"🧬 ChEMBL  (活性): {cb_total} 条 (已上市药物: {cb_approved})")
    
    # 3. PDB 统计
    col_pdb = db["refined_pdb"]
    pdb_total = col_pdb.count_documents({})
    pdb_hq = col_pdb.count_documents({"is_high_quality": True})
    print(f"💎 PDB     (靶点): {pdb_total} 条 (高质量结构: {pdb_hq})")
    
    print("-------------------------------------------\n")

if __name__ == "__main__":
    # 1. 先展示当前状态
    show_refined_stats()
    
    # 2. 交互式确认
    opt = input("⚠️ 是否要【清空】以上所有精炼后的数据？(y/n): ")
    if opt.lower() == 'y':
        confirm = input("❗ 请再次输入 'REFINED' 以确认彻底删除: ")
        if confirm == "REFINED":
            clear_refined_db()
        else:
            print("❌ 验证失败，操作取消。")
    else:
        print("操作已取消。")