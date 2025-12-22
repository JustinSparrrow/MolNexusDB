import pymongo
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import requests
import time
import sys
import os

# 引入配置
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# --- 连接数据库 ---
client = pymongo.MongoClient(config.MONGO_URI)
db_staging = client[config.DB_STAGING]
db_raw = client[config.DB_RAW]

col_queue = db_staging["task_queue"]
col_pubchem = db_raw["raw_pubchem"]
col_chembl = db_raw["raw_chembl"]  
col_pdb = db_raw["raw_pdb"]

# 初始化 ChEMBL 客户端
chembl_mol = new_client.molecule
chembl_act = new_client.activity

def pubchem_to_dict(compound):
    """
    专门解决序列化问题：将 PubChem 对象手动转为 MongoDB 可存的纯字典
    """
    c_dict = compound.to_dict()
    # 处理 atoms: 将对象转为普通字典列表
    if compound.atoms:
        c_dict['atoms'] = [
            {'aid': a.aid, 'element': a.element, 'x': a.x, 'y': a.y, 'z': a.z} 
            for a in compound.atoms
        ]
    # 处理 bonds: 同理
    if compound.bonds:
        c_dict['bonds'] = [
            {'aid1': b.aid1, 'aid2': b.aid2, 'order': b.order} 
            for b in compound.bonds
        ]
    return c_dict

def fetch_pubchem_data(name):
    """全量抓取 PubChem 数据"""
    try:
        compounds = pcp.get_compounds(name, namespace='name')
        if compounds:
            # 使用我们写的转换函数，确保数据是“纯净”的 JSON 格式
            return pubchem_to_dict(compounds[0])
    except Exception as e:
        print(f"      ⚠️ PubChem 抓取异常 ({name}): {e}")
    return None

def fetch_chembl_data(name, smiles=None):
    """全量抓取 ChEMBL 数据"""
    try:
        res = None
        if smiles:
            # 使用 flexmatch 寻找最匹配的结构
            res = chembl_mol.filter(molecule_structures__canonical_smiles__flexmatch=smiles)

        if not res or len(res) == 0:
            res = chembl_mol.search(name)
            
        if res and len(res) > 0:
            best_match = res[0] 
            chembl_id = best_match['molecule_chembl_id']
            
            # 抓取活性数据并强制转化为 list，确保 MongoDB 能存
            activities = list(chembl_act.filter(molecule_chembl_id=chembl_id)[0:50])
            
            return {
                "molecule_full_record": best_match,
                "all_activities": activities 
            }
    except Exception as e:
        print(f"      ⚠️ ChEMBL 抓取异常 ({name}): {e}")
    return None

def process_tasks():
    # 只拿 pending 状态的任务
    tasks = list(col_queue.find({"status": "pending"}))
    total = len(tasks)
    
    if total == 0:
        print(" Ø 没有发现等待处理的任务。请运行 loader 或重置任务状态。")
        return

    print(f"📦 发现 {total} 个待处理任务，开始全量抓取...\n")

    for i, task in enumerate(tasks):
        name = task["search_name"]
        category = task["category"]
        task_id = task["_id"]
        
        print(f"[{i+1}/{total}] 正在抓取: {name} ({category})...")

        try:
            success_any = False # 标记是否至少从一个源拿到了数据

            if category == "MOL":
                # 1. 尝试 PubChem
                pc_data = fetch_pubchem_data(name)
                current_smiles = None
                
                if pc_data:
                    col_pubchem.update_one(
                        {"query_name": name}, 
                        {"$set": {"data": pc_data, "updated_at": time.time()}}, 
                        upsert=True
                    )
                    print("      ✅ PubChem: 成功")
                    current_smiles = pc_data.get('isomeric_smiles')
                    success_any = True
                
                # 2. 尝试 ChEMBL
                cb_data = fetch_chembl_data(name, smiles=current_smiles)
                if cb_data:
                    col_chembl.update_one(
                        {"query_name": name},
                        {"$set": {"data": cb_data, "updated_at": time.time()}},
                        upsert=True
                    )
                    print(f"      ✅ ChEMBL: 成功 (ID: {cb_data['molecule_full_record']['molecule_chembl_id']})")
                    success_any = True

            elif category == "PDB":
                url = f"https://data.rcsb.org/rest/v1/core/entry/{name.lower()}"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    col_pdb.update_one(
                        {"query_id": name},
                        {"$set": {"data": response.json(), "updated_at": time.time()}},
                        upsert=True
                    )
                    print("      ✅ RCSB PDB: 成功")
                    success_any = True

            # 更新任务状态
            if success_any:
                col_queue.update_one({"_id": task_id}, {"$set": {"status": "done"}})
            else:
                col_queue.update_one({"_id": task_id}, {"$set": {"status": "failed", "reason": "No data found in any source"}})

            # 稍微停顿，防止被封 IP
            time.sleep(0.8)

        except Exception as e:
            print(f"   ❌ 运行时严重错误 ({name}): {e}")
            col_queue.update_one({"_id": task_id}, {"$set": {"status": "error", "error_msg": str(e)}})

    print("\n🎉 处理循环结束。请在 MongoDB Compass 中检查数据量。")

if __name__ == "__main__":
    process_tasks()

