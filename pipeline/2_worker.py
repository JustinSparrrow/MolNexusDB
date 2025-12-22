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

def fetch_pubchem_data(name):
    """抓取 PubChem 基础数据"""
    try:
        # namespace='name' 允许通过俗名搜索 (如 Aspirin)
        compounds = pcp.get_compounds(name, namespace='name')
        if compounds:
            c = compounds[0]
            data = c.to_dict(properties=['atoms', 'bonds', 'isomeric_smiles', 'molecular_weight'])
            # 补充存一个 CID，方便以后索引
            data['cid'] = c.cid 
            return data
    except Exception as e:
        print(f"      ⚠️ PubChem 报错: {e}")
    return None

def fetch_chembl_data(name, smiles=None):
    """抓取 ChEMBL 活性数据 (比较慢，耐心等)"""
    try:
        # 策略：优先用 SMILES 搜 (更准)，如果没有 SMILES (比如 PubChem 没查到) 再用名字搜
        res = None
        if smiles:
            res = chembl_mol.filter(molecule_structures__canonical_smiles__flexmatch=smiles).only(['molecule_chembl_id', 'molecule_properties'])
        
        if not res or len(res) == 0:
            # 备选：用名字模糊搜索
            res = chembl_mol.search(name).only(['molecule_chembl_id', 'molecule_properties'])
            
        if res and len(res) > 0:
            best_match = res[0]
            chembl_id = best_match['molecule_chembl_id']
            
            # 进阶：顺便把这个分子的活性数据 (IC50 等) 也抓几条回来
            # 这里的 filter 意思是：查这个分子 ID，且类型是 IC50 的数据
            activities = chembl_act.filter(molecule_chembl_id=chembl_id, standard_type="IC50").only(['standard_value', 'standard_units', 'target_chembl_id'])[0:5]
            
            return {
                "molecule_info": best_match,
                "bio_activities": list(activities) # 转成列表存起来
            }
    except Exception as e:
        print(f"      ⚠️ ChEMBL 报错: {e}")
    return None

def process_tasks():
    # 查找待处理任务
    tasks = list(col_queue.find({"status": "pending"}))
    total = len(tasks)
    print(f"📦 发现 {total} 个待处理任务，工人开始干活...\n")

    for i, task in enumerate(tasks):
        name = task["search_name"]
        category = task["category"]
        task_id = task["_id"]
        
        print(f"[{i+1}/{total}] 正在处理: {name} ({category})...")

        try:
            # === 分支 A: 小分子 (MOL) -> 同时抓 PubChem 和 ChEMBL ===
            if category == "MOL":
                # 1. 先去 PubChem
                pc_data = fetch_pubchem_data(name)
                current_smiles = None
                
                if pc_data:
                    col_pubchem.update_one(
                        {"query_name": name}, 
                        {"$set": {"data": pc_data, "updated_at": time.time()}}, 
                        upsert=True
                    )
                    print("      ✅ PubChem: 获取成功")
                    current_smiles = pc_data.get('isomeric_smiles') # 拿到 SMILES 给 ChEMBL 用
                else:
                    print("      ❌ PubChem: 未找到")

                # 2. 拿着 SMILES 去 ChEMBL 查活性 (实现多源数据融合的第一步)
                # 即使 PubChem 没查到，也尝试用名字去 ChEMBL 碰碰运气
                cb_data = fetch_chembl_data(name, smiles=current_smiles)
                
                if cb_data:
                    col_chembl.update_one(
                        {"query_name": name},
                        {"$set": {"data": cb_data, "updated_at": time.time()}},
                        upsert=True
                    )
                    print(f"      ✅ ChEMBL: 获取成功 (ID: {cb_data['molecule_info']['molecule_chembl_id']})")
                else:
                    print("      ❌ ChEMBL: 未找到")

            # === 分支 B: 蛋白质 (PDB) -> 抓 RCSB ===
            elif category == "PDB":
                url = f"https://data.rcsb.org/rest/v1/core/entry/{name.lower()}"
                response = requests.get(url)
                if response.status_code == 200:
                    col_pdb.update_one(
                        {"query_id": name},
                        {"$set": {"data": response.json(), "updated_at": time.time()}},
                        upsert=True
                    )
                    print("      ✅ RCSB PDB: 获取成功")
                else:
                    print(f"      ❌ RCSB PDB: 下载失败 {response.status_code}")

            # 标记完成
            col_queue.update_one({"_id": task_id}, {"$set": {"status": "done"}})
            
            # 稍微睡 0.5 秒，防止请求太快被封 IP
            time.sleep(0.5)

        except Exception as e:
            print(f"   ❌ 严重错误: {e}")
            col_queue.update_one({"_id": task_id}, {"$set": {"status": "error", "error_msg": str(e)}})

    print("\n🎉 所有任务处理完毕！请去 MongoDB 查看 raw_pubchem, raw_chembl, raw_pdb 集合。")

if __name__ == "__main__":
    process_tasks()