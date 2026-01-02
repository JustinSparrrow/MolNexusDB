import pymongo
import pubchempy as pcp
import requests
import time
import sys
import os
import logging

# --- 1. 环境与配置 ---
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 2. 数据库连接 ---
client = pymongo.MongoClient(config.MONGO_URI)
db_staging = client[config.DB_STAGING]
db_raw = client[config.DB_RAW]

col_queue = db_staging["task_queue"]
col_pubchem = db_raw["raw_pubchem"]
col_chembl = db_raw["raw_chembl"]
col_pdb = db_raw["raw_pdb"]

# --- 3. 核心抓取逻辑 ---

def pubchem_to_dict(compound):
    """PubChem 序列化处理"""
    c_dict = compound.to_dict()
    if compound.atoms:
        c_dict['atoms'] = [{'aid': a.aid, 'element': a.element, 'x': a.x, 'y': a.y, 'z': a.z} for a in compound.atoms]
    if compound.bonds:
        c_dict['bonds'] = [{'aid1': b.aid1, 'aid2': b.aid2, 'order': b.order} for b in compound.bonds]
    return c_dict

def fetch_chembl_via_rest(name):
    """
    通过底层 REST API 深度抓取 ChEMBL 数据。
    不再使用官方库，直接请求 URL 以保证数据最全。
    """
    try:
        # 第一步：根据名字搜索获取 ChEMBL ID
        search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
        search_res = requests.get(search_url, timeout=20).json()
        
        molecules = search_res.get('molecules', [])
        if not molecules:
            return None
            
        # 找到匹配度最高的一个 ID
        chembl_id = molecules[0]['molecule_chembl_id']
        logging.info(f"      🔗 识别到 ChEMBL ID: {chembl_id}")

        # 第二步：直接请求该 ID 的全量详情页面 (核心修复点)
        full_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
        full_record = requests.get(full_url, timeout=20).json()
        
        # 第三步：获取该分子的所有活性实验记录
        act_url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={chembl_id}&limit=100"
        activities_res = requests.get(act_url, timeout=20).json()
        activities = activities_res.get('activities', [])

        return {
            "molecule_full_record": full_record, # 这里包含 cross_references
            "all_activities": activities
        }
    except Exception as e:
        logging.warning(f"ChEMBL REST 抓取异常 ({name}): {e}")
    return None

def fetch_pubchem_deep(name):
    try:
        compounds = pcp.get_compounds(name, namespace='name')
        return pubchem_to_dict(compounds[0]) if compounds else None
    except: return None

def fetch_pdb_full(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        r = requests.get(url, timeout=15)
        return r.json() if r.status_code == 200 else None
    except: return None

# --- 4. 执行逻辑 ---

def start_worker():
    tasks = list(col_queue.find({"status": "pending"}))
    if not tasks:
        logging.info("🏁 没有待处理任务。")
        return

    for task in tasks:
        name = task["search_name"]
        category = task["category"]
        logging.info(f"🔎 挖掘中: {name}")
        
        try:
            has_data = False
            if category == "MOL":
                # 1. 抓 PubChem
                pc_data = fetch_pubchem_deep(name)
                if pc_data:
                    col_pubchem.update_one({"query_name": name}, {"$set": {"data": pc_data, "updated_at": time.time()}}, upsert=True)
                    has_data = True
                
                # 2. 直接用 REST 接口抓 ChEMBL (绕过官方库)
                cb_data = fetch_chembl_via_rest(name)
                if cb_data:
                    col_chembl.update_one({"query_name": name}, {"$set": {"data": cb_data, "updated_at": time.time()}}, upsert=True)
                    has_data = True

            elif category == "PDB":
                pdb_data = fetch_pdb_full(name)
                if pdb_data:
                    col_pdb.update_one({"query_id": name}, {"$set": {"data": pdb_data, "updated_at": time.time()}}, upsert=True)
                    has_data = True

            col_queue.update_one({"_id": task["_id"]}, {"$set": {"status": "done" if has_data else "failed"}})
            time.sleep(1.5) # 稍微加长间隔，防止被 REST API 屏蔽

        except Exception as e:
            logging.error(f"💥 崩溃 ({name}): {e}")
            col_queue.update_one({"_id": task["_id"]}, {"$set": {"status": "error"}})

if __name__ == "__main__":
    start_worker()