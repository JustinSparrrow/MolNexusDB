import pymongo
import sys
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from pipeline import config

def clear_raw_db():
    """清空原始数据库"""
    client = pymongo.MongoClient(config.MONGO_URI)
    db = client[config.DB_RAW]
    cols = db.list_collection_names()
    for c in cols:
        db[c].drop()
        print(f"🗑️ [Raw] 已删除集合: {c}")

def show_raw_stats():
    """统计下载了多少原始数据"""
    client = pymongo.MongoClient(config.MONGO_URI)
    db = client[config.DB_RAW]
    print("\n📦 原始库存储统计:")
    for col_name in ["raw_pubchem", "raw_chembl", "raw_pdb"]:
        count = db[col_name].count_documents({})
        print(f"  - {col_name}: {count} 条记录")

if __name__ == "__main__":
    show_raw_stats()
    opt = input("\n是否清空所有原始数据(Raw)？(y/n): ")
    if opt.lower() == 'y':
        clear_raw_db()