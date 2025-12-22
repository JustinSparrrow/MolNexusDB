import pymongo
import sys
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from pipeline import config

def clear_core_db():
    """清空核心数据库"""
    client = pymongo.MongoClient(config.MONGO_URI)
    db = client["drug_core"]
    db["merged_rough_data"].drop()
    print("🗑️ [Core] 核心金标库已清空。")

def show_core_stats():
    """核心库深度统计"""
    client = pymongo.MongoClient(config.MONGO_URI)
    col = client["myplan_core"]["merged_rough_data"]
    
    total = col.count_documents({})
    high_quality = col.count_documents({"analysis.is_high_quality": True})
    with_3d = col.count_documents({"structure.mol_block_3d": {"$ne": None}})
    
    print("\n🏆 核心库(MolNexus_Core)统计:")
    print(f"  - 总分子数: {total}")
    print(f"  - 高质量分子(HQ): {high_quality}")
    print(f"  - 包含3d构象: {with_3d}")

if __name__ == "__main__":
    show_core_stats()
    opt = input("\n是否清空核心金标库(Core)？(y/n): ")
    if opt.lower() == 'y':
        clear_core_db()