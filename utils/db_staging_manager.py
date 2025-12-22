import pymongo
import sys
import os

# 确保能找到配置
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from pipeline import config

def reset_tasks():
    """将所有已完成或出错的任务重置为 pending"""
    client = pymongo.MongoClient(config.MONGO_URI)
    col = client[config.DB_STAGING]["task_queue"]
    res = col.update_many(
        {"status": {"$in": ["done", "error"]}}, 
        {"$set": {"status": "pending"}}
    )
    print(f"🔄 [Staging] 已重置 {res.modified_count} 条任务为等待处理。")

def show_queue_status():
    """查看任务队列详情"""
    client = pymongo.MongoClient(config.MONGO_URI)
    col = client[config.DB_STAGING]["task_queue"]
    counts = {
        "Total": col.count_documents({}),
        "Pending": col.count_documents({"status": "pending"}),
        "Done": col.count_documents({"status": "done"}),
        "Error": col.count_documents({"status": "error"})
    }
    print("\n📋 任务队列状态:")
    for k, v in counts.items():
        print(f"  - {k}: {v}")

if __name__ == "__main__":
    show_queue_status()
    opt = input("\n是否重置所有任务？(y/n): ")
    if opt.lower() == 'y':
        reset_tasks()