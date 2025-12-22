# 装载器
import pymongo
import os
import sys

# --- 让脚本能找到同级目录下的 config.py ---
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config  # 导入配置

client = pymongo.MongoClient(config.MONGO_URI)
db = client[config.DB_STAGING]
col_queue = db["task_queue"]

def load_seeds_to_mongo():
    # 使用配置里的路径，不管在哪里运行脚本都能找到文件
    filepath = config.SEEDS_FILE_PATH 
    
    print(f"📄 正在读取种子文件: {filepath}")
    
    if not os.path.exists(filepath):
        print("❌ 错误：找不到 seeds.txt 文件，请检查路径！")
        return

    col_queue.drop() 
    print("🧹 已清空旧的任务队列，准备装载新货物...")

    new_tasks = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split("|")
            name = parts[0]
            category = parts[1] if len(parts) > 1 else "MOL"

            task = {
                "search_name": name,
                "category": category,
                "status": "pending",
                "created_at": os.times().elapsed
            }
            new_tasks.append(task)

    if new_tasks:
        col_queue.insert_many(new_tasks)
        print(f"✅ 成功装载 {len(new_tasks)} 个种子任务！")
    else:
        print("⚠️ 文件是空的。")

if __name__ == "__main__":
    load_seeds_to_mongo()
