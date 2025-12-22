# 这一步是为了让工人重新干活
import pymongo
import sys, os
sys.path.append(os.path.abspath("pipeline")) # 确保能找到 config
import config

client = pymongo.MongoClient(config.MONGO_URI)
col = client[config.DB_STAGING]["task_queue"]

# 把所有任务重置为 pending
col.update_many({}, {"$set": {"status": "pending"}})
print("任务已重置，请重新运行 worker")