# 🚀 MolNexus: Multi-modal Molecular Gold Standard Dataset
## 🛠️ Environment Setup & Database Initialization | 环境配置与数据库启动
本项目基于 Docker 容器化部署，集成了 MongoDB (数据存储) 和 Mongo Express (可视化管理后台)。在运行 Python 脚本前，请务必确保数据库服务已正常启动。

### 1. Prerequisites (前置准备)
确保你的电脑已安装 Docker Desktop。

### 2. Start Services (启动服务)
Step 1: 启动 Docker 引擎
- Mac / Windows: 双击打开 Docker Desktop 应用程序，等待左下角状态变为绿色的 Engine running。
- Linux: 终端运行 sudo systemctl start docker。

Step 2: 拉起数据库容器  
在项目根目录下（即 docker-compose.yml 所在目录），打开终端运行：
```
docker-compose up -d
```
参数说明：-d (Detached) 表示在后台静默运行，不会占用当前终端窗口。

### 3. Health Check (状态检测)
启动命令执行后，请通过以下任一方式检查服务是否正常：

方式 A: 命令行检测 (推荐).  
运行以下命令查看容器状态：
```
docker ps
```

预期输出： 

你应该看到两个容器 my_drug_db 和 my_db_admin，且 STATUS 栏显示为 Up (例如 Up 2 minutes)。
- ❌ 如果显示 Restarting，说明配置有误，请检查 docker-compose.yml 或端口占用情况。
- ❌ 如果列表为空，说明容器未启动，请运行 docker-compose logs 查看报错。

方式 B: 可视化后台检测

打开浏览器访问：
👉 http://localhost:8081

如果你能看到 Mongo Express 的登录或管理界面，说明数据库服务已完美运行。
- Web 登录账号: admin
- Web 登录密码: password

