# FASTQ 混样拆分服务（Demultiplexer）

基于 FastAPI 的网页服务，把混样测序的 R1 / R2 / I1 / I2 四行制
`FASTQ.gz` 按样本索引拆回各样本，产出可下载的 ZIP 交给下游比对流水线。

## 拆分规则

- 上传文件：`R1`、`R2`、`I1`、`I2` 四个 gzip 压缩的四行制 FASTQ，
  以及三列样本表 CSV（表头必须为 `sample_id,i7,i5`）。
- 样本索引等长、`i7+i5` 组合唯一；`i5` 按**匹配方向**填写（服务端不再做反向互补）。
- 四路读段名取**首个空白符之前**的部分，并去掉末尾 `/1`、`/2` 后对齐；
  任一条目四路名字不一致即整批失败。
- `I1` 比对样本 `i7`、`I2` 比对样本 `i5`，每端汉明距离 **≤ 1** 才进入候选；
  在候选中取**总距离最小且唯一**的样本。无候选或最优并列时，该读段进入
  `undetermined`（未分配）组。
- 以下情况**整批失败**（状态 `failed`，不提供下载）：
  - 四路记录条数不同；
  - 读段名错序 / 对不齐；
  - 序列与质量字符串不等长、FASTQ 结构损坏；
  - gzip 流损坏（魔数错误、截断、CRC/ISIZE 校验失败）；
  - 样本表不合法。
- 逐条流式处理：读一条、判一条、写一条，不在内存攒整批读段。
- 只有**全部读段校验并保存成功**后才生成并开放 ZIP。ZIP 内保留各组
  R1、R2 的**原始字节与组内顺序**，并附 `assignment_counts.csv` 分配数。

## 作业与并发

- 最多 **2 批并行**；已满时新批次直接返回 `409`（不排队）。
- 网页可启动、查看（每 2 秒自动刷新）、取消作业。
- 取消接口会阻塞到本批的实际读写线程真正停止后才返回，并立即释放并行名额；
  另一批照常运行，新批可立刻接替。被取消的任务不会再写文件、不发布结果，
  其工作目录与上传件全部删除。
- 失败或取消的作业一律不提供下载（下载接口返回 `404`）。
- 任务状态**只存内存**，进程重启即清空，启动时同时清空磁盘暂存目录。

## 目录结构

```
app/
  main.py        FastAPI 路由（上传 / 列表 / 详情 / 取消 / 下载 / 首页）
  jobs.py        内存态作业注册表、双批信号量、取消与清理
  pipeline.py    流式拆分主流程、分组 gzip 输出、ZIP 打包
  fastq_io.py    gzip FASTQ 流式读取器与损坏检测
  demux.py       样本表解析、汉明距离、候选/并列判定
  errors.py      领域异常
  templates/     单页前端
scripts/
  generate_samples.py   生成三套演示数据
samples/
  normal/      正常样例（全部唯一命中，含 1 错配）
  ambiguous/   歧义样例（含无候选与并列最优，各 1 条进 undetermined）
  corrupt/     末尾截断样例（R2 gzip 尾部损坏，整批失败）
tests/         pytest 测试（25 项）
```

## 本地运行（无 Docker）

```bash
python3 -m pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
# 可选：export DEMUX_JOBS_ROOT=/data/jobs   # 暂存目录，默认 /data/jobs
```

浏览器打开 <http://localhost:8000>，从 `samples/<场景>/` 选择五个文件上传。

重新生成演示数据：

```bash
python3 -I scripts/generate_samples.py samples
```

## Docker Compose 运行

```bash
docker compose up --build
# 浏览器打开 http://localhost:8000
```

容器以非 root 用户运行，暂存文件写入命名卷 `demux-jobs`（挂载于
`/data/jobs`），每次启动自动清空；容器内置健康检查。

## HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `GET`  | `/` | 上传与作业监控网页 |
| `GET`  | `/api/health` | 健康检查 |
| `POST` | `/api/jobs` | multipart 上传 5 个文件并启动作业，`201`；并发已满返回 `409` |
| `GET`  | `/api/jobs` | 全部作业快照 |
| `GET`  | `/api/jobs/{id}` | 单个作业快照（状态、进度、分配数、错误信息） |
| `POST` | `/api/jobs/{id}/cancel` | 取消并等待实际停止后返回 |
| `GET`  | `/api/jobs/{id}/download` | 仅完成态可下载 ZIP，否则 `404` |

## 测试

```bash
python3 -m pytest tests/ -q
```

覆盖：索引匹配（完美/单端错配/双端各一错配/并列/无候选/长度不符）、
样本表校验、gzip 与 FASTQ 损坏、条数不一致、错序、序质不等长、
ZIP 内容与原始字节顺序、双批上限 `409`、取消后名额释放与新批接替、
失败/取消禁止下载、内存态 404、首页渲染。
