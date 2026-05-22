**设计概要：目标 · 技术方案 · 数据准备 · 代码架构 · 验证方案**

**1 目标（Objectives）**
- 明确目标：构建一套可复用的 Xenium → scSpace 数据处理与分析流水线，输出伪空间（pseudo-space）表征并提供可评估的聚类与空间一致性指标。
- 可重用性：模块化设计，便于替换输入来源（细胞级 Xenium 原始表或已有 sc-like / st-like 表视图）。
- 可配置性：暴露关键超参数（如 bin-size、n_features、训练轮数、学习率、聚类参数），支持重现实验与敏感性分析。

**2 技术方案（Technical Approach）**
- 输入转换：实现从 Xenium 细胞级计数与坐标到 sc-like 与 st-like bundle 的转换，支持按像素/spot 聚合（可选的 binning/aggregation 算法）。
- 预处理：统一计数矩阵方向（基因×细胞/spot）、标准化、选择高变基因（HVG）或交集策略（`select_hvg`）。
- Pseudo-space 构建：基于 `scSpace.construct_pseudo_space`，使用小批量训练（`batch_size`）、激活函数（如 sigmoid）、优化器（指定学习率）训练网络以映射细胞到伪空间。
- 空间聚类：在生成的伪空间上运行空间聚类（`spatial_cluster`），使用 `Ks`, `Kg`, `alpha`, `res` 等参数调整宏/微观簇结构。
- 评估与可视化：计算距离分布、簇内一致性、spot-level assignment 准确度等指标，并生成可比图表（例如与 `half_split_benchmark` 的可视化结果做对比）。

**3 数据准备（Data Preparation）**
- 输入格式：需要四个文件：`<prefix>_sc_data.csv`, `<prefix>_sc_meta.csv`, `<prefix>_st_data.csv`, `<prefix>_st_meta.csv`（基因×样本矩阵与以样本 id 为索引的元数据）。
- 从 Xenium 原始表准备：解析细胞级基因计数表与坐标表，进行质量控制（细胞/spot 最小计数阈值、去除低表达基因），可选聚合到 Visium-like spots：指定 `bin-size` 与 `min-cells-per-bin`。
- 元数据对齐：确保元数据索引与计数矩阵列/行完全一致；空间坐标键（如 `xcoord`, `ycoord`）应在 meta 中明确。
- 小样本/稀疏处理：记录行/列稀疏度，若存在大量零样本则在预处理阶段过滤或补偿（normalize/scale 策略）。

**4 代码架构（Code Architecture）**
- 目录划分（高层）：
	- `xenium_scspace/cli.py`：对外命令行接口，提供 `prepare` 与 `run` 两个子命令。
	- `xenium_scspace/bundle.py`：bundle 路径管理与验证（`ScSpaceBundlePaths`, `bundle_files`, `validate_bundle`）。
	- `xenium_scspace/pipeline.py`：核心流水线（`run_scspace_pipeline`），封装加载、预处理、pseudo-space 构建、聚类与结果输出。
	- `xenium_scspace/advanced_analysis.py`：评估脚本与可视化函数，生成 benchmark 图表。
	- `tests/`：单元与集成测试（数据格式校验、端到端小数据运行）。
- 模块接口（关键函数与约定）：
	- `prepare(sc_counts, sc_meta, st_counts?, st_meta?, output_dir, bin_size, min_cells_per_bin)`：生成 bundle。
	- `run(bundle_dir, output_dir, st_type, **pipeline_args)`：触发 `run_scspace_pipeline` 并保存结果。
	- `run_scspace_pipeline(...)`（`pipeline.py`）返回字典包含 `sc_adata`, `st_adata`, `output_dir`（如写出）。
- 配置与可重复性：采用 YAML/CLI 参数传递超参数，记录随机种子与运行环境依赖快照（`pip freeze` 或 `pyproject.toml`）。

**5 验证方案（Validation Plan）**
- 输出核对：验证 `pseudo_space` CSV 存在且索引/列符合预期（IDs 对齐、坐标维度正确）。
- 指标计算：实现并记录以下指标
	- spot/cell-level距离分布（query vs. truth）
	- 聚类一致性（例如 ARI/NMI）与空间连贯性指标
	- spot-assignment 精度（若有 ground-truth）
- 可视化对比：重现 `half_split_benchmark/evaluation` 中的关键图表（距离分布图、assignment scatter、pseudo-space 可视化），并将数值与基线 JSON 比较（例如 `evaluation_metrics.json`）。
- 回归测试：在小规模合成数据上建立端到端测试，确保修改不会破坏主流程。

**7 交付物与文档**
- 更新 `design.md`（本文件）说明核心设计。
- 在 `_l_result/` 下保存示例运行结果与评估报告。
- 在 `README.md` 中增加快速入门与复现实验的最小命令集合。
