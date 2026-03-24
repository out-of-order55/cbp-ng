# 项目导入与执行规范

## 规则优先级

- `Must`：必须遵守，冲突时优先级最高。
- `Should`：强烈建议，若偏离需在计划中说明原因。
- `May`：可选优化项。

## 基本信息（Must）

1. 请用中文回答，英文写方案到文件。
2. 初次进入项目需按顺序阅读：
   - `@README.md`（编译、运行、项目概述）
   - `@cbp.hpp`、`@cbp.cpp`（模拟器执行路径）
   - `@harcom.hpp`、`@docs/harcom.pdf`（`val/ram/rom/reg/fanout/energy/floorplan`）
   - `@predictors/my_bp_v1.hpp`、`@predictors/tage.hpp`（预测器实现）
   - `@docs/vfs.pdf`、`@vfs.py`（评分机制）
3. 完成上述阅读后，回复：`我现在可以开始工作`。

## 工作流程（Must）

1. 每次开始先进入 plan 流程；若运行环境不支持显式 plan mode，也必须先给出计划。
2. 计划必须写入工作目录 `doc/` 子目录，文件内容用英文。
3. 修改后必须先编译，再运行指定 trace；若失败必须继续修复直到无错误。

## 变更类型与验证矩阵（Must）

### A. 修改 `my_bp_v1` 且属于计数器类改动

```bash
./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE -DPERF_COUNTERS -DCHEATING_MODE -DFREE_FANOUT
```

### B. 修改 `my_bp_v1` 且属于功能类改动（非计数器）

```bash
./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE
```

### C. 修改其他预测器文件

- 将 `-DPREDICTOR` 改为对应预测器名（例如 `tage<>`）。
- 编译输出名与后续运行命令保持一致。

### D. 统一运行命令（A/B/C 后都要执行）

```bash
./my_bp /nfs/home/guobing1/cbp-ng/trace/cbp-ng_training_traces/502-gcc-all_16112_trace.gz test 1000000 2000000 --format human
```

## Feature 开发要求（Should）

1. 对比修改前后：性能、功耗、延时。
2. 必须增加与 `tage` 的对比结果。
3. 报告中明确给出关键指标变化。

## 调试策略（Must）

1. 先复现最小失败命令。
2. 使用 `-DVERBOSE` 或相关计数器选项定位。
3. 若仍无法定位，使用 `gdb` 并记录关键 backtrace。

## 编码规范（Must）

1. 按硬件思维编程，显式考虑信号扇出；同一信号多处消费时必须检查是否需要设置 `fanout`。

```cpp
// 例子1
val<1> a = b;
// fanout = 1
c = a;
// fanout = 10
for (int i = 0; i < 10; i++) {
  c = a;
}

// 例子2
arr<val<1>, 10> a;
// fanout = 1
for (int i = 0; i < 10; i++) {
  b = a[i];
}
```

2. SRAM 读写条件需要放在一起：

```cpp
execute_if(was_used | do_alloc, [&](){
    gctr[w][j].write(gindex[j], new_ctr, extra_cycle);
});
```

3. `val` 不可二次赋值。
4. 不允许 `arr<arr<...>>` 嵌套。
5. 默认不改已有函数的外部行为与签名；若必须修改，需在计划中说明影响与理由。
6. 生成文件不得写入 `/tmp`，只能放在工作目录及子目录。
7. 测试unit文件放入test文件夹，测试完成后删除可执行文件

## 交付清单（Must）

1. 本次改动文件列表。
2. 实际执行的编译命令与结果。
3. 实际执行的运行命令与结果。
4. 若为 feature：提供 before/after 与 `tage` 对比。

## WB 可见周期统计规范（Must）

1. 讨论“几周期后可见”时，必须明确并同时报告两种口径：
   - `probe_latency_cycles`：从读探测开始到首次命中的周期数。
   - `e2e_latency_cycles`：从写请求锚点周期到首次命中的端到端周期数。
2. 若只报告一种口径，结论无效。
3. 若使用随机延迟抖动，必须固定随机种子并记录在报告中。
4. 随机延迟实验至少覆盖：
   - 地址延迟 `0~4`
   - 数据延迟 `0~4`
   - 连续 fill (`noconflict=0`) + 连续 drain (`noconflict=1`)
5. 当回答“最大延迟”时，不允许引用单次运行；必须多种子统计并报告：
   - `max`
   - `p99`
   - `mean`
   - 样本数 `N`
6. 默认推荐最小统计规模：`N >= 256`（可由多 seed * 多请求组成）。
7. 报告字段至少包含：
   - RAM 类型（`wb_ram/wb_mask_ram/wb_rwram`）
   - 测试模式（direct/buffered/continuous/jitter）
   - 种子
   - `probe_latency_cycles` 与 `e2e_latency_cycles`
   - 是否出现 stall/drop/full
8. 若 `probe_latency` 与 `e2e_latency` 差异过大，优先排查：
   - 探测起点定义是否一致
   - drain 启动时机
   - 读探测是否违反 `single RAM read per cycle`

## WB 测试输出规范（Should）

1. 每条关键写请求建议打印：`idx / addr / data / write_anchor / first_visible / probe_latency / e2e_latency`。
2. 输出中建议附加 `queue_depth` 峰值与 drain 完成周期，便于快速定位瓶颈。
3. 回归对比时，优先比较同一 seed、同一请求序列，避免随机差异误导结论。
