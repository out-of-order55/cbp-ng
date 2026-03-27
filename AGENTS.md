# 项目导入与执行规范

## 规则优先级

- `Must`：必须遵守，冲突时优先级最高。
- `Should`：强烈建议，若偏离需在计划中说明原因。
- `May`：可选优化项。

## 基本信息（Must）

1. 请用中文回答，每次开发需要明确方案。
2. 初次进入项目需按顺序阅读：
   - `@README.md`（编译、运行、项目概述）
   - `@cbp.hpp`、`@cbp.cpp`（模拟器执行路径）
   - `@harcom.hpp`、`@docs/harcom.pdf`（`val/ram/rom/reg/fanout/energy/floorplan`）
   - `@predictors/my_bp_v1.hpp`、`@predictors/tage.hpp`（预测器实现）
   - `@docs/vfs.pdf`、`@vfs.py`（评分机制）
   - @CBP2016.h、`@docs/RR-9561.pdf`
3. 完成上述阅读后，回复：`我现在可以开始工作`。

## 工作流程（Must）

1. 每次开始先进入 plan 流程；若运行环境不支持显式 plan mode，也必须先给出计划。
2. 计划必须写入工作目录 `doc/` 子目录，文件内容用英文。
3. 修改后必须先编译，再运行指定 trace；若失败必须继续修复直到无错误。
4. 提交不需要提交doc，只提交更改的代码
5. 添加的性能计数器要使用函数封装

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
7. 测试unit文件放入test文件夹，测试完成后删除可执行文件，无论性能测试还是unit测试

## 交付清单（Must）

1. 本次改动文件列表。
2. 实际执行的编译命令与结果。
3. 实际执行的运行命令与结果。
4. 若为 feature：提供 before/after 与 `tage` 对比。

## BHYST 聚合 WB RAM 经验沉淀（Must）

1. 讨论“延迟是否变大”时，必须同时给出三种口径，缺一不可：
   - `macro_latency_cycles`：SRAM 宏模型 `static_ram::LATENCY`。
   - `intrinsic_latency_cycles`：地址预连接到目标 RAM 后的读延迟（尽量剥离线网影响）。
   - `e2e_floorplan_latency_cycles`：含 floorplan/wire 的端到端读延迟。
2. 将分散 SRAM 改成 `wb_ram` 后，必须避免“同周期后端 RAM 读+写”，否则会触发 `single RAM access per cycle`。
3. `wb_ram` 使用规则：有读且要写时优先入队（`noconflict=0`），无读周期再 drain（`noconflict=1`），不要在同周期强制直写后端。
4. `fo1()` 只能消费一次；同一信号被多个分支/统计复用时，先 `fanout` 再分别读取，避免 `misuse of fo1()`。
5. 结构改造必须读写双侧一起改：改了读取组织（如按行聚合）后，更新路径也必须同步改为按行读改写，禁止只改读不改写。
6. 提交前必须执行：
   - `git add <目标文件>` 精确暂存；
   - `git show --name-only -n 1` 复核提交内容；
   - 若混入非目标文件（如 `AGENTS.md`/`doc`），必须在提交前清理。

## table1_hyst SRAM 选型经验（Must）

1. 对 `table1_hyst`（1-bit hyst，按行更新）优先选 `agg + wb_ram`，不优先 `agg_mask + wb_mask_ram`。
2. 本项目实测（`LOGLB=6, LOGP1=14`）：
   - `split`：`READ=217ps`，`WRITE_READY=33ps`，`WRITE_PURE=41ps`，`PWR=5.672665mW`
   - `agg(wb_ram)`：`READ=166ps`，`WRITE_READY=177ps`，`WRITE_PURE=39ps`，`PWR=1.492450mW`
   - `agg_mask(wb_mask_ram)`：`READ=170ps`，`WRITE_READY=226ps`，`WRITE_PURE=60ps`，`PWR=2.494258mW`
3. 结论：
   - `agg(wb_ram)` 相对 `split`：读更快（`-23.50%`），功耗显著下降（`-73.69%`）。
   - `agg_mask(wb_mask_ram)` 相对 `agg(wb_ram)`：读略慢（`+2.41%`），写更慢（`WRITE_READY +27.68%`），功耗更高（`+67.13%`）。
4. 写延迟解读必须拆分两项并同时报告：
   - `WRITE_READY`：含读依赖（RMW 口径）
   - `WRITE_PURE`：不含读依赖（纯写逻辑口径）
5. 何时考虑 `wb_mask_ram`：
   - 仅在“宽向量+稀疏 lane 更新+队列内高概率同址合并”场景下再评估；
   - 若是 1-bit/低位宽覆盖式更新，`wb_mask_ram` 常为过度设计。
