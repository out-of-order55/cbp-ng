# 项目导入

## 基本信息

1. please use english to answer promt
2. 阅读@README.md，了解项目如何编译和运行，以及项目的概述
3. 阅读@cbp.hpp,@cbp.cpp,了解模拟器如何工作
4. 阅读@harcom.hpp和@docs/harcom.pdf，了解harcom机制，包括val，ram，rom，reg，fanout，energy，floorplan
5. 阅读@predictors/my_bp.hpp和@predictors/tage.hpp和@predictors/my_bp_designtage.hpp.md,了解预测器是如何工作的
6. 阅读@docs/vfs.pdf和@vfs.py，了解评分机制
7. 如果上述步骤完成，请回复:我现在可以开始工作

## 工作流程

每次开始请以plan模式开始，先给出具体计划，写入文件，并且模型为Haiku,做完计划后代码请用sonnet4.6编写

每次修改文件后，需要运行./compile cbp -DPREDICTOR="my_bp<>" -DVERBOSE -DPERF_COUNTERS -DCHEATING_MODE -DUPDATEALTONWEAKMISP，保证没有错误，否则继续修改

之后运行./cbp ./gcc_test_trace.gz test 100 400 --format human，保证没有错误，否则继续修改

如果遇到错误，请使用gdb调试

## 编码规范

1. 每次编写时需要按照硬件思维编程，如一个val<1> a; b=a ,那么需要使得a的扇出为1：a.fanout{hard<1>{}},需要注意信号是否要设置fanout。
2. SRAM读写需要将读写条件放在一起，比如

```
                execute_if(was_used | do_alloc, [&](){
                    gctr[w][j].write(gindex[j], new_ctr, extra_cycle);
                });
```

3. val无法再次赋值
