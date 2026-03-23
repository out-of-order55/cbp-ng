# Plan: Remove WB drain threshold T=2

## Goal
Disable the T=2 drain accumulation behavior for the active FGEHL write buffer.

## Scope
- Minimal functional change in `predictors/my_bp_v1.hpp` only.
- Change FGEHL WB instantiation from `T=2` to `T=1`.

## Steps
1. Update FGEHL declaration to use `wb_ram<...,4,1>`.
2. Build with:
   `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE`
3. Run fixed trace:
   `./my_bp /nfs/home/guobing1/cbp-ng/trace/cbp-ng_training_traces/502-gcc-all_16112_trace.gz test 1000000 2000000 --format human`
4. Report metrics and files changed.
