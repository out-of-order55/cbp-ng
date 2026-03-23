# Plan: Run one PERF_COUNTERS validation

## Goal
Run one counter-mode build and trace execution for current `my_bp_v1` state.

## Steps
1. Compile with counter flags:
   `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE -DPERF_COUNTERS -DCHEATING_MODE -DFREE_FANOUT`
2. Run fixed trace:
   `./my_bp /nfs/home/guobing1/cbp-ng/trace/cbp-ng_training_traces/502-gcc-all_16112_trace.gz test 1000000 2000000 --format human`
3. Save full output and extract key counters.
