# Plan: WB Request Lifetime Distribution (enqueue -> RAM write)

## Goal
Measure how many cycles a buffered write request stays in WB before being drained to RAM.

## Implementation
1. Add PERF-only per-entry enqueue timestamp tracking in `wb_ram`.
2. On dequeue, compute request lifetime in cycles and update a histogram.
3. Expose histogram getters from `wb_ram`.
4. Aggregate and print FGEHL WB lifetime histogram in `my_bp_v1::print_perf_counters()`.

## Validation
1. Compile with:
   `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE -DPERF_COUNTERS -DCHEATING_MODE -DFREE_FANOUT`
2. Run fixed trace:
   `./my_bp /nfs/home/guobing1/cbp-ng/trace/cbp-ng_training_traces/502-gcc-all_16112_trace.gz test 1000000 2000000 --format human`
3. Report lifetime distribution and key summary statistics.
