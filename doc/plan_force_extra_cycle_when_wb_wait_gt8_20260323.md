# Plan: Force extra cycle when WB wait exceeds 8 cycles

## Goal
If a buffered write request cannot drain for more than 8 cycles, force `need_extra_cycle` and allow WB drain.

## Implementation
1. In `wb_ram`, track head wait cycles with a small saturating counter.
2. Expose `drain_urgent()` signal: queue non-empty and head-wait > 8.
3. In `my_bp_v1::update_cycle`, OR `drain_urgent()` (from FGEHL WB instances) into `extra_cycle`.

## Validation
1. Compile functional build:
   `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE`
2. Run fixed trace and save output.
3. Compile counter build:
   `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE -DPERF_COUNTERS -DCHEATING_MODE -DFREE_FANOUT`
4. Run fixed trace and check WB lifetime distribution.
