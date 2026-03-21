# Mispred Reason Counter Plan (2026-03-21)

## Objective
Add performance counters that explain which predictor path is responsible for mispredictions.

## Scope
- File: `predictors/my_bp_v1.hpp`
- Keep existing predictor behavior unchanged.
- Add only counter bookkeeping and report printing.

## Counter Design
1. Total mispredictions: `perf_predictions - perf_correct`.
2. Source-level wrong distribution (by selected source before SC override logic):
   - Provider wrong (per-table and total)
   - Alternate wrong (per-table and total)
   - Bimodal wrong
3. Reason-level blame (exclusive partition of final wrong outcomes):
   - `P1` blame (when `GATE` selects P1 and final prediction is wrong)
   - `SC` blame (when `MY_SC` flips a correct TAGE decision into a wrong final decision)
   - `TAGE` blame (all remaining wrong outcomes)
4. SC mispred detail breakdown (diagnostic, non-exclusive to source table):
   - Wrong with `use_sc=0`
   - Wrong with `use_sc=1` and keep TAGE direction
   - Wrong with `use_sc=1` and flip direction
   - Within wrong flips: `SC harmful flip` (TAGE correct, SC wrong) vs `both wrong`.

## Validation
1. Build:
   - `./compile my_bp -DPREDICTOR="my_bp_v1<>" -DVERBOSE`
2. Run:
   - `./my_bp /nfs/home/guobing1/cbp-ng/trace/cbp-ng_training_traces/502-gcc-all_16112_trace.gz test 1000000 2000000 --format human`
3. Confirm report appears and counters are internally consistent.
