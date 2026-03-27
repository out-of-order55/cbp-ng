#pragma once

#include <cstdint>
#include <ostream>

#include "my_bp_v1_perf_state.hpp"

namespace my_bp_v1_perf {

inline void print_phase_c1_summary(
    std::ostream &os,
    const PerfState &state,
    std::uint64_t legacy_predictions,
    std::uint64_t legacy_correct,
    std::uint64_t legacy_mispredictions)
{
    os << "\n┌─ Scheme C Phase C1 (Event API Baseline) ───────────────────────┐\n";
    os << "│ Event Predictions/Correct/Mispred: "
       << state.predictions << " / " << state.correct << " / " << state.mispredictions << "│\n";
    os << "│ Legacy Predictions/Correct/Mispred: "
       << legacy_predictions << " / " << legacy_correct << " / " << legacy_mispredictions << "│\n";
    const bool pass = (state.predictions == legacy_predictions) &&
                      (state.correct == legacy_correct) &&
                      (state.mispredictions == legacy_mispredictions);
    os << "│ Baseline parity check: " << (pass ? "PASS" : "FAIL") << "                              │\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";
}

inline void print_phase_c2_ec_gate_bgehl_summary(
    std::ostream &os,
    const PerfState &state,
    std::uint64_t legacy_extra_cycle_total,
    std::uint64_t legacy_extra_cycle_badpred,
    std::uint64_t legacy_extra_cycle_mispredict,
    std::uint64_t legacy_extra_cycle_p1_update,
    std::uint64_t legacy_extra_cycle_sc_update,
    std::uint64_t legacy_gate_count,
    std::uint64_t legacy_gate_mispred,
    std::uint64_t legacy_gate_update_inc,
    std::uint64_t legacy_gate_update_dec,
    std::uint64_t legacy_bgehl_read_ops,
    std::uint64_t legacy_bgehl_write_ops,
    std::uint64_t legacy_bgehl_branch_slots,
    std::uint64_t legacy_bgehl_filter_allow,
    std::uint64_t legacy_bgehl_filter_block,
    std::uint64_t legacy_bgehl_nonzero_contrib,
    std::uint64_t legacy_bgehl_update_candidate,
    std::uint64_t legacy_bgehl_update_allow,
    std::uint64_t legacy_bgehl_update_block)
{
    const bool extra_cycle_pass =
        state.extra_cycle_total == legacy_extra_cycle_total &&
        state.extra_cycle_badpred == legacy_extra_cycle_badpred &&
        state.extra_cycle_mispredict == legacy_extra_cycle_mispredict &&
        state.extra_cycle_p1_update == legacy_extra_cycle_p1_update &&
        state.extra_cycle_sc_update == legacy_extra_cycle_sc_update;

    const bool gate_pass =
        state.gate_count == legacy_gate_count &&
        state.gate_mispred == legacy_gate_mispred &&
        state.gate_update_inc == legacy_gate_update_inc &&
        state.gate_update_dec == legacy_gate_update_dec;

    const bool bgehl_pass =
        state.bgehl_read_ops == legacy_bgehl_read_ops &&
        state.bgehl_write_ops == legacy_bgehl_write_ops &&
        state.bgehl_branch_slots == legacy_bgehl_branch_slots &&
        state.bgehl_filter_allow == legacy_bgehl_filter_allow &&
        state.bgehl_filter_block == legacy_bgehl_filter_block &&
        state.bgehl_nonzero_contrib == legacy_bgehl_nonzero_contrib &&
        state.bgehl_update_candidate == legacy_bgehl_update_candidate &&
        state.bgehl_update_allow == legacy_bgehl_update_allow &&
        state.bgehl_update_block == legacy_bgehl_update_block;

    os << "\n┌─ Scheme C Phase C2 (ExtraCycle/Gate/BGEHL) ───────────────────┐\n";
    os << "│ ExtraCycle parity: " << (extra_cycle_pass ? "PASS" : "FAIL") << "                               │\n";
    os << "│ Gate parity:       " << (gate_pass ? "PASS" : "FAIL") << "                               │\n";
    os << "│ BGEHL parity:      " << (bgehl_pass ? "PASS" : "FAIL") << "                               │\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";
}

} // namespace my_bp_v1_perf
