#pragma once

#include <cstdint>

namespace my_bp_v1_perf {

struct ResolveEvent {
    std::uint64_t cycle = 0;
    std::uint8_t offset = 0;
    bool is_branch = false;
    bool correct = false;
};

struct ExtraCycleEvent {
    bool extra_cycle = false;
    bool badpred = false;
    bool mispredict = false;
    bool p1_update = false;
    bool sc_update = false;
};

struct GateEvent {
    bool gating = false;
    bool mispredict = false;
    bool gate_inc = false;
    bool gate_dec = false;
};

struct BgehlSlotEvent {
    bool is_branch = false;
    bool filter_open = false;
    bool nonzero_contrib = false;
    bool update_candidate = false;
    bool update_allow = false;
    bool update_block = false;
};

} // namespace my_bp_v1_perf
