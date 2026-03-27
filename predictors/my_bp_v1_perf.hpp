#pragma once

#include "my_bp_v1_perf_events.hpp"
#include "my_bp_v1_perf_report.hpp"
#include "my_bp_v1_perf_state.hpp"

namespace my_bp_v1_perf {

inline void on_resolve(PerfState &state, const ResolveEvent &event)
{
    state.on_resolve(event);
}

inline void on_extra_cycle(PerfState &state, const ExtraCycleEvent &event)
{
    state.on_extra_cycle(event);
}

inline void on_gate(PerfState &state, const GateEvent &event)
{
    state.on_gate(event);
}

inline void on_bgehl_read_ops(PerfState &state, std::uint64_t read_ops)
{
    state.on_bgehl_read_ops(read_ops);
}

inline void on_bgehl_write_ops(PerfState &state, std::uint64_t write_ops)
{
    state.on_bgehl_write_ops(write_ops);
}

inline void on_bgehl_slot(PerfState &state, const BgehlSlotEvent &event)
{
    state.on_bgehl_slot(event);
}

} // namespace my_bp_v1_perf
