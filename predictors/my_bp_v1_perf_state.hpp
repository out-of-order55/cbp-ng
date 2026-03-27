#pragma once

#include <cstdint>

#include "my_bp_v1_perf_events.hpp"

namespace my_bp_v1_perf {

struct PerfState {
    std::uint64_t predictions = 0;
    std::uint64_t correct = 0;
    std::uint64_t mispredictions = 0;
    std::uint64_t extra_cycle_total = 0;
    std::uint64_t extra_cycle_badpred = 0;
    std::uint64_t extra_cycle_mispredict = 0;
    std::uint64_t extra_cycle_p1_update = 0;
    std::uint64_t extra_cycle_sc_update = 0;
    std::uint64_t gate_count = 0;
    std::uint64_t gate_mispred = 0;
    std::uint64_t gate_update_inc = 0;
    std::uint64_t gate_update_dec = 0;
    std::uint64_t bgehl_read_ops = 0;
    std::uint64_t bgehl_write_ops = 0;
    std::uint64_t bgehl_branch_slots = 0;
    std::uint64_t bgehl_filter_allow = 0;
    std::uint64_t bgehl_filter_block = 0;
    std::uint64_t bgehl_nonzero_contrib = 0;
    std::uint64_t bgehl_update_candidate = 0;
    std::uint64_t bgehl_update_allow = 0;
    std::uint64_t bgehl_update_block = 0;

    void on_resolve(const ResolveEvent &event)
    {
        if (!event.is_branch) {
            return;
        }
        ++predictions;
        if (event.correct) {
            ++correct;
        } else {
            ++mispredictions;
        }
    }

    void on_extra_cycle(const ExtraCycleEvent &event)
    {
        extra_cycle_total += static_cast<std::uint64_t>(event.extra_cycle);
        extra_cycle_badpred += static_cast<std::uint64_t>(event.badpred);
        extra_cycle_mispredict += static_cast<std::uint64_t>(event.mispredict);
        extra_cycle_p1_update += static_cast<std::uint64_t>(event.p1_update);
        extra_cycle_sc_update += static_cast<std::uint64_t>(event.sc_update);
    }

    void on_gate(const GateEvent &event)
    {
        gate_count += static_cast<std::uint64_t>(event.gating);
        gate_mispred += static_cast<std::uint64_t>(event.gating && event.mispredict);
        gate_update_inc += static_cast<std::uint64_t>(event.gate_inc);
        gate_update_dec += static_cast<std::uint64_t>(event.gate_dec);
    }

    void on_bgehl_read_ops(std::uint64_t read_ops)
    {
        bgehl_read_ops += read_ops;
    }

    void on_bgehl_write_ops(std::uint64_t write_ops)
    {
        bgehl_write_ops += write_ops;
    }

    void on_bgehl_slot(const BgehlSlotEvent &event)
    {
        if (!event.is_branch) {
            return;
        }
        ++bgehl_branch_slots;
        if (event.filter_open) {
            ++bgehl_filter_allow;
            bgehl_nonzero_contrib += static_cast<std::uint64_t>(event.nonzero_contrib);
        } else {
            ++bgehl_filter_block;
        }
        bgehl_update_candidate += static_cast<std::uint64_t>(event.update_candidate);
        bgehl_update_allow += static_cast<std::uint64_t>(event.update_allow);
        bgehl_update_block += static_cast<std::uint64_t>(event.update_block);
    }
};

} // namespace my_bp_v1_perf
