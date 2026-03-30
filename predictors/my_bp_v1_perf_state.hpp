#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace my_bp_v1_perf {

template<std::size_t NUMG>
struct PerfState {
    static_assert(NUMG % 2 == 0, "my_bp_v1 perf state expects even NUMG");
    static constexpr std::size_t NUMP = NUMG / 2;

    std::uint64_t perf_predictions = 0;
    std::uint64_t perf_correct = 0;
    std::array<std::uint64_t, NUMG> perf_provider_used = {};
    std::array<std::uint64_t, NUMG> perf_provider_correct = {};
    std::array<std::uint64_t, NUMG> perf_provider_wrong = {};
    std::array<std::uint64_t, NUMG> perf_alt_used = {};
    std::array<std::uint64_t, NUMG> perf_alt_correct = {};
    std::array<std::uint64_t, NUMG> perf_alt_wrong = {};
    std::uint64_t perf_bimodal_used = 0;
    std::uint64_t perf_bimodal_correct = 0;
    std::uint64_t perf_bimodal_wrong = 0;
    std::uint64_t perf_mispred_blame_tage = 0;
    std::uint64_t perf_mispred_blame_sc = 0;
    std::uint64_t perf_mispred_blame_p1 = 0;
    std::array<std::uint64_t, NUMG> perf_table_reads = {};
    std::array<std::uint64_t, NUMG> perf_table_hits = {};
    std::array<std::uint64_t, NUMG> perf_table_alloc = {};
    std::uint64_t perf_alloc_failures = 0;
    std::uint64_t perf_alloc_fail_highest = 0;
    std::uint64_t perf_alloc_fail_noubit = 0;
    std::uint64_t perf_extra_cycle_total = 0;
    std::uint64_t perf_extra_cycle_badpred = 0;
    std::uint64_t perf_extra_cycle_mispredict = 0;
    std::uint64_t perf_extra_cycle_p1_update = 0;
    std::uint64_t perf_extra_cycle_sc_update = 0;
    std::uint64_t perf_sc_override = 0;
    std::uint64_t perf_sc_override_correct = 0;
    std::uint64_t perf_sc_use = 0;
    std::uint64_t perf_sc_use_correct = 0;
    std::uint64_t perf_sc_use_taken = 0;
    std::uint64_t perf_sc_use_nottaken = 0;
    std::uint64_t perf_sc_use_same_as_tage = 0;
    std::uint64_t perf_sc_use_flip_tage = 0;
    std::uint64_t perf_sc_use_weak = 0;
    std::uint64_t perf_sc_use_mid = 0;
    std::uint64_t perf_sc_use_sat = 0;
    std::uint64_t perf_sc_stage_prov_hit = 0;
    std::uint64_t perf_sc_stage_do_update = 0;
    std::uint64_t perf_sc_stage_candidate = 0;
    std::uint64_t perf_sc_stage_guard_pass = 0;
    std::uint64_t perf_sc_skip_no_provider = 0;
    std::uint64_t perf_sc_skip_no_do_update = 0;
    std::uint64_t perf_sc_skip_guard = 0;
    std::uint64_t perf_global_thre_update = 0;
    std::uint64_t perf_global_thre_inc = 0;
    std::uint64_t perf_global_thre_dec = 0;
    std::uint64_t perf_thre_update = 0;
    std::uint64_t perf_thre_update_inc = 0;
    std::uint64_t perf_thre_update_dec = 0;
    std::uint64_t perf_mispred_sc_not_used = 0;
    std::uint64_t perf_mispred_sc_keep = 0;
    std::uint64_t perf_mispred_sc_flip = 0;
    std::uint64_t perf_mispred_sc_flip_harmful = 0;
    std::uint64_t perf_mispred_sc_flip_both_wrong = 0;
    std::uint64_t perf_bgehl_read_ops = 0;
    std::uint64_t perf_bgehl_write_ops = 0;
    std::uint64_t perf_bgehl_branch_slots = 0;
    std::uint64_t perf_bgehl_filter_allow = 0;
    std::uint64_t perf_bgehl_filter_block = 0;
    std::uint64_t perf_bgehl_nonzero_contrib = 0;
    std::uint64_t perf_bgehl_update_candidate = 0;
    std::uint64_t perf_bgehl_update_allow = 0;
    std::uint64_t perf_bgehl_update_block = 0;
    std::array<std::uint64_t, 3> perf_share_route_zero = {};
    std::array<std::uint64_t, 3> perf_share_route_one = {};
    std::array<std::uint64_t, NUMP> perf_group_idx_mismatch = {};
    std::uint64_t perf_group_idx_samples = 0;
    std::array<std::array<std::uint64_t, NUMP>, NUMG> perf_logical_tag_hit_phys = {};
    std::array<std::array<std::uint64_t, NUMP>, NUMG> perf_logical_sel_phys = {};
    std::array<std::uint64_t, NUMP> perf_phys_slot0_hit = {};
    std::array<std::uint64_t, NUMP> perf_phys_slot1_hit = {};
    std::array<std::uint64_t, NUMP> perf_phys_any_hit = {};
    std::array<std::uint64_t, NUMP> perf_phys_dual_slot_hit = {};
    std::array<std::uint64_t, NUMP> perf_phys_tag_w = {};
    std::array<std::uint64_t, NUMP> perf_phys_pred_w = {};
    std::array<std::uint64_t, NUMP> perf_phys_hyst_w = {};
    std::array<std::uint64_t, NUMP> perf_phys_u_w = {};
    std::array<std::array<std::uint64_t, 4>, NUMG> perf_conf = {};

    // P1 micro-TAGE counters (P1 correctness is defined against P2 direction).
    std::uint64_t perf_mt_p1_slots = 0;
    std::uint64_t perf_mt_p1_match_p2 = 0;
    std::uint64_t perf_mt_p1_disagree_p2 = 0;
    std::uint64_t perf_mt_p1_from_micro = 0;
    std::uint64_t perf_mt_p1_from_gshare = 0;
    std::uint64_t perf_mt_p1_micro_match_p2 = 0;
    std::uint64_t perf_mt_p1_micro_disagree_p2 = 0;
    std::uint64_t perf_mt_p1_gshare_match_p2 = 0;
    std::uint64_t perf_mt_p1_gshare_disagree_p2 = 0;
    std::uint64_t perf_mt_provider_hit = 0;
    std::array<std::uint64_t, 4> perf_mt_provider_hit_table = {};
    std::array<std::uint64_t, 4> perf_mt_provider_use_table = {};
    std::array<std::uint64_t, 4> perf_mt_provider_correct_table = {};
    std::array<std::uint64_t, 4> perf_mt_provider_wrong_table = {};
    std::array<std::uint64_t, 4> perf_mt_table_reads = {};
    std::array<std::uint64_t, 4> perf_mt_table_hits = {};
    std::array<std::uint64_t, 4> perf_mt_table_alloc = {};
    std::array<std::array<std::uint64_t, 4>, 4> perf_mt_conf = {};
    std::uint64_t perf_mt_provider_wrong = 0;
    std::uint64_t perf_mt_uclr_ctr_gate = 0;
    std::uint64_t perf_mt_uclr_ctr_inc = 0;
    std::uint64_t perf_mt_uclr_ctr_dec = 0;
    std::uint64_t perf_mt_uclr_sat = 0;
    std::uint64_t perf_mt_u_reset = 0;
    std::array<std::uint64_t, 4> perf_mt_alloc_pick_table = {};
    std::array<std::uint64_t, 4> perf_mt_alloc_req_table = {};
    std::array<std::uint64_t, 4> perf_mt_pred_req_table = {};
    std::array<std::uint64_t, 4> perf_mt_hyst_req_table = {};
    std::array<std::uint64_t, 4> perf_mt_u_req_table = {};
    std::array<std::uint64_t, 4> perf_mt_uclear_req_table = {};
    std::array<std::uint64_t, 4> perf_mt_tag_w_table = {};
    std::array<std::uint64_t, 4> perf_mt_pred_w_table = {};
    std::array<std::uint64_t, 4> perf_mt_hyst_w_table = {};
    std::array<std::uint64_t, 4> perf_mt_u_w_table = {};
};

} // namespace my_bp_v1_perf
