#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ostream>
#include <string>

#include "my_bp_v1_perf_state.hpp"

namespace my_bp_v1_perf {

template<std::size_t NUMG>
struct PerfReportContext {
    std::array<std::uint64_t, NUMG> hist_len = {};
    const char *ghyst_impl_name = "";
    std::uint64_t ghyst_reads = 0;
    std::uint64_t ghyst_writes = 0;
    std::uint64_t ghyst_stale_reads = 0;
    std::uint64_t ghyst_drop_writes = 0;
    std::uint64_t ghyst_write_hit_events = 0;
};

struct PerfReportTotals {
    std::uint64_t total_all = 0;
    std::uint64_t total_mispred = 0;
};

template<std::size_t NUMG>
inline PerfReportTotals print_main_perf_report(
    std::ostream &os,
    const PerfState<NUMG> &state,
    const PerfReportContext<NUMG> &ctx)
{
    os << "\n╔════════════════════════════════════════════════════════════════╗\n";
    os << "║           TAGE PREDICTOR PERFORMANCE COUNTERS                   ║\n";
    os << "╚════════════════════════════════════════════════════════════════╝\n";

    os << "\n┌─ Overall Statistics ────────────────────────────────────────────┐\n";
    os << "│ Total Predictions: " << std::setw(43) << std::left << state.perf_predictions << "│\n";
    os << "│ Correct:           " << std::setw(43) << std::left << state.perf_correct << "│\n";
    if (state.perf_predictions > 0) {
        double accuracy = (100.0 * state.perf_correct) / state.perf_predictions;
        os << "│ Accuracy:          " << std::fixed << std::setprecision(2)
           << std::setw(40) << std::left << accuracy << "% │\n";
    }
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ TAGE Table Statistics ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
    os << "│ Tbl │ HistLen │ Reads │ Hits │ Hit%  │ Prov │ PrvAcc% │ Alt │ AltAcc% │ Total │ TotAcc% │ Alloc │ C0(weak) │ C1(mwk) │ C2(mst) │ C3(str) │\n";
    os << "├─────┼─────────┼───────┼──────┼───────┼──────┼─────────┼─────┼─────────┼───────┼─────────┼───────┼──────────┼─────────┼─────────┼─────────┤\n";

    for (std::size_t j = 0; j < NUMG; j++) {
        std::uint64_t reads = state.perf_table_reads[j];
        std::uint64_t hits = state.perf_table_hits[j];
        std::uint64_t prov = state.perf_provider_used[j];
        std::uint64_t prov_ok = state.perf_provider_correct[j];
        std::uint64_t alt = state.perf_alt_used[j];
        std::uint64_t alt_ok = state.perf_alt_correct[j];
        std::uint64_t total = prov + alt;
        std::uint64_t total_ok = prov_ok + alt_ok;
        std::uint64_t alloc = state.perf_table_alloc[j];
        std::uint64_t c0 = state.perf_conf[j][0], c1 = state.perf_conf[j][1], c2 = state.perf_conf[j][2], c3 = state.perf_conf[j][3];
        std::uint64_t ct = c0 + c1 + c2 + c3;

        os << "│ " << std::setw(3) << std::left << j << " │ ";
        os << std::setw(7) << std::right << ctx.hist_len[j] << " │ ";
        os << std::setw(5) << std::right << reads << " │ ";
        os << std::setw(4) << std::right << hits << " │ ";
        if (reads > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(5) << std::right << (100.0 * hits / reads) << "% │ ";
        } else {
            os << "  N/A │ ";
        }
        os << std::setw(4) << std::right << prov << " │ ";
        if (prov > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0 * prov_ok / prov) << "% │ ";
        } else {
            os << "    N/A │ ";
        }
        os << std::setw(3) << std::right << alt << " │ ";
        if (alt > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0 * alt_ok / alt) << "% │ ";
        } else {
            os << "    N/A │ ";
        }
        os << std::setw(5) << std::right << total << " │ ";
        if (total > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0 * total_ok / total) << "% │ ";
        } else {
            os << "    N/A │ ";
        }
        os << std::setw(5) << std::right << alloc << " │ ";
        if (ct > 0) {
            os << std::setw(5) << c0 << "(" << std::fixed << std::setprecision(0) << std::setw(2) << (100.0 * c0 / ct) << "%) │ ";
            os << std::setw(5) << c1 << "(" << std::setw(2) << (100.0 * c1 / ct) << "%) │ ";
            os << std::setw(5) << c2 << "(" << std::setw(2) << (100.0 * c2 / ct) << "%) │ ";
            os << std::setw(5) << c3 << "(" << std::setw(2) << (100.0 * c3 / ct) << "%) │\n";
        } else {
            os << "     N/A │      N/A │      N/A │      N/A │\n";
        }
    }
    os << "└─────┴─────────┴───────┴──────┴───────┴──────┴─────────┴─────┴─────────┴───────┴─────────┴───────┴──────────┴─────────┴─────────┴─────────┘\n";

    os << "\n┌─ Prediction Source Distribution ────────────────────────────────┐\n";
    os << "│ Source   │ Count │ Correct │ Accuracy │ % of Total │\n";
    os << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";

    std::uint64_t total_prov = 0, total_prov_ok = 0;
    std::uint64_t total_alt = 0, total_alt_ok = 0;
    std::uint64_t total_prov_wrong = 0, total_alt_wrong = 0;
    for (std::size_t j = 0; j < NUMG; j++) {
        total_prov += state.perf_provider_used[j];
        total_prov_ok += state.perf_provider_correct[j];
        total_prov_wrong += state.perf_provider_wrong[j];
        total_alt += state.perf_alt_used[j];
        total_alt_ok += state.perf_alt_correct[j];
        total_alt_wrong += state.perf_alt_wrong[j];
    }

    std::uint64_t total_all = total_prov + total_alt + state.perf_bimodal_used;
    std::uint64_t total_mispred = state.perf_predictions - state.perf_correct;

    os << "│ Provider │ " << std::setw(5) << std::right << total_prov << " │ ";
    os << std::setw(7) << std::right << total_prov_ok << " │ ";
    if (total_prov > 0) {
        os << std::fixed << std::setprecision(2) << std::setw(8) << std::right << (100.0 * total_prov_ok / total_prov) << "% │ ";
    } else {
        os << "    N/A │ ";
    }
    if (total_all > 0) {
        os << std::fixed << std::setprecision(1) << std::setw(9) << std::right << (100.0 * total_prov / total_all) << "% │\n";
    } else {
        os << "    N/A │\n";
    }

    os << "│ Alternate│ " << std::setw(5) << std::right << total_alt << " │ ";
    os << std::setw(7) << std::right << total_alt_ok << " │ ";
    if (total_alt > 0) {
        os << std::fixed << std::setprecision(2) << std::setw(8) << std::right << (100.0 * total_alt_ok / total_alt) << "% │ ";
    } else {
        os << "    N/A │ ";
    }
    if (total_all > 0) {
        os << std::fixed << std::setprecision(1) << std::setw(9) << std::right << (100.0 * total_alt / total_all) << "% │\n";
    } else {
        os << "    N/A │\n";
    }

    os << "│ Bimodal  │ " << std::setw(5) << std::right << state.perf_bimodal_used << " │ ";
    os << std::setw(7) << std::right << state.perf_bimodal_correct << " │ ";
    if (state.perf_bimodal_used > 0) {
        os << std::fixed << std::setprecision(2) << std::setw(8) << std::right
           << (100.0 * state.perf_bimodal_correct / state.perf_bimodal_used) << "% │ ";
    } else {
        os << "    N/A │ ";
    }
    if (total_all > 0) {
        os << std::fixed << std::setprecision(1) << std::setw(9) << std::right << (100.0 * state.perf_bimodal_used / total_all) << "% │\n";
    } else {
        os << "    N/A │\n";
    }

#ifdef MY_SC
    os << "│ SC Overr │ " << std::setw(5) << std::right << state.perf_sc_override << " │ ";
    os << std::setw(7) << std::right << state.perf_sc_override_correct << " │ ";
    if (state.perf_sc_override > 0) {
        os << std::fixed << std::setprecision(2) << std::setw(8) << std::right
           << (100.0 * state.perf_sc_override_correct / state.perf_sc_override) << "% │ ";
    } else {
        os << "    N/A │ ";
    }
    if (total_all > 0) {
        os << std::fixed << std::setprecision(1) << std::setw(9) << std::right << (100.0 * state.perf_sc_override / total_all) << "% │\n";
    } else {
        os << "    N/A │\n";
    }
    os << "│ ThreUpd  │ " << std::setw(5) << std::right << state.perf_thre_update << " │ ";
    os << std::setw(7) << std::right << state.perf_thre_update_inc << " inc │ ";
    os << std::setw(7) << std::right << state.perf_thre_update_dec << " dec │\n";
#endif

    os << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";
    os << "│ Total    │ " << std::setw(5) << std::right << total_all << " │ ";
    std::uint64_t total_correct = total_prov_ok + total_alt_ok + state.perf_bimodal_correct;
    os << std::setw(7) << std::right << total_correct << " │ ";
    if (total_all > 0) {
        os << std::fixed << std::setprecision(2) << std::setw(8) << std::right << (100.0 * total_correct / total_all) << "% │ ";
        os << std::fixed << std::setprecision(1) << std::setw(9) << std::right << 100.0 << "% │\n";
    } else {
        os << "    N/A │     N/A │\n";
    }
    os << "└──────────┴───────┴─────────┴──────────┴────────────┘\n";

#ifdef MY_SC
    os << "\n┌─ SC Usage Statistics ────────────────────────────────────────────┐\n";
    os << "│ SC Use Count:           " << std::setw(36) << std::left << state.perf_sc_use << "│\n";
    os << "│ SC Use Correct:         " << std::setw(36) << std::left << state.perf_sc_use_correct << "│\n";
    if (state.perf_sc_use > 0) {
        os << "│ SC Use Accuracy:        " << std::setw(35) << std::left
           << (std::to_string((100.0 * state.perf_sc_use_correct / state.perf_sc_use))).substr(0, 5) + "%" << "│\n";
    } else {
        os << "│ SC Use Accuracy:        " << std::setw(36) << std::left << "N/A" << "│\n";
    }
    os << "│ SC Use Taken/NotTaken:  " << std::setw(36) << std::left
       << (std::to_string(state.perf_sc_use_taken) + " / " + std::to_string(state.perf_sc_use_nottaken)) << "│\n";
    os << "│ SC Keep/Flip TAGE:      " << std::setw(36) << std::left
       << (std::to_string(state.perf_sc_use_same_as_tage) + " / " + std::to_string(state.perf_sc_use_flip_tage)) << "│\n";
    os << "│ SC Use Weak/Mid/Sat:    " << std::setw(36) << std::left
       << (std::to_string(state.perf_sc_use_weak) + " / " + std::to_string(state.perf_sc_use_mid) + " / " + std::to_string(state.perf_sc_use_sat)) << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ SC Update Pipeline ─────────────────────────────────────────────┐\n";
    os << "│ Provider Hit:           " << std::setw(36) << std::left << state.perf_sc_stage_prov_hit << "│\n";
    os << "│ do_update Passed:       " << std::setw(36) << std::left << state.perf_sc_stage_do_update << "│\n";
    os << "│ Candidate (hit&update): " << std::setw(36) << std::left << state.perf_sc_stage_candidate << "│\n";
    os << "│ Guard Pass:             " << std::setw(36) << std::left << state.perf_sc_stage_guard_pass << "│\n";
    os << "│ Skip No Provider:       " << std::setw(36) << std::left << state.perf_sc_skip_no_provider << "│\n";
    os << "│ Skip No do_update:      " << std::setw(36) << std::left << state.perf_sc_skip_no_do_update << "│\n";
    os << "│ Skip Guard:             " << std::setw(36) << std::left << state.perf_sc_skip_guard << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ SC Threshold Updates ───────────────────────────────────────────┐\n";
    os << "│ Local Threshold (tot/inc/dec): " << std::setw(29) << std::left
       << (std::to_string(state.perf_thre_update) + " / " + std::to_string(state.perf_thre_update_inc) + " / " + std::to_string(state.perf_thre_update_dec)) << "│\n";
    os << "│ Global Threshold (tot/inc/dec): " << std::setw(28) << std::left
       << (std::to_string(state.perf_global_thre_update) + " / " + std::to_string(state.perf_global_thre_inc) + " / " + std::to_string(state.perf_global_thre_dec)) << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";
#endif

    os << "\n┌─ Mispred Reason Statistics ──────────────────────────────────────┐\n";
    os << "│ Total Mispred:           " << std::setw(36) << std::left << total_mispred << "│\n";
    os << "│ Src Wrong (Prov/Alt/Bim): " << std::setw(35) << std::left
       << (std::to_string(total_prov_wrong) + " / " + std::to_string(total_alt_wrong) + " / " + std::to_string(state.perf_bimodal_wrong)) << "│\n";
    os << "│ Blame TAGE/SC/P1:        " << std::setw(36) << std::left
       << (std::to_string(state.perf_mispred_blame_tage) + " / " + std::to_string(state.perf_mispred_blame_sc) + " / " + std::to_string(state.perf_mispred_blame_p1)) << "│\n";
#ifdef MY_SC
    os << "│ SC Wrong use_sc=0:       " << std::setw(36) << std::left << state.perf_mispred_sc_not_used << "│\n";
    os << "│ SC Wrong Keep/Flip:      " << std::setw(36) << std::left
       << (std::to_string(state.perf_mispred_sc_keep) + " / " + std::to_string(state.perf_mispred_sc_flip)) << "│\n";
    os << "│ SC Flip Harm/BothWrong:  " << std::setw(36) << std::left
       << (std::to_string(state.perf_mispred_sc_flip_harmful) + " / " + std::to_string(state.perf_mispred_sc_flip_both_wrong)) << "│\n";
#endif
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ Extra Cycle Statistics ──────────────────────────────────────┐\n";
    os << "│ Total Extra Cycles:     " << std::setw(36) << std::left << state.perf_extra_cycle_total << "│\n";
    os << "│   Weak & Wrong (TAGE):  " << std::setw(36) << std::left << state.perf_extra_cycle_badpred << "│\n";
    os << "│   Misprediction:        " << std::setw(36) << std::left << state.perf_extra_cycle_mispredict << "│\n";
    os << "│   P1 Update:            " << std::setw(36) << std::left << state.perf_extra_cycle_p1_update << "│\n";
#ifdef MY_SC
    os << "│   SC Update:            " << std::setw(36) << std::left << state.perf_extra_cycle_sc_update << "│\n";
#endif
    os << "└─────────────────────────────────────────────────────────────────┘\n";

#ifdef MY_SC
#ifdef SC_BGEHL
    os << "\n┌─ BGEHL Statistics ──────────────────────────────────────────────┐\n";
    os << "│ ReadOps/WriteOps:        " << std::setw(36) << std::left
       << (std::to_string(state.perf_bgehl_read_ops) + " / " + std::to_string(state.perf_bgehl_write_ops)) << "│\n";
    os << "│ Branch Slots:            " << std::setw(36) << std::left << state.perf_bgehl_branch_slots << "│\n";
    os << "│ Filter Allow/Block:      " << std::setw(36) << std::left
       << (std::to_string(state.perf_bgehl_filter_allow) + " / " + std::to_string(state.perf_bgehl_filter_block)) << "│\n";
    os << "│ Nonzero Contribution:    " << std::setw(36) << std::left
       << (std::to_string(state.perf_bgehl_nonzero_contrib) +
           (state.perf_bgehl_branch_slots ? " (" + std::to_string(100.0 * state.perf_bgehl_nonzero_contrib / state.perf_bgehl_branch_slots).substr(0, 6) + "%)" : "")) << "│\n";
    os << "│ Update Cand/Allow/Block: " << std::setw(36) << std::left
       << (std::to_string(state.perf_bgehl_update_candidate) + " / " +
           std::to_string(state.perf_bgehl_update_allow) + " / " +
           std::to_string(state.perf_bgehl_update_block)) << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";
#endif
#endif

    os << "\n┌─ GHyst RAM Statistics ──────────────────────────────────────────┐\n";
    os << "│ Impl:                   " << std::setw(36) << std::left << ctx.ghyst_impl_name << "│\n";
    os << "│ Reads/Writes:           " << std::setw(36) << std::left
       << (std::to_string(ctx.ghyst_reads) + " / " + std::to_string(ctx.ghyst_writes)) << "│\n";
    os << "│ Stale Reads:            " << std::setw(36) << std::left
       << (std::to_string(ctx.ghyst_stale_reads) + (ctx.ghyst_reads ? " (" + std::to_string(100.0 * ctx.ghyst_stale_reads / ctx.ghyst_reads).substr(0, 6) + "%)" : "")) << "│\n";
    os << "│ Dropped Writes:         " << std::setw(36) << std::left
       << (std::to_string(ctx.ghyst_drop_writes) + (ctx.ghyst_writes ? " (" + std::to_string(100.0 * ctx.ghyst_drop_writes / ctx.ghyst_writes).substr(0, 6) + "%)" : "")) << "│\n";
    os << "│ Write Hit Events:       " << std::setw(36) << std::left
       << (std::to_string(ctx.ghyst_write_hit_events) + (ctx.ghyst_writes ? " (" + std::to_string(100.0 * ctx.ghyst_write_hit_events / ctx.ghyst_writes).substr(0, 6) + "%)" : "")) << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ Allocation Statistics ─────────────────────────────────────────┐\n";
    os << "│ Allocation Failures: " << std::setw(41) << std::left << state.perf_alloc_failures << "│\n";
    os << "│   - Already at highest table: " << std::setw(32) << std::left << state.perf_alloc_fail_highest << "│\n";
    os << "│   - No ubit=0 victim found: " << std::setw(34) << std::left << state.perf_alloc_fail_noubit << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    constexpr std::size_t NUMP = PerfState<NUMG>::NUMP;
    os << "\n┌─ Physical Table Write Load ─────────────────────────────────────┐\n";
    os << "│ Phys │ TagW │ PredW │ HystW │ UW │\n";
    os << "├──────┼──────┼───────┼───────┼────┤\n";
    for (std::size_t p = 0; p < NUMP; p++) {
        os << "│ p" << p << "   │ "
           << std::setw(4) << std::right << state.perf_phys_tag_w[p] << " │ "
           << std::setw(5) << std::right << state.perf_phys_pred_w[p] << " │ "
           << std::setw(5) << std::right << state.perf_phys_hyst_w[p] << " │ "
           << std::setw(2) << std::right << state.perf_phys_u_w[p] << " │\n";
    }
    os << "└──────┴──────┴───────┴───────┴────┘\n";

    os << "\n┌─ Hyst Skip (Saturated + Correct) ───────────────────────────────┐\n";
    os << "│ P1 Skip Count:            " << std::setw(36) << std::left << state.perf_hyst_skip_p1 << "│\n";
    os << "│ Bimodal Skip Count:       " << std::setw(36) << std::left << state.perf_hyst_skip_bim << "│\n";
    os << "│ TAGE Skip Count:          " << std::setw(36) << std::left << state.perf_hyst_skip_tage << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ Weak & Wrong Events ───────────────────────────────────────────┐\n";
    os << "│ P1 Weak&Wrong:            " << std::setw(36) << std::left << state.perf_weak_wrong_p1 << "│\n";
    os << "│ Bim Weak&Wrong:           " << std::setw(36) << std::left << state.perf_weak_wrong_bim << "│\n";
    os << "│ TAGE Weak&Wrong:          " << std::setw(36) << std::left << state.perf_weak_wrong_tage << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ P1 Micro-TAGE (P1 vs P2) ──────────────────────────────────────┐\n";
    os << "│ P1 Branch Slots:         " << std::setw(36) << std::left << state.perf_mt_p1_slots << "│\n";
    os << "│ P1 Match/Disagree P2:    " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_p1_match_p2) + " / " + std::to_string(state.perf_mt_p1_disagree_p2)) << "│\n";
    if (state.perf_mt_p1_slots > 0) {
        os << "│ P1 Match Rate vs P2:     " << std::setw(35) << std::left
           << (std::to_string(100.0 * state.perf_mt_p1_match_p2 / state.perf_mt_p1_slots)).substr(0, 6) + "%" << "│\n";
    } else {
        os << "│ P1 Match Rate vs P2:     " << std::setw(36) << std::left << "N/A" << "│\n";
    }
    os << "│ P1 Source Micro/Gshare:  " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_p1_from_micro) + " / " + std::to_string(state.perf_mt_p1_from_gshare)) << "│\n";
    os << "│ Micro Match/Disagree:    " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_p1_micro_match_p2) + " / " + std::to_string(state.perf_mt_p1_micro_disagree_p2)) << "│\n";
    os << "│ Gshare Match/Disagree:   " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_p1_gshare_match_p2) + " / " + std::to_string(state.perf_mt_p1_gshare_disagree_p2)) << "│\n";
    os << "│ Provider Hit/ Wrong:     " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_provider_hit) + " / " + std::to_string(state.perf_mt_provider_wrong)) << "│\n";
    os << "│ UClr Gate Inc/Dec:       " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_uclr_ctr_gate) + " / " +
           std::to_string(state.perf_mt_uclr_ctr_inc) + " / " +
           std::to_string(state.perf_mt_uclr_ctr_dec)) << "│\n";
    os << "│ UClr Sat / UReset:       " << std::setw(36) << std::left
       << (std::to_string(state.perf_mt_uclr_sat) + " / " + std::to_string(state.perf_mt_u_reset)) << "│\n";
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    os << "\n┌─ P1 Predictor Accuracy (vs P2) ────────────────────────────────┐\n";
    os << "│ Predictor │ Count │ Match │ Disagree │ Accuracy │\n";
    os << "├───────────┼───────┼───────┼──────────┼──────────┤\n";
    os << "│ Gshare    │ " << std::setw(5) << std::right << state.perf_mt_p1_from_gshare << " │ "
       << std::setw(5) << std::right << state.perf_mt_p1_gshare_match_p2 << " │ "
       << std::setw(8) << std::right << state.perf_mt_p1_gshare_disagree_p2 << " │ ";
    if (state.perf_mt_p1_from_gshare > 0) {
        os << std::fixed << std::setprecision(2)
           << std::setw(8) << std::right
           << (100.0 * state.perf_mt_p1_gshare_match_p2 / state.perf_mt_p1_from_gshare) << "% │\n";
    } else {
        os << "    N/A │\n";
    }
    os << "│ MicroTot  │ " << std::setw(5) << std::right << state.perf_mt_p1_from_micro << " │ "
       << std::setw(5) << std::right << state.perf_mt_p1_micro_match_p2 << " │ "
       << std::setw(8) << std::right << state.perf_mt_p1_micro_disagree_p2 << " │ ";
    if (state.perf_mt_p1_from_micro > 0) {
        os << std::fixed << std::setprecision(2)
           << std::setw(8) << std::right
           << (100.0 * state.perf_mt_p1_micro_match_p2 / state.perf_mt_p1_from_micro) << "% │\n";
    } else {
        os << "    N/A │\n";
    }
    os << "│ P1 Total  │ " << std::setw(5) << std::right << state.perf_mt_p1_slots << " │ "
       << std::setw(5) << std::right << state.perf_mt_p1_match_p2 << " │ "
       << std::setw(8) << std::right << state.perf_mt_p1_disagree_p2 << " │ ";
    if (state.perf_mt_p1_slots > 0) {
        os << std::fixed << std::setprecision(2)
           << std::setw(8) << std::right
           << (100.0 * state.perf_mt_p1_match_p2 / state.perf_mt_p1_slots) << "% │\n";
    } else {
        os << "    N/A │\n";
    }
    os << "└───────────┴───────┴───────┴──────────┴──────────┘\n";

    os << "\n┌─ P1 Micro-TAGE Table Stats (TAGE-style) ───────────────────────────────────────────────────────────────┐\n";
    os << "│ Tbl │ Reads │ Hits │ Hit%  │ Use │ Correct │ UseAcc% │ Alloc │ C0(weak) │ C1(mwk) │ C2(mst) │ C3(str) │\n";
    os << "├─────┼───────┼──────┼───────┼─────┼─────────┼─────────┼───────┼──────────┼─────────┼─────────┼─────────┤\n";
    for (std::size_t t = 0; t < 4; t++) {
        std::uint64_t reads = state.perf_mt_table_reads[t];
        std::uint64_t hits = state.perf_mt_table_hits[t];
        std::uint64_t use = state.perf_mt_provider_use_table[t];
        std::uint64_t ok = state.perf_mt_provider_correct_table[t];
        std::uint64_t alloc = state.perf_mt_table_alloc[t];
        std::uint64_t c0 = state.perf_mt_conf[t][0], c1 = state.perf_mt_conf[t][1], c2 = state.perf_mt_conf[t][2], c3 = state.perf_mt_conf[t][3];
        std::uint64_t ct = c0 + c1 + c2 + c3;

        os << "│ " << std::setw(3) << std::left << t << " │ "
           << std::setw(5) << std::right << reads << " │ "
           << std::setw(4) << std::right << hits << " │ ";
        if (reads > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(5) << std::right << (100.0 * hits / reads) << "% │ ";
        } else {
            os << "  N/A │ ";
        }
        os << std::setw(3) << std::right << use << " │ "
           << std::setw(7) << std::right << ok << " │ ";
        if (use > 0) {
            os << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0 * ok / use) << "% │ ";
        } else {
            os << "    N/A │ ";
        }
        os << std::setw(5) << std::right << alloc << " │ ";
        if (ct > 0) {
            os << std::setw(5) << c0 << "(" << std::fixed << std::setprecision(0) << std::setw(2) << (100.0 * c0 / ct) << "%) │ ";
            os << std::setw(5) << c1 << "(" << std::setw(2) << (100.0 * c1 / ct) << "%) │ ";
            os << std::setw(5) << c2 << "(" << std::setw(2) << (100.0 * c2 / ct) << "%) │ ";
            os << std::setw(5) << c3 << "(" << std::setw(2) << (100.0 * c3 / ct) << "%) │\n";
        } else {
            os << "     N/A │      N/A │      N/A │      N/A │\n";
        }
    }
    os << "└─────┴───────┴──────┴───────┴─────┴─────────┴─────────┴───────┴──────────┴─────────┴─────────┴─────────┘\n";

    os << "\n┌─ P1 Micro-TAGE Per-Table Events ────────────────────────────────┐\n";
    os << "│ Tbl │ Hit │ Use │ AllocPk │ AllocReq │ PredReq │ HystReq │ UReq │ UClrReq │ TagW │ PredW │ HystW │ UW │\n";
    os << "├─────┼─────┼─────┼─────────┼──────────┼─────────┼─────────┼──────┼─────────┼──────┼───────┼───────┼────┤\n";
    for (std::size_t t = 0; t < 4; t++) {
        os << "│ " << std::setw(3) << std::left << t << " │ "
           << std::setw(3) << std::right << state.perf_mt_provider_hit_table[t] << " │ "
           << std::setw(3) << std::right << state.perf_mt_provider_use_table[t] << " │ "
           << std::setw(7) << std::right << state.perf_mt_alloc_pick_table[t] << " │ "
           << std::setw(8) << std::right << state.perf_mt_alloc_req_table[t] << " │ "
           << std::setw(7) << std::right << state.perf_mt_pred_req_table[t] << " │ "
           << std::setw(7) << std::right << state.perf_mt_hyst_req_table[t] << " │ "
           << std::setw(4) << std::right << state.perf_mt_u_req_table[t] << " │ "
           << std::setw(7) << std::right << state.perf_mt_uclear_req_table[t] << " │ "
           << std::setw(4) << std::right << state.perf_mt_tag_w_table[t] << " │ "
           << std::setw(5) << std::right << state.perf_mt_pred_w_table[t] << " │ "
           << std::setw(5) << std::right << state.perf_mt_hyst_w_table[t] << " │ "
           << std::setw(2) << std::right << state.perf_mt_u_w_table[t] << " │\n";
    }
    os << "└─────┴─────┴─────┴─────────┴──────────┴─────────┴─────────┴──────┴─────────┴──────┴───────┴───────┴────┘\n";

    os << "\n┌─ Verification ──────────────────────────────────────────────────┐\n";
    os << "│ Source total matches predictions: ";
    if (total_all == state.perf_predictions) {
        os << "✓ PASS                      │\n";
    } else {
        os << "✗ FAIL (" << total_all << " vs " << state.perf_predictions << ")    │\n";
    }
    os << "│ Mispred blame sum check:      ";
    if (total_mispred == (state.perf_mispred_blame_tage + state.perf_mispred_blame_sc + state.perf_mispred_blame_p1)) {
        os << "✓ PASS                      │\n";
    } else {
        os << "✗ FAIL (" << total_mispred << " vs "
           << (state.perf_mispred_blame_tage + state.perf_mispred_blame_sc + state.perf_mispred_blame_p1)
           << ")    │\n";
    }
    os << "└─────────────────────────────────────────────────────────────────┘\n";

    return PerfReportTotals{.total_all = total_all, .total_mispred = total_mispred};
}

} // namespace my_bp_v1_perf
