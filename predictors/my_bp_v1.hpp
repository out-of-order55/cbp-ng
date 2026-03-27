// this is a basic TAGE, not necessarily well optimized

#define USE_ALT
#define RESET_UBITS

// #define GATE
#define MY_SC
#if !defined(SC_DISABLE_FGEHL)
#define SC_FGEHL
#endif
// #define SC_DISABLE_BIAS
// #define SC_DISABLE_GEHL
// #define SC_DISABLE_FGEHL

#if !defined(SC_DISABLE_BIAS)
#define SC_USE_BIAS
#endif

#if !defined(SC_DISABLE_GEHL)
#define SC_USE_GEHL
#endif

#ifndef MY_BP_V1_GHYST_USE_WB_RWRAM
#define MY_BP_V1_GHYST_USE_WB_RWRAM 1
#endif

#ifndef MY_BP_V1_GHYST_HIT_UPDATE
#define MY_BP_V1_GHYST_HIT_UPDATE 0
#endif

#ifndef SC_GLOBAL_THRE_INIT
#define SC_GLOBAL_THRE_INIT 23
#endif

#ifndef SC_GLOBAL_THRE_PIPE_STAGES
#define SC_GLOBAL_THRE_PIPE_STAGES 2
#endif

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#ifdef DEBUG_ENERGY
#include "tutorial/energy_monitor.hpp"
#endif
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <sqlite3.h>

using namespace hcm;
#ifdef DEBUG_ENERGY
    struct energy_monitor monitor;
#endif
template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=100, u64 LOGP1=14, u64 GHIST1=6,u64 LOGBANKS = 1,u64 LOGBIAS = 11>
struct my_bp_v1 : predictor {
    // provides 2^(LOGLB-2) predictions per cycle
    // P2 is a TAGE, P1 is a gshare
    static_assert(LOGLB>2);
    static_assert(NUMG>0);

#ifdef DEBUG_ENERGY
    void energy_cycle_reset()
    {
        monitor.last_file = "reset";
        monitor.last_line = 0;
        monitor.last_energy_fJ = panel.energy_fJ();
    }

    void energy_mark(const char *label)
    {
        monitor.record(label, 0);
    }
#endif

    //TODO:need review
    static constexpr u64 PERCWIDTH = 6;
    static constexpr u64 NUMGEHL = 2;
    static constexpr u64 LOGGEHL = 10;
    static_assert(SC_GLOBAL_THRE_PIPE_STAGES >= 1);
#ifdef SC_FGEHL
    static constexpr u64 LOGFGEHL = LOGGEHL;
    static constexpr u64 FHIST_BITS = 48;
#endif
    static constexpr u64 TOTAL_THREBITS = 10;
    static constexpr u64 GLOBAL_THREBITS = 9;
    static constexpr u64 PRE_PC_THREBITS = 9;
    // static constexpr u64 LOGTHREBITS = 5;

    static constexpr u64 MINHIST = 2;
    static constexpr u64 METABITS = 4;
    static constexpr u64 UCTRBITS = 8;
    static constexpr u64 PATHBITS = 6;

    static constexpr u64 LOGGATEBITS = 5;
    static constexpr u64 TAGGATEBITS = 8;
    static constexpr u64 CTRGATEBITS = 3;
#if MY_BP_V1_GHYST_USE_WB_RWRAM
  #if WB_RWRAM_DROP_OLDEST_ON_CONFLICT
    #if MY_BP_V1_GHYST_HIT_UPDATE
    static constexpr const char *GHYST_IMPL_NAME = "wb_rwram_d1_drop_oldest_hit_update";
    #else
    static constexpr const char *GHYST_IMPL_NAME = "wb_rwram_d1_drop_oldest_hit_overwrite";
    #endif
  #else
    #if MY_BP_V1_GHYST_HIT_UPDATE
    static constexpr const char *GHYST_IMPL_NAME = "wb_rwram_d1_keep_oldest_hit_update";
    #else
    static constexpr const char *GHYST_IMPL_NAME = "wb_rwram_d1_keep_oldest_hit_overwrite";
    #endif
  #endif
#else
    static constexpr const char *GHYST_IMPL_NAME = "rwram";
#endif

#ifdef USE_ALT
    static constexpr u64 METAPIPE = 2;
#endif
    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGB > LOGLINEINST);
    static_assert(LOGGEHL > LOGLINEINST);
#ifdef SC_FGEHL
    static_assert(LOGFGEHL > LOGLINEINST);
#endif
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 bindex_bits = LOGB-LOGLINEINST;
    static_assert(TAGW > LOGLINEINST); // the unhashed line offset is part of the tag

    static constexpr u64 HTAGBITS = TAGW-LOGLINEINST; // hashed tag bits


    geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;

    reg<1> true_block = 1;

    // for P1
    reg<GHIST1> global_history1;
    reg<index1_bits> index1;
    arr<reg<1>,LINEINST> readp1; // prediction bits read from P1 table for each offset
    reg<LINEINST> p1; // P1 predictions





    arr<reg<1>,LINEINST> pred2; // alternate P2 prediction for each offset
    arr<reg<1>,LINEINST> use_sc; // alternate P2 prediction for each offset
    reg<LINEINST> p2; // final P2 predictions
    reg<LINEINST> tage_p2; // final P2 predictions
    reg<LINEINST> sc_p2; // final P2 predictions
    arr<reg<1>,LINEINST> sc_pred;
    arr<reg<NUMG+1>,LINEINST> match; // all matches for each offset
    arr<reg<NUMG+1>,LINEINST> match1; // longest match for each offset
    arr<reg<NUMG+1>,LINEINST> match2; // second longest match for each offset
    arr<reg<1>,LINEINST> prov_weak;
    arr<reg<1>,LINEINST> prov_mid; 
    arr<reg<1>,LINEINST> prov_sat;  
    // for P2
    reg<bindex_bits> bindex; // bimodal table index
    arr<reg<LOGG>,NUMG> gindex; // global tables indexes
    arr<reg<HTAGBITS>,NUMG> htag; // computed hashed tags
    arr<reg<1>,LINEINST> readb; // read bimodal prediction bit for each offset
    arr<reg<TAGW>,NUMG> readt; // read tags
    arr<reg<1>,NUMG> readc; // read predictions
    arr<reg<2>,NUMG> readh; // read hysteresis
    arr<reg<1>,NUMG> readu; // read u bits
    reg<NUMG> notumask; // read u bits, inverted 

    arr<reg<1>,LINEINST> pred1; // primary P2 prediction for each offset


#ifdef USE_ALT
    arr<reg<METABITS,i64>,METAPIPE> meta; // select between pred1 and pred2
    arr<reg<1>,LINEINST> newly_alloc;
#endif

#ifdef RESET_UBITS
    reg<UCTRBITS> uctr; // u bits counter (reset u bits when counter saturates)
#endif

    // simulation artifacts, hardware cost may not be real
    u64 num_branch = 0;
    u64 block_size = 0;
    arr<reg<LOGLINEINST>,LINEINST> branch_offset;
    arr<reg<1>,LINEINST> branch_dir;
    arr<reg<64>,LINEINST> branch_pc;
    reg<LINEINST> block_entry; // one-hot vector

#ifdef PERF_COUNTERS
    // Store prediction source per offset (set in predict2, read in update_cycle)
    // 0=bimodal, 1=provider, 2=alt
    arr<reg<2>,LINEINST> pred_source_stored;
    arr<reg<NUMG+1>,LINEINST> pred_match1_stored; // stored match1 for hit info
    arr<reg<NUMG+1>,LINEINST> pred_match2_stored; // stored match2 for alt info

    // Overall statistics
    u64 perf_predictions = 0;
    u64 perf_correct = 0;

    // Per-table source tracking
    u64 perf_provider_used[NUMG] = {};
    u64 perf_provider_correct[NUMG] = {};
    u64 perf_provider_wrong[NUMG] = {};
    u64 perf_alt_used[NUMG] = {};
    u64 perf_alt_correct[NUMG] = {};
    u64 perf_alt_wrong[NUMG] = {};

    // Bimodal tracking
    u64 perf_bimodal_used = 0;
    u64 perf_bimodal_correct = 0;
    u64 perf_bimodal_wrong = 0;

    // Final misprediction blame (exclusive partition)
    u64 perf_mispred_blame_tage = 0;
    u64 perf_mispred_blame_sc = 0;
    u64 perf_mispred_blame_p1 = 0;

    // Table statistics
    u64 perf_table_reads[NUMG] = {};
    u64 perf_table_hits[NUMG] = {};
    u64 perf_table_alloc[NUMG] = {};

    // Allocation failures
    u64 perf_alloc_failures = 0;
    u64 perf_alloc_fail_highest = 0;  // already at highest table
    u64 perf_alloc_fail_noubit = 0;   // no ubit=0 victim found

    // Extra cycle condition counters
    u64 perf_extra_cycle_total = 0;
    u64 perf_extra_cycle_badpred = 0;   // weak & wrong TAGE prediction
    u64 perf_extra_cycle_mispredict = 0;
    u64 perf_extra_cycle_p1_update = 0; // P1 disagree
#ifdef GATE
    u64 perf_gate_count = 0;         // times gating==1 (TAGE skipped, fall back to P1)
    u64 perf_gate_mispred = 0;       // gating==1 and block mispredicted
    u64 perf_gate_update_inc = 0;    // gate_inc triggered (no-disagree, correct, hit)
    u64 perf_gate_update_dec = 0;    // gate_dec triggered (mispred, hit)
#endif
#ifdef MY_SC
    u64 perf_extra_cycle_sc_update = 0; // SC needs update
    // SC prediction source counters
    u64 perf_sc_override = 0;        // SC overrode TAGE (p2 != tage_p2)
    u64 perf_sc_override_correct = 0; // SC override was correct
    u64 perf_sc_use = 0;             // SC selected (use_sc=1)
    u64 perf_sc_use_correct = 0;     // selected SC prediction correct
    u64 perf_sc_use_taken = 0;       // selected SC predicted taken
    u64 perf_sc_use_nottaken = 0;    // selected SC predicted not-taken
    u64 perf_sc_use_same_as_tage = 0; // use_sc but SC dir == TAGE dir
    u64 perf_sc_use_flip_tage = 0;   // use_sc and SC dir != TAGE dir
    u64 perf_sc_use_weak = 0;        // use_sc on weak provider
    u64 perf_sc_use_mid = 0;         // use_sc on mid provider
    u64 perf_sc_use_sat = 0;         // use_sc on high-confidence provider
    u64 perf_sc_stage_prov_hit = 0;  // branch has provider hit
    u64 perf_sc_stage_do_update = 0; // do_update stage passed
    u64 perf_sc_stage_candidate = 0; // provider-hit & do_update
    u64 perf_sc_stage_guard_pass = 0; // candidate & threshold guard
    u64 perf_sc_skip_no_provider = 0; // blocked by no provider
    u64 perf_sc_skip_no_do_update = 0; // blocked by do_update=0
    u64 perf_sc_skip_guard = 0;      // blocked by threshold guard
    u64 perf_global_thre_update = 0; // delayed global threshold update events
    u64 perf_global_thre_inc = 0;    // global threshold incremented
    u64 perf_global_thre_dec = 0;    // global threshold decremented
    // threshold update counters
    u64 perf_thre_update = 0;        // thre1 updated
    u64 perf_thre_update_inc = 0;    // thre1 incremented (sc wrong)
    u64 perf_thre_update_dec = 0;    // thre1 decremented (sc correct)
    // SC misprediction diagnostics
    u64 perf_mispred_sc_not_used = 0;       // use_sc=0 and final wrong
    u64 perf_mispred_sc_keep = 0;           // use_sc=1, sc==tage, final wrong
    u64 perf_mispred_sc_flip = 0;           // use_sc=1, sc!=tage, final wrong
    u64 perf_mispred_sc_flip_harmful = 0;   // flip, tage correct, sc wrong
    u64 perf_mispred_sc_flip_both_wrong = 0; // flip, tage wrong, sc wrong

#endif

    // Confidence distribution per table: [table][0..3] = ctr value buckets
    // 0=weakest(0), 1=weak(1), 2=strong(2), 3=strongest(3)
    u64 perf_conf[NUMG][4] = {};

    // SQLite streaming trace
    sqlite3 *trace_db = nullptr;
    sqlite3_stmt *trace_stmt = nullptr;
    u64 trace_seq = 0;

    void open_trace_db() {
        sqlite3_open("trace_v1.db", &trace_db);
        sqlite3_exec(trace_db, "PRAGMA journal_mode=WAL;", nullptr, nullptr, nullptr);
        sqlite3_exec(trace_db, "PRAGMA synchronous=NORMAL;", nullptr, nullptr, nullptr);
        sqlite3_exec(trace_db,
            "CREATE TABLE IF NOT EXISTS trace("
            "seq INTEGER, cycle INTEGER, pc INTEGER, offset INTEGER,"
            "actual_dir INTEGER, predicted_dir INTEGER, mispredict INTEGER,"
            "pred_source TEXT, pred_table INTEGER, bim_index INTEGER,"
            "hit INTEGER, hit_table INTEGER, hit_gtag INTEGER, hit_gindex INTEGER,"
            "alloc INTEGER, alloc_table INTEGER, alloc_gindex INTEGER, alloc_tag INTEGER,"
            "conf INTEGER, sc_override INTEGER, sc_dir INTEGER, sc_sum INTEGER, threshold INTEGER"
            ");",
            nullptr, nullptr, nullptr);
        sqlite3_exec(trace_db, "BEGIN;", nullptr, nullptr, nullptr);
        sqlite3_prepare_v2(trace_db,
            "INSERT INTO trace VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
            -1, &trace_stmt, nullptr);
    }

    void flush_trace_db() {
        sqlite3_exec(trace_db, "COMMIT;", nullptr, nullptr, nullptr);
        sqlite3_exec(trace_db, "BEGIN;", nullptr, nullptr, nullptr);
    }

    void close_trace_db() {
        sqlite3_exec(trace_db, "COMMIT;", nullptr, nullptr, nullptr);
        sqlite3_finalize(trace_stmt);
        sqlite3_close(trace_db);
        trace_db = nullptr; trace_stmt = nullptr;
    }

    // Misprediction trace
    struct MispredRecord {
        u64 count = 0;
        u64 actual_dir = 0;
        u64 hit = 0;
        u64 hit_table = 0;
        u64 hit_gtag = 0;
        u64 hit_gindex = 0;
    };
    std::unordered_map<u64, MispredRecord> mispred_db;

    // SQLite streaming trace
    void insert_trace(u64 cycle, u64 pc, u64 offset,
                      u64 actual_dir, u64 predicted_dir, u64 mispredict,
                      u64 pred_source, u64 pred_table, u64 bim_index,
                      u64 hit, u64 hit_table, u64 hit_gtag, u64 hit_gindex,
                      u64 alloc, u64 alloc_table, u64 alloc_gindex, u64 alloc_tag,
                      u64 conf, u64 sc_override, i64 sc_dir, i64 sc_sum_val, u64 threshold_val) {
        if (!trace_stmt) return;
        const char *src_str = (pred_source == 0) ? "bimodal" :
                              (pred_source == 1) ? "tage_prov" : "tage_alt";
        sqlite3_reset(trace_stmt);
        sqlite3_bind_int64(trace_stmt,  1, static_cast<i64>(trace_seq++));
        sqlite3_bind_int64(trace_stmt,  2, static_cast<i64>(cycle));
        sqlite3_bind_int64(trace_stmt,  3, static_cast<i64>(pc));
        sqlite3_bind_int64(trace_stmt,  4, static_cast<i64>(offset));
        sqlite3_bind_int64(trace_stmt,  5, static_cast<i64>(actual_dir));
        sqlite3_bind_int64(trace_stmt,  6, static_cast<i64>(predicted_dir));
        sqlite3_bind_int64(trace_stmt,  7, static_cast<i64>(mispredict));
        sqlite3_bind_text (trace_stmt,  8, src_str, -1, SQLITE_STATIC);
        sqlite3_bind_int64(trace_stmt,  9, pred_source != 0 ? static_cast<i64>(pred_table) : -1);
        sqlite3_bind_int64(trace_stmt, 10, static_cast<i64>(bim_index));
        sqlite3_bind_int64(trace_stmt, 11, static_cast<i64>(hit));
        sqlite3_bind_int64(trace_stmt, 12, hit ? static_cast<i64>(hit_table)  : -1);
        sqlite3_bind_int64(trace_stmt, 13, hit ? static_cast<i64>(hit_gtag)   : -1);
        sqlite3_bind_int64(trace_stmt, 14, hit ? static_cast<i64>(hit_gindex) : -1);
        sqlite3_bind_int64(trace_stmt, 15, static_cast<i64>(alloc));
        sqlite3_bind_int64(trace_stmt, 16, alloc ? static_cast<i64>(alloc_table)  : -1);
        sqlite3_bind_int64(trace_stmt, 17, alloc ? static_cast<i64>(alloc_gindex) : -1);
        sqlite3_bind_int64(trace_stmt, 18, alloc ? static_cast<i64>(alloc_tag)    : -1);
        sqlite3_bind_int64(trace_stmt, 19, static_cast<i64>(conf));
        sqlite3_bind_int64(trace_stmt, 20, static_cast<i64>(sc_override));
        sqlite3_bind_int64(trace_stmt, 21, sc_dir);
        sqlite3_bind_int64(trace_stmt, 22, sc_sum_val);
        sqlite3_bind_int64(trace_stmt, 23, static_cast<i64>(threshold_val));
        sqlite3_step(trace_stmt);
        if (trace_seq % 100000 == 0) {
            sqlite3_exec(trace_db, "COMMIT;", nullptr, nullptr, nullptr);
            sqlite3_exec(trace_db, "BEGIN;",  nullptr, nullptr, nullptr);
        }
    }

    void print_perf_counters() {
        std::cerr << "\n╔════════════════════════════════════════════════════════════════╗\n";
        std::cerr << "║           TAGE PREDICTOR PERFORMANCE COUNTERS                   ║\n";
        std::cerr << "╚════════════════════════════════════════════════════════════════╝\n";

        // Overall statistics
        std::cerr << "\n┌─ Overall Statistics ────────────────────────────────────────────┐\n";
        std::cerr << "│ Total Predictions: " << std::setw(43) << std::left << perf_predictions << "│\n";
        std::cerr << "│ Correct:           " << std::setw(43) << std::left << perf_correct << "│\n";
        if (perf_predictions > 0) {
            double accuracy = (100.0 * perf_correct) / perf_predictions;
            std::cerr << "│ Accuracy:          " << std::fixed << std::setprecision(2)
                      << std::setw(40) << std::left << accuracy << "% │\n";
        }
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Per-table statistics with confidence distribution
        std::cerr << "\n┌─ TAGE Table Statistics ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
        std::cerr << "│ Tbl │ HistLen │ Reads │ Hits │ Hit%  │ Prov │ PrvAcc% │ Alt │ AltAcc% │ Total │ TotAcc% │ Alloc │ C0(weak) │ C1(mwk) │ C2(mst) │ C3(str) │\n";
        std::cerr << "├─────┼─────────┼───────┼──────┼───────┼──────┼─────────┼─────┼─────────┼───────┼─────────┼───────┼──────────┼─────────┼─────────┼─────────┤\n";

        for (u64 j=0; j<NUMG; j++) {
            u64 reads = perf_table_reads[j];
            u64 hits = perf_table_hits[j];
            u64 prov = perf_provider_used[j];
            u64 prov_ok = perf_provider_correct[j];
            u64 alt = perf_alt_used[j];
            u64 alt_ok = perf_alt_correct[j];
            u64 total = prov + alt;
            u64 total_ok = prov_ok + alt_ok;
            u64 alloc = perf_table_alloc[j];
            u64 c0=perf_conf[j][0], c1=perf_conf[j][1], c2=perf_conf[j][2], c3=perf_conf[j][3];
            u64 ct = c0+c1+c2+c3;

            std::cerr << "│ " << std::setw(3) << std::left << j << " │ ";
            std::cerr << std::setw(7) << std::right << gfolds.HLEN[j] << " │ ";
            std::cerr << std::setw(5) << std::right << reads << " │ ";
            std::cerr << std::setw(4) << std::right << hits << " │ ";
            if (reads > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(5) << std::right << (100.0*hits/reads) << "% │ ";
            else
                std::cerr << "  N/A │ ";
            std::cerr << std::setw(4) << std::right << prov << " │ ";
            if (prov > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0*prov_ok/prov) << "% │ ";
            else
                std::cerr << "    N/A │ ";
            std::cerr << std::setw(3) << std::right << alt << " │ ";
            if (alt > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0*alt_ok/alt) << "% │ ";
            else
                std::cerr << "    N/A │ ";
            std::cerr << std::setw(5) << std::right << total << " │ ";
            if (total > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right << (100.0*total_ok/total) << "% │ ";
            else
                std::cerr << "    N/A │ ";
            std::cerr << std::setw(5) << std::right << alloc << " │ ";
            if (ct > 0) {
                std::cerr << std::setw(5) << c0 << "(" << std::fixed << std::setprecision(0) << std::setw(2) << (100.0*c0/ct) << "%) │ ";
                std::cerr << std::setw(5) << c1 << "(" << std::setw(2) << (100.0*c1/ct) << "%) │ ";
                std::cerr << std::setw(5) << c2 << "(" << std::setw(2) << (100.0*c2/ct) << "%) │ ";
                std::cerr << std::setw(5) << c3 << "(" << std::setw(2) << (100.0*c3/ct) << "%) │\n";
            } else {
                std::cerr << "     N/A │      N/A │      N/A │      N/A │\n";
            }
        }
        std::cerr << "└─────┴─────────┴───────┴──────┴───────┴──────┴─────────┴─────┴─────────┴───────┴─────────┴───────┴──────────┴─────────┴─────────┴─────────┘\n";

        // Prediction source distribution
        std::cerr << "\n┌─ Prediction Source Distribution ────────────────────────────────┐\n";
        std::cerr << "│ Source   │ Count │ Correct │ Accuracy │ % of Total │\n";
        std::cerr << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";

        u64 total_prov = 0, total_prov_ok = 0;
        u64 total_alt = 0, total_alt_ok = 0;
        u64 total_prov_wrong = 0, total_alt_wrong = 0;
        for (u64 j=0; j<NUMG; j++) {
            total_prov += perf_provider_used[j];
            total_prov_ok += perf_provider_correct[j];
            total_prov_wrong += perf_provider_wrong[j];
            total_alt += perf_alt_used[j];
            total_alt_ok += perf_alt_correct[j];
            total_alt_wrong += perf_alt_wrong[j];
        }

        u64 total_all = total_prov + total_alt + perf_bimodal_used;
        u64 total_mispred = perf_predictions - perf_correct;

        // Provider row
        std::cerr << "│ Provider │ " << std::setw(5) << std::right << total_prov << " │ ";
        std::cerr << std::setw(7) << std::right << total_prov_ok << " │ ";
        if (total_prov > 0)
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right
                      << (100.0*total_prov_ok/total_prov) << "% │ ";
        else
            std::cerr << "    N/A │ ";
        if (total_all > 0)
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right
                      << (100.0*total_prov/total_all) << "% │\n";
        else
            std::cerr << "    N/A │\n";

        // Alternate row
        std::cerr << "│ Alternate│ " << std::setw(5) << std::right << total_alt << " │ ";
        std::cerr << std::setw(7) << std::right << total_alt_ok << " │ ";
        if (total_alt > 0)
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right
                      << (100.0*total_alt_ok/total_alt) << "% │ ";
        else
            std::cerr << "    N/A │ ";
        if (total_all > 0)
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right
                      << (100.0*total_alt/total_all) << "% │\n";
        else
            std::cerr << "    N/A │\n";

        // Bimodal row
        std::cerr << "│ Bimodal  │ " << std::setw(5) << std::right << perf_bimodal_used << " │ ";
        std::cerr << std::setw(7) << std::right << perf_bimodal_correct << " │ ";
        if (perf_bimodal_used > 0)
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right
                      << (100.0*perf_bimodal_correct/perf_bimodal_used) << "% │ ";
        else
            std::cerr << "    N/A │ ";
        if (total_all > 0)
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right
                      << (100.0*perf_bimodal_used/total_all) << "% │\n";
        else
            std::cerr << "    N/A │\n";

#ifdef MY_SC
        // SC Override row
        std::cerr << "│ SC Overr │ " << std::setw(5) << std::right << perf_sc_override << " │ ";
        std::cerr << std::setw(7) << std::right << perf_sc_override_correct << " │ ";
        if (perf_sc_override > 0)
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right
                      << (100.0*perf_sc_override_correct/perf_sc_override) << "% │ ";
        else
            std::cerr << "    N/A │ ";
        if (total_all > 0)
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right
                      << (100.0*perf_sc_override/total_all) << "% │\n";
        else
            std::cerr << "    N/A │\n";
        // Threshold update row
        std::cerr << "│ ThreUpd  │ " << std::setw(5) << std::right << perf_thre_update << " │ ";
        std::cerr << std::setw(7) << std::right << perf_thre_update_inc << " inc │ ";
        std::cerr << std::setw(7) << std::right << perf_thre_update_dec << " dec │\n";
#endif
#ifdef GATE
        // Gate counters row
        std::cerr << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";
        std::cerr << "│ Gate cnt │ " << perf_gate_count << " │ ";
        std::cerr << "misp=" << perf_gate_mispred;
        if (perf_gate_count > 0)
            std::cerr << " (" << std::fixed << std::setprecision(1)
                      << (100.0*perf_gate_mispred/perf_gate_count) << "%)";
        std::cerr << " │ inc=" << perf_gate_update_inc
                  << " dec=" << perf_gate_update_dec << " │\n";
#endif

        // Total row
        std::cerr << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";
        std::cerr << "│ Total    │ " << std::setw(5) << std::right << total_all << " │ ";
        u64 total_correct = total_prov_ok + total_alt_ok + perf_bimodal_correct;
        std::cerr << std::setw(7) << std::right << total_correct << " │ ";
        if (total_all > 0) {
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right
                      << (100.0*total_correct/total_all) << "% │ ";
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right
                      << 100.0 << "% │\n";
        } else {
            std::cerr << "    N/A │     N/A │\n";
        }
        std::cerr << "└──────────┴───────┴─────────┴──────────┴────────────┘\n";
#ifdef MY_SC
        std::cerr << "\n┌─ SC Usage Statistics ────────────────────────────────────────────┐\n";
        std::cerr << "│ SC Use Count:           " << std::setw(36) << std::left << perf_sc_use << "│\n";
        std::cerr << "│ SC Use Correct:         " << std::setw(36) << std::left << perf_sc_use_correct << "│\n";
        if (perf_sc_use > 0) {
            std::cerr << "│ SC Use Accuracy:        " << std::setw(35) << std::left
                      << (std::to_string((100.0 * perf_sc_use_correct / perf_sc_use))).substr(0,5) + "%" << "│\n";
        } else {
            std::cerr << "│ SC Use Accuracy:        " << std::setw(36) << std::left << "N/A" << "│\n";
        }
        std::cerr << "│ SC Use Taken/NotTaken:  " << std::setw(36) << std::left
                  << (std::to_string(perf_sc_use_taken) + " / " + std::to_string(perf_sc_use_nottaken)) << "│\n";
        std::cerr << "│ SC Keep/Flip TAGE:      " << std::setw(36) << std::left
                  << (std::to_string(perf_sc_use_same_as_tage) + " / " + std::to_string(perf_sc_use_flip_tage)) << "│\n";
        std::cerr << "│ SC Use Weak/Mid/Sat:    " << std::setw(36) << std::left
                  << (std::to_string(perf_sc_use_weak) + " / " + std::to_string(perf_sc_use_mid) + " / " + std::to_string(perf_sc_use_sat)) << "│\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ SC Update Pipeline ─────────────────────────────────────────────┐\n";
        std::cerr << "│ Provider Hit:           " << std::setw(36) << std::left << perf_sc_stage_prov_hit << "│\n";
        std::cerr << "│ do_update Passed:       " << std::setw(36) << std::left << perf_sc_stage_do_update << "│\n";
        std::cerr << "│ Candidate (hit&update): " << std::setw(36) << std::left << perf_sc_stage_candidate << "│\n";
        std::cerr << "│ Guard Pass:             " << std::setw(36) << std::left << perf_sc_stage_guard_pass << "│\n";
        std::cerr << "│ Skip No Provider:       " << std::setw(36) << std::left << perf_sc_skip_no_provider << "│\n";
        std::cerr << "│ Skip No do_update:      " << std::setw(36) << std::left << perf_sc_skip_no_do_update << "│\n";
        std::cerr << "│ Skip Guard:             " << std::setw(36) << std::left << perf_sc_skip_guard << "│\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ SC Threshold Updates ───────────────────────────────────────────┐\n";
        std::cerr << "│ Local Threshold (tot/inc/dec): " << std::setw(29) << std::left
                  << (std::to_string(perf_thre_update) + " / " + std::to_string(perf_thre_update_inc) + " / " + std::to_string(perf_thre_update_dec)) << "│\n";
        std::cerr << "│ Global Threshold (tot/inc/dec): " << std::setw(28) << std::left
                  << (std::to_string(perf_global_thre_update) + " / " + std::to_string(perf_global_thre_inc) + " / " + std::to_string(perf_global_thre_dec)) << "│\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

#endif

        std::cerr << "\n┌─ Mispred Reason Statistics ──────────────────────────────────────┐\n";
        std::cerr << "│ Total Mispred:           " << std::setw(36) << std::left << total_mispred << "│\n";
        std::cerr << "│ Src Wrong (Prov/Alt/Bim): " << std::setw(35) << std::left
                  << (std::to_string(total_prov_wrong) + " / " + std::to_string(total_alt_wrong) + " / " + std::to_string(perf_bimodal_wrong)) << "│\n";
        std::cerr << "│ Blame TAGE/SC/P1:        " << std::setw(36) << std::left
                  << (std::to_string(perf_mispred_blame_tage) + " / " + std::to_string(perf_mispred_blame_sc) + " / " + std::to_string(perf_mispred_blame_p1)) << "│\n";
#ifdef MY_SC
        std::cerr << "│ SC Wrong use_sc=0:       " << std::setw(36) << std::left << perf_mispred_sc_not_used << "│\n";
        std::cerr << "│ SC Wrong Keep/Flip:      " << std::setw(36) << std::left
                  << (std::to_string(perf_mispred_sc_keep) + " / " + std::to_string(perf_mispred_sc_flip)) << "│\n";
        std::cerr << "│ SC Flip Harm/BothWrong:  " << std::setw(36) << std::left
                  << (std::to_string(perf_mispred_sc_flip_harmful) + " / " + std::to_string(perf_mispred_sc_flip_both_wrong)) << "│\n";
#endif
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Extra cycle statistics
        std::cerr << "\n┌─ Extra Cycle Statistics ──────────────────────────────────────┐\n";
        std::cerr << "│ Total Extra Cycles:     " << std::setw(36) << std::left << perf_extra_cycle_total << "│\n";
        std::cerr << "│   Weak & Wrong (TAGE):  " << std::setw(36) << std::left << perf_extra_cycle_badpred << "│\n";
        std::cerr << "│   Misprediction:        " << std::setw(36) << std::left << perf_extra_cycle_mispredict << "│\n";
        std::cerr << "│   P1 Update:            " << std::setw(36) << std::left << perf_extra_cycle_p1_update << "│\n";
#ifdef MY_SC
        std::cerr << "│   SC Update:            " << std::setw(36) << std::left << perf_extra_cycle_sc_update << "│\n";
#endif
        std::cerr << "└────────────────────────────────────────────────────────���──────┘\n";

        u64 ghyst_reads = 0;
        u64 ghyst_writes = 0;
        u64 ghyst_stale_reads = 0;
        u64 ghyst_drop_writes = 0;
        u64 ghyst_write_hit_events = 0;
#if MY_BP_V1_GHYST_USE_WB_RWRAM
        u64 ghyst_read_pending_hits = 0;
        u64 ghyst_read_pending_depth1 = 0;
#endif
        for (u64 j = 0; j < NUMG; j++) {
            ghyst_reads += ghyst[j].perf_total_reads();
            ghyst_writes += ghyst[j].perf_total_writes();
            ghyst_stale_reads += ghyst[j].perf_total_stale_reads();
            ghyst_drop_writes += ghyst[j].perf_total_drop_writes();
            ghyst_write_hit_events += ghyst[j].perf_total_write_updates();
#if MY_BP_V1_GHYST_USE_WB_RWRAM
            ghyst_read_pending_hits += ghyst[j].perf_total_read_pending_hits();
            ghyst_read_pending_depth1 += ghyst[j].perf_total_read_pending_hit_depth(0);
#endif
        }
        std::cerr << "\n┌─ GHyst RAM Statistics ──────────────────────────────────────────┐\n";
        std::cerr << "│ Impl:                   " << std::setw(36) << std::left << GHYST_IMPL_NAME << "│\n";
        std::cerr << "│ Reads/Writes:           " << std::setw(36) << std::left
                  << (std::to_string(ghyst_reads) + " / " + std::to_string(ghyst_writes)) << "│\n";
        std::cerr << "│ Stale Reads:            " << std::setw(36) << std::left
                  << (std::to_string(ghyst_stale_reads) + (ghyst_reads ? " (" + std::to_string(100.0 * ghyst_stale_reads / ghyst_reads).substr(0,6) + "%)" : "")) << "│\n";
        std::cerr << "│ Dropped Writes:         " << std::setw(36) << std::left
                  << (std::to_string(ghyst_drop_writes) + (ghyst_writes ? " (" + std::to_string(100.0 * ghyst_drop_writes / ghyst_writes).substr(0,6) + "%)" : "")) << "│\n";
#if MY_BP_V1_GHYST_HIT_UPDATE
        std::cerr << "│ Write Updates:          " << std::setw(36) << std::left
                  << (std::to_string(ghyst_write_hit_events) + (ghyst_writes ? " (" + std::to_string(100.0 * ghyst_write_hit_events / ghyst_writes).substr(0,6) + "%)" : "")) << "│\n";
#else
        std::cerr << "│ Write Hit Events:       " << std::setw(36) << std::left
                  << (std::to_string(ghyst_write_hit_events) + (ghyst_writes ? " (" + std::to_string(100.0 * ghyst_write_hit_events / ghyst_writes).substr(0,6) + "%)" : "")) << "│\n";
#endif
#if MY_BP_V1_GHYST_USE_WB_RWRAM
        std::cerr << "│ Read Pending Hits:      " << std::setw(36) << std::left
                  << (std::to_string(ghyst_read_pending_hits) + (ghyst_reads ? " (" + std::to_string(100.0 * ghyst_read_pending_hits / ghyst_reads).substr(0,6) + "%)" : "")) << "│\n";
        std::cerr << "│ Pending Depth1:         " << std::setw(36) << std::left
                  << (std::to_string(ghyst_read_pending_depth1) +
                      (ghyst_read_pending_hits ? " (" + std::to_string(100.0 * ghyst_read_pending_depth1 / ghyst_read_pending_hits).substr(0,6) + "%)" : ""))
                  << "│\n";
#endif
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Allocation failures
        std::cerr << "\n┌─ Allocation Statistics ─────────────────────────────────────────┐\n";
        std::cerr << "│ Allocation Failures: " << std::setw(41) << std::left
                  << perf_alloc_failures << "│\n";
        std::cerr << "│   - Already at highest table: " << std::setw(32) << std::left
                  << perf_alloc_fail_highest << "│\n";
        std::cerr << "│   - No ubit=0 victim found: " << std::setw(34) << std::left
                  << perf_alloc_fail_noubit << "│\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Verification
        std::cerr << "\n┌─ Verification ──────────────────────────────────────────────────┐\n";
        std::cerr << "│ Source total matches predictions: ";
        if (total_all == perf_predictions)
            std::cerr << "✓ PASS                      │\n";
        else
            std::cerr << "✗ FAIL (" << total_all << " vs " << perf_predictions << ")    │\n";
        std::cerr << "│ Mispred blame sum check:      ";
        if (total_mispred == (perf_mispred_blame_tage + perf_mispred_blame_sc + perf_mispred_blame_p1))
            std::cerr << "✓ PASS                      │\n";
        else
            std::cerr << "✗ FAIL (" << total_mispred << " vs "
                      << (perf_mispred_blame_tage + perf_mispred_blame_sc + perf_mispred_blame_p1)
                      << ")    │\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Finalize SQLite trace
        close_trace_db();
        std::cerr << "│ Full exec trace written to: trace_v1.db                         │\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Write per-PC mispred summary + top-20 to stderr
        {
            std::vector<std::pair<u64,MispredRecord>> sorted_db(mispred_db.begin(), mispred_db.end());
            std::sort(sorted_db.begin(), sorted_db.end(),
                [](const auto &a, const auto &b){ return a.second.count > b.second.count; });

            std::cerr << "\n┌─ Top 20 Most Mispredicted PCs (full trace -> trace_v1.db) ���────────────────────────────────┐\n";
            std::cerr << "│ Rank │       PC       │ Mispreds │ Dir │ Hit │ Tbl │    GTTag     │  GIndex  │\n";
            std::cerr << "├──────┼────────────────┼──────────┼─────┼─────┼─────┼──────────────┼──────────┤\n";
            u64 rank = 0;
            for (auto &[pc, rec] : sorted_db) {
                if (rank >= 20) break;
                std::cerr << "│ " << std::setw(4) << std::right << (rank+1) << " │ ";
                std::cerr << "0x" << std::hex << std::setw(12) << std::setfill('0') << pc << std::dec << std::setfill(' ') << " │ ";
                std::cerr << std::setw(8) << rec.count << " │ ";
                std::cerr << std::setw(3) << rec.actual_dir << " │ ";
                std::cerr << std::setw(3) << rec.hit << " │ ";
                if (rec.hit) {
                    std::cerr << std::setw(3) << rec.hit_table << " │ ";
                    std::cerr << "0x" << std::hex << std::setw(10) << std::setfill('0') << rec.hit_gtag << std::dec << std::setfill(' ') << " │ ";
                    std::cerr << std::setw(8) << rec.hit_gindex << " │\n";
                } else {
                    std::cerr << "N/A │          N/A │      N/A │\n";
                }
                rank++;
            }
            std::cerr << "└──────┴────────────────┴──────────┴─────┴─────┴─────┴──────────────┴──────────┘\n";
        }
    }

#ifdef CHEATING_MODE
#ifdef USE_ALT
    void perf_store_prediction_source(val<1> metasign) {
#else
    void perf_store_prediction_source() {
#endif
        for (u64 offset=0; offset<LINEINST; offset++) {
#ifdef USE_ALT
            val<1> use_alt_o = metasign & newly_alloc[offset] & (match2[offset]!=hard<0>{});
            val<1> has_tage = (match1[offset] & hard<(1<<NUMG)-1>{}) != hard<0>{};
            val<2> src = select(use_alt_o.fo1(), val<2>{2},
                         select(has_tage, val<2>{1}, val<2>{0}));
#else
            val<1> has_tage = (match1[offset] & hard<(1<<NUMG)-1>{}) != hard<0>{};
            val<2> src = select(has_tage, val<2>{1}, val<2>{0});
#endif
            pred_source_stored[offset] = src;
            pred_match1_stored[offset] = match1[offset];
            pred_match2_stored[offset] = match2[offset];
        }
    }

    void perf_count_table_reads_hits(arr<val<1>,LINEINST> is_branch) {
        u64 num_branches = 0;
        for (u64 offset=0; offset<LINEINST; offset++) {
            if (static_cast<bool>(is_branch[offset])) num_branches++;
        }

        for (u64 j=0; j<NUMG; j++) {
            perf_table_reads[j] += num_branches;
            for (u64 offset=0; offset<LINEINST; offset++) {
                if (!static_cast<bool>(is_branch[offset])) continue;
                val<NUMG> all_hits = val<NUMG+1>{match[offset]} >> 1;
                val<1> hit_j = (all_hits >> j) & val<1>{1};
                if (static_cast<bool>(hit_j)) {
                    perf_table_hits[j]++;
                }
            }
        }
    }

    void perf_count_table_alloc(arr<val<1>,NUMG> allocate) {
        for (u64 j=0; j<NUMG; j++) {
            if (static_cast<bool>(allocate[j])) {
                perf_table_alloc[j]++;
            }
        }
    }

#ifdef MY_SC
    void perf_count_extra_cycle(val<1> extra_cycle, val<1> some_badpred1, val<1> mispredict, val<1> p1_update, val<1> sc_need_update) {
#else
    void perf_count_extra_cycle(val<1> extra_cycle, val<1> some_badpred1, val<1> mispredict, val<1> p1_update) {
#endif
        perf_extra_cycle_total      += static_cast<u64>(extra_cycle);
        perf_extra_cycle_badpred    += static_cast<u64>(some_badpred1);
        perf_extra_cycle_mispredict += static_cast<u64>(mispredict);
        perf_extra_cycle_p1_update  += static_cast<u64>(p1_update);
#ifdef MY_SC
        perf_extra_cycle_sc_update  += static_cast<u64>(sc_need_update);
#endif
    }

#ifdef GATE
    void perf_count_gate(val<1> gating_local, val<1> mispredict_local, val<1> gate_inc, val<1> gate_dec) {
        perf_gate_count      += static_cast<u64>(gating_local);
        perf_gate_mispred    += static_cast<u64>(gating_local & mispredict_local);
        perf_gate_update_inc += static_cast<u64>(gate_inc);
        perf_gate_update_dec += static_cast<u64>(gate_dec);
    }
#endif

    void perf_count_end_of_cycle(
        arr<val<1>,LINEINST> is_branch,
        arr<val<1>,LINEINST> branch_taken,
        val<1> mispredict,
        arr<val<1>,NUMG> allocate,
        val<NUMG> postmask,
        val<1> noalloc) {
        if (static_cast<bool>(noalloc) && static_cast<bool>(mispredict)) {
            perf_alloc_failures++;
            if (static_cast<bool>(postmask == hard<0>{})) {
                perf_alloc_fail_highest++;
            } else {
                perf_alloc_fail_noubit++;
            }
        }

        for (u64 offset=0; offset<LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

            perf_predictions++;

            u64 src = static_cast<u64>(val<2>{pred_source_stored[offset]});
            val<NUMG> prov_mask = val<NUMG+1>{pred_match1_stored[offset]} & hard<(1<<NUMG)-1>{};
            val<NUMG> alt_mask  = val<NUMG+1>{pred_match2_stored[offset]} & hard<(1<<NUMG)-1>{};

            if (src == 2) {
                for (u64 j=0; j<NUMG; j++) {
                    if (static_cast<u64>(alt_mask >> j) & 1) { perf_alt_used[j]++; break; }
                }
            } else if (src == 1) {
                for (u64 j=0; j<NUMG; j++) {
                    if (static_cast<u64>(prov_mask >> j) & 1) { perf_provider_used[j]++; break; }
                }
            } else {
                perf_bimodal_used++;
            }

            val<1> actual = branch_taken[offset];
            val<1> predicted = (p2 >> offset) & val<1>{1};
            val<1> correct = (predicted == actual);

            if (static_cast<bool>(correct)) {
                perf_correct++;
                if (src == 2) {
                    for (u64 j=0; j<NUMG; j++) {
                        if (static_cast<u64>(alt_mask >> j) & 1) { perf_alt_correct[j]++; break; }
                    }
                } else if (src == 1) {
                    for (u64 j=0; j<NUMG; j++) {
                        if (static_cast<u64>(prov_mask >> j) & 1) { perf_provider_correct[j]++; break; }
                    }
                } else {
                    perf_bimodal_correct++;
                }
            } else {
                if (src == 2) {
                    for (u64 j=0; j<NUMG; j++) {
                        if (static_cast<u64>(alt_mask >> j) & 1) { perf_alt_wrong[j]++; break; }
                    }
                } else if (src == 1) {
                    for (u64 j=0; j<NUMG; j++) {
                        if (static_cast<u64>(prov_mask >> j) & 1) { perf_provider_wrong[j]++; break; }
                    }
                } else {
                    perf_bimodal_wrong++;
                }

#ifdef GATE
                if (static_cast<bool>(gating)) {
                    perf_mispred_blame_p1++;
                } else
#endif
                {
#ifdef MY_SC
                    val<1> tage_p2_bit = (tage_p2 >> offset) & val<1>{1};
                    val<1> use_sc_bit = use_sc[offset];
                    val<1> sc_pred_bit = sc_pred[offset];

                    if (static_cast<bool>(use_sc_bit)) {
                        if (static_cast<bool>(sc_pred_bit == tage_p2_bit)) {
                            perf_mispred_sc_keep++;
                            perf_mispred_blame_tage++;
                        } else {
                            perf_mispred_sc_flip++;
                            if (static_cast<bool>(tage_p2_bit == actual)) {
                                perf_mispred_sc_flip_harmful++;
                                perf_mispred_blame_sc++;
                            } else {
                                perf_mispred_sc_flip_both_wrong++;
                                perf_mispred_blame_tage++;
                            }
                        }
                    } else {
                        perf_mispred_sc_not_used++;
                        perf_mispred_blame_tage++;
                    }
#else
                    perf_mispred_blame_tage++;
#endif
                }
            }

            u64 pc_val = 0;
            for (u64 bi = 0; bi < num_branch; bi++) {
                if (static_cast<u64>(branch_offset[bi]) == offset) {
                    pc_val = static_cast<u64>(branch_pc[bi]);
                    break;
                }
            }

            u64 hit_found = 0, hit_table = NUMG, hit_gtag = 0, hit_gindex = 0;
            u64 pmask = static_cast<u64>(val<NUMG+1>{pred_match1_stored[offset]}) & ((1ULL << NUMG) - 1);
            if (pmask != 0) {
                for (u64 j = 0; j < NUMG; j++) {
                    if ((pmask >> j) & 1) {
                        hit_found  = 1;
                        hit_table  = j;
                        hit_gtag   = static_cast<u64>(readt[j]);
                        hit_gindex = static_cast<u64>(gindex[j]);
                        break;
                    }
                }
            }

            u64 pred_source = src;
            u64 pred_table  = NUMG;
            if (pred_source == 2) {
                u64 amask = static_cast<u64>(alt_mask);
                for (u64 j = 0; j < NUMG; j++) {
                    if ((amask >> j) & 1) { pred_table = j; break; }
                }
            } else if (pred_source == 1) {
                u64 pmask2 = static_cast<u64>(prov_mask);
                for (u64 j = 0; j < NUMG; j++) {
                    if ((pmask2 >> j) & 1) { pred_table = j; break; }
                }
            }

            u64 is_last = (offset == static_cast<u64>(branch_offset[num_branch-1]));
            u64 is_misp = static_cast<u64>(mispredict);
            u64 alloc_found = 0, alloc_table = 0, alloc_gindex_val = 0, alloc_tag_val = 0;
            if (is_misp && is_last) {
                for (u64 j = 0; j < NUMG; j++) {
                    if (static_cast<bool>(allocate[j])) {
                        alloc_found     = 1;
                        alloc_table     = j;
                        alloc_gindex_val = static_cast<u64>(gindex[j]);
                        alloc_tag_val   = static_cast<u64>(concat(branch_offset[num_branch-1], htag[j]));
                        break;
                    }
                }
            }

            {
                u64 conf_val = 0;
                if (pred_source == 1 || pred_source == 2) {
                    u64 tbl = (pred_source == 1) ? pred_table :
                              static_cast<u64>(static_cast<u64>(alt_mask) ? [&]{ for(u64 j=0;j<NUMG;j++) if((static_cast<u64>(alt_mask)>>j)&1) return j; return (u64)0; }() : 0);
                    if (tbl < NUMG) conf_val = static_cast<u64>(readh[tbl]);
                }

                u64 sc_override_val = 0;
                i64 sc_dir_val = -1;
                i64 sc_sum_val = 0;
                u64 threshold_val = 0;
#ifdef MY_SC
                {
                    u64 p2_bit     = static_cast<u64>((p2     >> offset) & val<1>{1});
                    u64 tage_p2_bit= static_cast<u64>((tage_p2>> offset) & val<1>{1});
                    sc_override_val = (p2_bit != tage_p2_bit) ? 1 : 0;
                    if (sc_override_val) sc_dir_val = static_cast<i64>(p2_bit);
                    sc_sum_val   = static_cast<i64>(sc_sum[offset]);
                    threshold_val = static_cast<u64>(threshold[offset]);
                }
#endif
                insert_trace(
                    static_cast<u64>(panel.cycle), pc_val, offset,
                    static_cast<u64>(actual), static_cast<u64>(predicted),
                    is_misp & is_last,
                    pred_source, pred_table, static_cast<u64>(bindex),
                    hit_found, hit_table, hit_gtag, hit_gindex,
                    is_misp & is_last & alloc_found,
                    alloc_table, alloc_gindex_val, alloc_tag_val,
                    conf_val, sc_override_val, sc_dir_val, sc_sum_val, threshold_val);
            }

            if (pred_source == 1 || pred_source == 2) {
                u64 tbl = (pred_source == 1) ? pred_table :
                          static_cast<u64>(static_cast<u64>(alt_mask) ? [&]{ for(u64 j=0;j<NUMG;j++) if((static_cast<u64>(alt_mask)>>j)&1) return j; return (u64)0; }() : 0);
                if (tbl < NUMG) {
                    u64 hval = static_cast<u64>(readh[tbl]);
                    if (hval < 4) perf_conf[tbl][hval]++;
                }
            }

            if (!static_cast<bool>(correct)) {
                auto &rec = mispred_db[pc_val];
                rec.count++;
                rec.actual_dir  = static_cast<u64>(actual);
                rec.hit         = hit_found;
                rec.hit_table   = hit_table;
                rec.hit_gtag    = hit_gtag;
                rec.hit_gindex  = hit_gindex;
            }
        }

#ifdef MY_SC
        for (u64 offset=0; offset<LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;
            val<1> p2_bit = (p2 >> offset) & val<1>{1};
            val<1> tage_p2_bit = (tage_p2 >> offset) & val<1>{1};
            val<1> actual = branch_taken[offset];
            val<1> sc_override = (p2_bit != tage_p2_bit);
            val<1> use_sc_bit = use_sc[offset];
            val<1> sc_pred_bit = sc_pred[offset];
            val<1> prov_hit = prov_hit_arr[offset];
            val<1> do_update = do_update_arr[offset];
            val<1> thre_guard = thre_guard_arr[offset];

            if (static_cast<bool>(use_sc_bit)) {
                perf_sc_use++;
                perf_sc_use_correct += static_cast<u64>(sc_pred_bit == actual);
                if (static_cast<bool>(sc_pred_bit))
                    perf_sc_use_taken++;
                else
                    perf_sc_use_nottaken++;
                if (static_cast<bool>(sc_pred_bit == tage_p2_bit))
                    perf_sc_use_same_as_tage++;
                else
                    perf_sc_use_flip_tage++;
                if (static_cast<bool>(prov_weak[offset]))
                    perf_sc_use_weak++;
                if (static_cast<bool>(prov_mid[offset]))
                    perf_sc_use_mid++;
                if (static_cast<bool>(prov_sat[offset]))
                    perf_sc_use_sat++;
            }

            if (static_cast<bool>(prov_hit))
                perf_sc_stage_prov_hit++;
            if (static_cast<bool>(do_update))
                perf_sc_stage_do_update++;
            if (static_cast<bool>(prov_hit & do_update))
                perf_sc_stage_candidate++;
            if (static_cast<bool>(prov_hit & do_update & thre_guard))
                perf_sc_stage_guard_pass++;

            if (!static_cast<bool>(prov_hit)) {
                perf_sc_skip_no_provider++;
            } else if (!static_cast<bool>(do_update)) {
                perf_sc_skip_no_do_update++;
            } else if (!static_cast<bool>(thre_guard)) {
                perf_sc_skip_guard++;
            }

            if (static_cast<bool>(sc_override)) {
                perf_sc_override++;
                perf_sc_override_correct += static_cast<u64>(p2_bit == actual);
            }
            if (static_cast<bool>(is_branch[offset] & do_update_arr[offset] & thre_guard_arr[offset] & prov_hit_arr[offset])) {
                perf_thre_update++;
                if (static_cast<bool>(sc_wrong_arr[offset]))
                    perf_thre_update_inc++;
                else
                    perf_thre_update_dec++;
            }
        }
#endif
    }
#endif
#endif

#ifdef GATE
    // ram<val<TAGGATEBITS>,(1<<LOGGATEBITS)> gate_tag {"gate_tag"};
    // ram<val<1>,(1<<LOGGATEBITS)> gate_pred {"gate_pred"};
    // ram<val<CTRGATEBITS-1>,(1<<LOGGATEBITS)> gate_hyst {"gate_hyst"};
    // reg<1> gating;
    // reg<TAGGATEBITS> atag[2];
    // reg<1> apred[2];

    arr<reg<TAGGATEBITS>,(1<<LOGGATEBITS)> gate_tag[3];
    arr<reg<CTRGATEBITS>,(1<<LOGGATEBITS)> gate_ctr[3];
    reg<LOGGATEBITS> gate_idx;
    reg<1> gate_hit;
    reg<TAGGATEBITS> line_tag;
    reg<CTRGATEBITS> read_gate_ctr;
    reg<1> gating;
    reg<TAGGATEBITS> atag[2];
    reg<1> apred[2];
    reg<LOGGATEBITS> aidx[2];//for update
#endif  
    
    // Cluster P1 + bpred + SC tables to reduce dominant SC-path wire distance.
    zone P1_SC_CLUSTER;
    // P1 (gshare)
    ram<val<1>,(1<<index1_bits)> table1_pred[LINEINST] {"P1 pred"}; // P1 prediction bit
    ram<val<1>,(1<<bindex_bits)> bim[LINEINST] {"bpred"}; // bimodal prediction bits
    #ifdef MY_SC
    //bias_pc is concat pc and tage
    //idx = (pc>>(LOGLB+2)) 2bit tage info(prov_weak,prov_taken)
    arr<reg<TOTAL_THREBITS>,LINEINST> threshold;
#ifdef SC_USE_BIAS
    ram<arr<val<PERCWIDTH,i64>,4>,(1<<(LOGBIAS-LOGLINEINST-2))> bias_pc[LINEINST] {"Bias pc"};
    arr<reg<PERCWIDTH,i64>,LINEINST> bias_bank_map[4];
    arr<reg<PERCWIDTH,i64>,4> bias_lmap;
    arr<reg<PERCWIDTH,i64>,LINEINST> bias_map;
    arr<reg<2>,LINEINST> tage_info;
    reg<LOGBIAS-LOGLINEINST-2> bias_high_idx;
#endif

    arr<reg<PRE_PC_THREBITS>,LINEINST> thre1;
    reg<GLOBAL_THREBITS> global_thre = val<GLOBAL_THREBITS>{SC_GLOBAL_THRE_INIT};
    arr<reg<2,i64>,SC_GLOBAL_THRE_PIPE_STAGES> global_thre_pipe;
    arr<reg<TOTAL_THREBITS,i64>,LINEINST> sc_sum;
#ifdef SC_USE_GEHL
    arr<reg<LOGGEHL-LOGLINEINST>,NUMGEHL> gehl_idx;
    ram<val<PERCWIDTH,i64>,(1<<(LOGGEHL-LOGLINEINST))> gehl[NUMGEHL][LINEINST] {"GEHL"};
    arr<reg<PERCWIDTH,i64>,LINEINST> gehl_map[NUMGEHL];
#endif
#ifdef SC_FGEHL
    reg<LOGFGEHL-LOGLINEINST> fgehl_idx;
    ram<val<PERCWIDTH,i64>,(1<<(LOGFGEHL-LOGLINEINST))> fgehl[LINEINST] {"FGEHL"};
    arr<reg<PERCWIDTH,i64>,LINEINST> fgehl_map;
    reg<FHIST_BITS> fhist;
#endif
    // intermediate update signals
    arr<reg<1>,LINEINST> sc_wrong_arr;
    arr<reg<1>,LINEINST> do_update_arr;
    arr<reg<1>,LINEINST> prov_hit_arr;
    arr<reg<1>,LINEINST> thre_guard_arr;
    // rwram<PERCWIDTH,(1<<LOGBIAS)/NUMBANKS,NUMBANKS> bias_pc {"Bias pc"};

    
    #endif
    // Put TAGE tables in a dedicated cluster.
    zone TAGE_CLUSTER;
    // P2 (TAGE)
    ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"tags"}; // tags
    // rwram<TAGW,(1<<LOGG),4> gtag[NUMG] {"tags"}; // tags
    ram<val<1>,(1<<LOGG)> gpred[NUMG] {"gpred"}; // predictions
#if MY_BP_V1_GHYST_USE_WB_RWRAM
    wb_rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"}; // hysteresis
#else
    rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"}; // hysteresis
#endif

    rwram<1,(1<<LOGG),4> ubit[NUMG] {"u"}; // "useful" bits
    zone UPDATE_ONLY;
    ram<val<1>,(1<<index1_bits)> table1_hyst[LINEINST] {"P1 hyst"}; // P1 hysteresis
    ram<val<1>,(1<<bindex_bits)> bhyst[LINEINST] {"bhyst"}; // bimodal hysteresis

// #endif
    my_bp_v1()
    {
#ifdef TAGE_VERBOSE
        std::cerr << "TAGE history lengths: ";
        for (u64 i=0; i<NUMG; i++) std::cerr << gfolds.HLEN[i] << " ";
        std::cerr << std::endl;
        if (LOGG == HTAGBITS) {
            std::cerr << "WARNING: the tag function and index function are not different enough\n";
        }
#endif
#ifdef PERF_COUNTERS
        open_trace_db();
#endif
    }


    ~my_bp_v1() {
#ifdef PERF_COUNTERS
        print_perf_counters();
#endif
#ifdef DEBUG_ENERGY
        monitor.report();
#endif
    }


    void new_block(val<64> inst_pc)
    {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        block_entry = offset.fo1().decode().concat();
        block_entry.fanout(hard<6*LINEINST>{});
        block_size = 1;
    }

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
#ifdef DEBUG_ENERGY
        energy_cycle_reset();
#endif
        inst_pc.fanout(hard<3>{});
        new_block(inst_pc);
#ifdef GATE
        gate_predict2(inst_pc);
#endif

        val<std::max(index1_bits,GHIST1)> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<2>{});
        if constexpr (GHIST1 <= index1_bits) {
            index1 = lineaddr ^ (val<index1_bits>{global_history1}<<(index1_bits-GHIST1));
        } else {
            index1 = global_history1.make_array(val<index1_bits>{}).append(lineaddr).fold_xor();
        }

        index1.fanout(hard<LINEINST>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            readp1[offset] = table1_pred[offset].read(index1);
        }
        readp1.fanout(hard<2>{});
        p1 = readp1.concat();
        p1.fanout(hard<LINEINST>{});
#ifdef DEBUG_ENERGY
        energy_mark("P1Predict");
#endif
        return (block_entry & p1) != hard<0>{};
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        return ((block_entry<<block_size) & p1) != hard<0>{};
    };
#ifdef GATE
    void gate_predict2([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<3>{});
        gate_idx = inst_pc >> LOGLB;
        gate_idx.fanout(hard<2*(1<<LOGGATEBITS)>{});
        line_tag = (inst_pc >> (LOGLB+LOGGATEBITS)) ^ (inst_pc >> (LOGLB+LOGGATEBITS+TAGGATEBITS));
        val<TAGGATEBITS> tag =  (inst_pc >> (LOGLB+LOGGATEBITS)) ^ (inst_pc >> (LOGLB+LOGGATEBITS+TAGGATEBITS));

        // aidx[0] = lineaddr;
        
        val<TAGGATEBITS> read_tag = arr<val<TAGGATEBITS>,(1<<LOGGATEBITS)>{[&](u64 i){
            return select(gate_idx == val<LOGGATEBITS>{i}, val<TAGGATEBITS>{gate_tag[2][i]}, val<TAGGATEBITS>{0});
        }}.fold_or();
        read_gate_ctr = arr<val<CTRGATEBITS>,(1<<LOGGATEBITS)>{[&](u64 i){
            return select(gate_idx == val<LOGGATEBITS>{i}, val<CTRGATEBITS>{gate_ctr[2][i]}, val<CTRGATEBITS>{0});
        }}.fold_or();
        
        // atag[1] = atag[0];
        // apred[1] = apred[0];
        // apred[0] = read_pred;
        // atag[0] = read_tag;
        read_tag.fanout(hard<2>{});
        tag.fanout(hard<2>{});
        gate_hit = (read_tag == tag);
        gating = ((read_tag == tag) & (read_gate_ctr == val<CTRGATEBITS>{7}));
    };
#endif
    void tage_pred(val<64> inst_pc){
        val<std::max(bindex_bits,LOGG)> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<1+NUMG*2>{});
        gfolds.fanout(hard<2>{});
        // val<1> gate = gate_predict2(inst_pc);

        // compute indexes
        bindex = lineaddr;
        bindex.fanout(hard<LINEINST>{});
        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = lineaddr ^ gfolds.template get<0>(i);
        }
      

        // compute hashed tags
        for (u64 i=0; i<NUMG; i++) {
            htag[i] = val<HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
        }
        htag.fanout(hard<2>{});

        // read tables
        for (u64 offset=0; offset<LINEINST; offset++) {
            readb[offset] = bim[offset].read(bindex);
        }
        readb.fanout(hard<2>{});
        gindex.fanout(hard<4>{});
#ifdef GATE

            
            gating.fanout(hard<4>{});
            for (u64 i=0; i<NUMG; i++) {
                execute_if(~gating, [&](){
                    readt[i] = gtag[i].read(gindex[i]);
                    readc[i] = gpred[i].read(gindex[i]);
                    readh[i] = ghyst[i].read(gindex[i]);
                    readu[i] = ubit[i].read(gindex[i]);
                });
            }

#else

        for (u64 i=0; i<NUMG; i++) {
            readt[i] = gtag[i].read(gindex[i]);
            readc[i] = gpred[i].read(gindex[i]);
            readh[i] = ghyst[i].read(gindex[i]);
            readu[i] = ubit[i].read(gindex[i]);
        }

#endif
        readt.fanout(hard<LINEINST+1>{});
        readc.fanout(hard<3>{});
#ifdef MY_SC
        readh.fanout(hard<1+4*LINEINST>{});
#endif
        readu.fanout(hard<2>{});
        notumask = ~readu.concat();
        notumask.fanout(hard<2>{});
#ifdef DEBUG_ENERGY
        energy_mark("TAGERead");
#endif

        // gather prediction bits for each offset
        val<NUMG> gpreds = readc.concat();
        gpreds.fanout(hard<LINEINST>{});
        arr<val<NUMG+1>,LINEINST> preds = [&](u64 offset){return concat(readb[offset],gpreds);};
        preds.fanout(hard<2*LINEINST>{});

        // hashed tags comparisons
        arr<val<1>,NUMG> htagcmp_split = [&](int i){return val<HTAGBITS>{readt[i]} == htag[i];};
        val<NUMG> htagcmp = htagcmp_split.fo1().concat();
        htagcmp.fanout(hard<LINEINST>{});

        // generate match mask for each offset
        static_loop<LINEINST>([&]<u64 offset>(){
            arr<val<1>,NUMG> tagcmp = [&](int i){return val<LOGLINEINST>{readt[i]>>HTAGBITS} == hard<offset>{};};
            match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp); // bimodal is default when no match

        });

        match.fanout(hard<2>{});

        // for each offset, find longest match and select primary prediction
        for (u64 offset=0; offset<LINEINST; offset++) {
            match1[offset] = match[offset].one_hot();

        }
#ifdef MY_SC
        match1.fanout(hard<3+3*NUMG>{});
#else
        match1.fanout(hard<3>{});
#endif
        for (u64 offset=0; offset<LINEINST; offset++) {
#ifdef MY_SC
            arr<val<1>,NUMG> weakctr = [&](int i) {return (readh[i]==hard<0>{}) & (val<NUMG>{match1[offset]}!=val<NUMG>{0});};
            arr<val<1>,NUMG> midctr = [&](int i) {return ((readh[i]==hard<1>{}) | (readh[i]==hard<2>{})) & (val<NUMG>{match1[offset]}!=val<NUMG>{0});};
            arr<val<1>,NUMG> highctr = [&](int i) {return (readh[i]==hard<3>{}) & (val<NUMG>{match1[offset]}!=val<NUMG>{0});};
            prov_weak[offset] = weakctr.fold_or();
            prov_mid[offset] = midctr.fold_or();
            prov_sat[offset] = highctr.fold_or();
#endif
            pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
        }
        pred1.fanout(hard<3>{});

        // for each offset, find second longest match and select secondary prediction
        for (u64 offset=0; offset<LINEINST; offset++) {
            match2[offset] = (match[offset]^match1[offset]).one_hot();
        }
        match2.fanout(hard<2>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
        }
        pred2.fanout(hard<2>{});


#ifdef USE_ALT
        meta.fanout(hard<2>{});
        arr<val<1>,NUMG> weakctr = [&](int i) {return readh[i]==hard<0>{};};
        val<NUMG> coldctr = notumask & weakctr.fo1().concat();
        coldctr.fanout(hard<LINEINST>{});
        val<1> metasign = (meta[METAPIPE-1] >= hard<0>{});
        metasign.fanout(hard<LINEINST>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
        }
        newly_alloc.fanout(hard<2>{});
        arr<val<1>,LINEINST> altsel = [&](u64 offset){
            arr<val<1>,3> inputs = {metasign, newly_alloc[offset], match2[offset]!=hard<0>{}};
            return inputs.fo1().fold_and();
        };
        tage_p2 = arr<val<1>,LINEINST> {[&](u64 offset){
            return select(altsel[offset].fo1(),pred2[offset],pred1[offset]);
        }}.concat();
#else
        tage_p2 = pred1.concat();
#endif
        // tage_p2.fanout(hard<2>{});
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
#ifdef USE_ALT
        perf_store_prediction_source(metasign);
#else
        perf_store_prediction_source();
#endif
#endif
#endif
#ifdef DEBUG_ENERGY
        energy_mark("TAGESelect");
#endif
    }
#ifdef MY_SC
    void sc_predict(val<64> inst_pc){

        // inst_pc.fanout(hard<5>{});
        match1.fanout(hard<2>{});
        prov_weak.fanout(hard<2>{});
        // bias
#ifdef SC_USE_BIAS
        val<LOGBIAS-LOGLINEINST-2> bias_pc_high_index =  inst_pc>>(LOGLB+2);
        bias_pc_high_index.fanout(hard<(LINEINST+1)>{});
        bias_map = arr<val<PERCWIDTH,i64>,LINEINST>{[&](u64 offset){
            arr<val<PERCWIDTH,i64>,4> bank_out = bias_pc[offset].read(bias_pc_high_index);
            for (u64 bank = 0; bank < 4; bank++) {
                bias_bank_map[bank][offset] = bank_out[bank];
            }
            val<1> prov_hit  = val<NUMG>{match1[offset]}!=val<NUMG>{0};
            val<1> prov_pred = pred1[offset] & prov_hit;
            val<1> is_prov_weak =  prov_weak[offset];
            tage_info[offset] = concat(prov_pred,is_prov_weak);
            val<PERCWIDTH,i64> result = bank_out.select(tage_info[offset]);
            return result;
        }};
        tage_info.fanout(hard<2>{});
#ifdef DEBUG_ENERGY
        energy_mark("SCBiasRead");
#endif
        bias_map.fanout(hard<2>{});
        bias_high_idx = bias_pc_high_index;
#endif
#ifdef SC_USE_GEHL
        val<LOGGEHL-LOGLINEINST> gehl_base_index = inst_pc >> LOGLB;
        gfolds.fanout(hard<NUMGEHL>{});
        gehl_base_index.fanout(hard<NUMGEHL>{});
        for (u64 k = 0; k < NUMGEHL; k++) {
            val<LOGGEHL-LOGLINEINST> idx = gehl_base_index ^ val<LOGGEHL-LOGLINEINST>{gfolds.template get<0>(k+1)};
            idx.fanout(hard<LINEINST+1>{});
            gehl_idx[k] = idx;
            gehl_map[k] = arr<val<PERCWIDTH,i64>,LINEINST>{[&](u64 offset){
                return gehl[k][offset].read(idx);
            }};
            gehl_map[k].fanout(hard<2>{});
        }
#endif
#ifdef SC_FGEHL
        val<LOGFGEHL-LOGLINEINST> fgehl_base_index = inst_pc >> LOGLB;
        val<FHIST_BITS> fh_mix = val<FHIST_BITS>{fhist} ^ (val<FHIST_BITS>{fhist} >> hard<13>{}) ^ (val<FHIST_BITS>{fhist} >> hard<29>{});
        val<LOGFGEHL-LOGLINEINST> fidx = fgehl_base_index ^ val<LOGFGEHL-LOGLINEINST>{fh_mix};
        fidx.fanout(hard<LINEINST+1>{});
        fgehl_idx = fidx;
        fgehl_map = arr<val<PERCWIDTH,i64>,LINEINST>{[&](u64 offset){
            return fgehl[offset].read(fidx);
        }};
        fgehl_map.fanout(hard<2>{});
#endif
#ifdef DEBUG_ENERGY
        energy_mark("SCHistRead");
#endif

        threshold = arr<val<TOTAL_THREBITS>,LINEINST>{[&](u64 offset){
            return global_thre + thre1[offset];
        }};
        threshold.fanout(hard<4>{});

        //100 ps 2.5 cycle
        sc_sum = arr<val<TOTAL_THREBITS,i64>,LINEINST>{[&](u64 offset){
#if !defined(SC_USE_BIAS) && !defined(SC_USE_GEHL) && !defined(SC_FGEHL)
            static_cast<void>(offset);
#endif
#ifdef SC_USE_BIAS
            val<PERCWIDTH,i64> bias = bias_map[offset];
            bias.fanout(hard<2>{});
            val<1> bias_sign = val<1>{bias >> hard<PERCWIDTH-1>{}};
            val<TOTAL_THREBITS,i64> bias_ext = concat(bias_sign.replicate(hard<TOTAL_THREBITS-PERCWIDTH>{}).concat(), val<PERCWIDTH>{bias});
#else
            val<TOTAL_THREBITS,i64> bias_ext = val<TOTAL_THREBITS,i64>{0};
#endif
#ifdef SC_USE_GEHL
            val<PERCWIDTH,i64> gehl0 = gehl_map[0][offset];
            val<PERCWIDTH,i64> gehl1 = gehl_map[1][offset];
            gehl0.fanout(hard<2>{});
            gehl1.fanout(hard<2>{});
            val<1> gehl0_sign = val<1>{gehl0 >> hard<PERCWIDTH-1>{}};
            val<1> gehl1_sign = val<1>{gehl1 >> hard<PERCWIDTH-1>{}};
            val<TOTAL_THREBITS,i64> gehl0_ext = concat(gehl0_sign.replicate(hard<TOTAL_THREBITS-PERCWIDTH>{}).concat(), val<PERCWIDTH>{gehl0});
            val<TOTAL_THREBITS,i64> gehl1_ext = concat(gehl1_sign.replicate(hard<TOTAL_THREBITS-PERCWIDTH>{}).concat(), val<PERCWIDTH>{gehl1});
#else
            val<TOTAL_THREBITS,i64> gehl0_ext = val<TOTAL_THREBITS,i64>{0};
            val<TOTAL_THREBITS,i64> gehl1_ext = val<TOTAL_THREBITS,i64>{0};
#endif
#ifdef SC_FGEHL
            val<PERCWIDTH,i64> fgehl_v = fgehl_map[offset];
            fgehl_v.fanout(hard<2>{});
            val<1> fgehl_sign = val<1>{fgehl_v >> hard<PERCWIDTH-1>{}};
            val<TOTAL_THREBITS,i64> fgehl_ext = concat(fgehl_sign.replicate(hard<TOTAL_THREBITS-PERCWIDTH>{}).concat(), val<PERCWIDTH>{fgehl_v});
#else
            val<TOTAL_THREBITS,i64> fgehl_ext = val<TOTAL_THREBITS,i64>{0};
#endif
            val<PERCWIDTH+1,i64> gehl_pair_sum = gehl0_ext + gehl1_ext;
            val<PERCWIDTH+2,i64> sc_hist_sum = gehl_pair_sum + fgehl_ext;
            return sc_hist_sum + bias_ext;
        }};
        
        sc_sum.fanout(hard<4>{});
        //100 ps
        use_sc = arr<val<1>,LINEINST>{[&](u64 offset){
            val<1> prov_hit  = val<NUMG>{match1[offset]}!=val<NUMG>{0};
            val<1> is_prov_weak =  prov_weak[offset]; 
            val<1> is_prov_mid = prov_mid[offset].fo1();
            val<1> is_prov_high = prov_sat[offset].fo1();
            // return prov_hit & (
            //         ((is_prov_weak) & (sc_sum[offset]>(threshold[offset]>>3))));
            val<TOTAL_THREBITS> abs_sc = select(sc_sum[offset]>val<TOTAL_THREBITS,i64>{0}, val<TOTAL_THREBITS>{sc_sum[offset]}, val<TOTAL_THREBITS>{val<TOTAL_THREBITS,i64>{0}-sc_sum[offset]});
            abs_sc.fanout(hard<3>{});
            return prov_hit & (((is_prov_high) & (abs_sc>(threshold[offset]>>hard<1>{})))|
                    ((is_prov_mid) & (abs_sc>(threshold[offset]>>hard<2>{})))|
                    ((is_prov_weak) & (abs_sc>(threshold[offset]>>hard<3>{}))));
        }};
        use_sc.fanout(hard<2>{});
        sc_pred = arr<val<1>,LINEINST>{[&](u64 offset){
            return sc_sum[offset] > val<TOTAL_THREBITS,i64>{0};
        }};
        sc_pred.fanout(hard<2>{});
        //70ps
        sc_p2 = arr<val<1>,LINEINST> {[&](u64 offset){
            val<1> tage_pred = tage_p2.make_array(val<1>{})[offset];
            return select(use_sc[offset],sc_pred[offset],tage_pred);
        }}.concat();
#ifdef DEBUG_ENERGY
        energy_mark("SCSelect");
#endif
    }
#endif
    val<1> predict2(val<64> inst_pc)
    {
        tage_pred(inst_pc);
#ifdef GATE
    gating.fanout(hard<2>{});
    p1.fanout(hard<2>{});        
    #ifndef MY_SC
        p2 = select(gating,p1,tage_p2);
    #else
        sc_predict(inst_pc);
        p2 = select(gating,p1,sc_p2);
    #endif
#else
    #ifndef MY_SC
        p2 = tage_p2;
    #else
        sc_predict(inst_pc);
        p2 = sc_p2;
    #endif
#endif
        // p2 = tage_p2;
        p2.fanout(hard<LINEINST>{});
        val<1> taken = (block_entry & p2) != hard<0>{};
        taken.fanout(hard<2>{});

        reuse_prediction(~val<1>{block_entry>>(LINEINST-1)});
#ifdef DEBUG_ENERGY
        energy_mark("Predict2Final");
#endif
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = ((block_entry<<block_size) & p2) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{block_entry>>(LINEINST-1-block_size)});
        block_size++;
        return taken;
    }

    void update_condbr(val<64> branch_pc_in, val<1> taken, [[maybe_unused]] val<64> next_pc)
    {
        assert(num_branch < LINEINST);
        branch_offset[num_branch] = branch_pc_in >> 2;
        branch_dir[num_branch] = taken.fo1();
        branch_pc[num_branch] = branch_pc_in.fo1();
        num_branch++;
    }

#ifdef MY_SC
    void pipe_global_thre_delta(val<2,i64> new_delta, val<1> flush = val<1>{0})
    {
        flush.fanout(hard<SC_GLOBAL_THRE_PIPE_STAGES + 2>{});
        val<2,i64> apply_delta = select(flush, val<2,i64>{0}, val<2,i64>{global_thre_pipe[SC_GLOBAL_THRE_PIPE_STAGES-1]});
        apply_delta.fanout(hard<2>{});
        val<1> apply_upd = apply_delta != val<2,i64>{0};
        val<1> apply_inc = apply_delta > val<2,i64>{0};
        execute_if(apply_upd, [&](){
            global_thre = update_ctr(global_thre, apply_inc);
        });
        for (u64 i = SC_GLOBAL_THRE_PIPE_STAGES-1; i != 0; i--) {
            global_thre_pipe[i] = select(flush, val<2,i64>{0}, val<2,i64>{global_thre_pipe[i-1]});
        }
        global_thre_pipe[0] = select(flush, val<2,i64>{0}, new_delta);
    }

#endif

    void update_cycle(instruction_info &block_end_info)
    {
        val<1> &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc = block_end_info.next_pc;
        // updates for all conditional branches in the predicted block
        if (num_branch == 0) {
            for (u64 i=0; i<NUMG; i++) {
#if MY_BP_V1_GHYST_USE_WB_RWRAM
                ghyst[i].write(gindex[i],val<2>{0},val<1>{0},gindex[i],val<1>{1},val<1>{0});
#endif
            }
            // no conditional branch in this block
            val<1> line_end = block_entry >> (LINEINST-block_size);
            // update global history if previous block ended on a mispredicted not-taken branch
            // (we are still in the same line, this is the last chunk)
            // or if the block ends before the line boundary (unconditional jump)
            val<1> actual_block = ~(true_block & line_end.fo1());
            actual_block.fanout(hard<GHIST+NUMG*2+2>{});
            execute_if(actual_block, [&](){
                next_pc.fanout(hard<2>{});
                global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
                gfolds.update(val<PATHBITS>{next_pc>>2});
                true_block = 1;
            });
#ifdef MY_SC
            pipe_global_thre_delta(val<2,i64>{0}, val<1>{0});
#endif
#ifdef DEBUG_ENERGY
            energy_mark("UpdateNoBranch");
            energy_mark("CycleEnd");
#endif
            return; // stop here
        }
        //TODO: gate
        mispredict.fanout(hard<3>{});
        val<1> correct_pred = ~mispredict;
        correct_pred.fanout(hard<NUMG+2>{});
        index1.fanout(hard<LINEINST+2>{});
        p2.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST+2>{});
        gindex.fanout(hard<4>{});
        htag.fanout(hard<3>{});
        readb.fanout(hard<2>{});
        readt.fanout(hard<4>{});
        readc.fanout(hard<2>{});
        readh.fanout(hard<2>{});
        //SC fix
        match1.fanout(hard<3>{});
        
        match2.fanout(hard<2>{});
        pred1.fanout(hard<2>{});
        pred2.fanout(hard<2+NUMG>{});
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<2>{});
        gfolds.fanout(hard<2>{});
#ifdef USE_ALT
        meta.fanout(hard<2>{});
#endif
        val<LOGLINEINST> last_offset = branch_offset[num_branch-1];
        last_offset.fanout(hard<4*NUMG+2>{});

        u64 update_valid = (u64(1)<<num_branch)-1;
        arr<val<LINEINST>,LINEINST> update_mask = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){return branch_offset[i] == offset;};
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});

        arr<val<1>,LINEINST> is_branch = [&](u64 offset){
            return update_mask[offset] != hard<0>{};
        };
        is_branch.fanout(hard<11>{});

        val<LINEINST> branch_mask = is_branch.concat();
        branch_mask.fanout(hard<2>{});

        val<LINEINST> actualdirs = branch_dir.concat();
        actualdirs.fanout(hard<LINEINST>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_table_reads_hits(is_branch);
#endif
#endif

        arr<val<1>,LINEINST> branch_taken = [&](u64 offset){
            return (actualdirs & update_mask[offset]) != hard<0>{};
        };

        //TODO: fix it
        branch_taken.fanout(hard<7>{});

        arr<val<NUMG+1>,LINEINST> actual_match1 = [&] (u64 offset) {
            return select(is_branch[offset],match1[offset],val<NUMG+1>{0});
        };
        actual_match1.fanout(hard<2>{});

        val<NUMG> primary_mask = actual_match1.fold_or();
        primary_mask.fanout(hard<2>{});
        arr<val<1>,NUMG> primary = primary_mask.make_array(val<1>{});
        primary.fanout(hard<3>{});

        arr<val<1>,LINEINST> primary_wrong = [&](u64 offset){
            return pred1[offset] != branch_taken[offset];
        };
        primary_wrong.fanout(hard<2>{});
        arr<val<1>,LINEINST> bhyst_bim_primary = [&](u64 offset){ return actual_match1[offset] >> NUMG; };

        // select some candidate entries for allocation
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();

        arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == concat(last_offset,htag[i]);};

        val<NUMG+1> last_match1 = last_tagcmp.fo1().append(1).concat().one_hot();
        last_match1.fanout(hard<2>{});
        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match1-1);
        postmask.fanout(hard<2>{});
        val<NUMG> candallocmask = postmask & notumask; // candidate post entries for allocation
        candallocmask.fanout(hard<2>{});
        // if multiple candidate entries, we select a single one, with some randomization
        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask1 = collamask.one_hot();
        collamask1.fanout(hard<3>{});
        val<NUMG> collamask2 = (collamask^collamask1).one_hot();
        val<NUMG> collamask12 = select(val<2>{std::rand()}==hard<0>{}, collamask2.fo1(), collamask1);
        arr<val<1>,NUMG> allocate = collamask12.fo1().reverse().make_array(val<1>{});
        allocate.fanout(hard<7>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_table_alloc(allocate);
#endif
#endif

        // associate a branch direction to each global table
        arr<val<1>,NUMG> bdir = [&](u64 i) {
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            val<LOGLINEINST> offset = select(allocate[i],last_offset,tag_offset.fo1());
            offset.fanout(hard<LINEINST>{});
            arr<val<1>,LINEINST> match_offset = [&](u64 j){return branch_offset[j] == offset;};
            return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
        };

        bdir.fanout(hard<2>{});

        // tell if global prediction is incorrect
        arr<val<1>,NUMG> badpred1 = [&](u64 i){
            return readc[i] != bdir[i];
        };
        badpred1.fanout(hard<3>{});

        // associate to each global table a bit telling if local prediction differs from secondary prediction
        arr<val<1>,NUMG> altdiffer = [&](u64 i){
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return readc[i] != pred2.select(tag_offset.fo1());
        };

        // associate to each global table a bit telling if prediction for owning branch is correct
        arr<val<1>,NUMG> goodpred = [&](u64 i){
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return (tag_offset.fo1() != last_offset) | correct_pred;
        };

        // do P1 and P2 agree?
        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask;
        //TODO:gate
        disagree_mask.fanout(hard<2+2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});

        // read the P1 hysteresis if P1 and P2 disagree
        arr<val<1>,LINEINST> p1_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff disagreement and hysteresis is weak
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(index1); // hyst=0 means weak
            });
        };



        //TODO:gate
        // read the bimodal hysteresis if bimodal caused a misprediction
        arr<val<1>,LINEINST> b_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff cause of misprediction and hysteresis is weak
#ifdef GATE
            return execute_if(bhyst_bim_primary[offset].fo1() & primary_wrong[offset] & (~gating), [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
#else
            return execute_if(bhyst_bim_primary[offset].fo1() & primary_wrong[offset], [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
#endif
        };

        // determine which primary global predictions are incorrect with a weak hysteresis
        arr<val<1>,NUMG> g_weak = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & badpred1[i] & (readh[i]==hard<0>{});
        };
#ifdef GATE
        arr<val<1>,NUMG> g_sat = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & (readh[i]==hard<3>{});
        };
#endif
        // need extra cycle for modifying prediction bits and for TAGE allocation
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        // some_badpred1 is consumed by extra_cycle and PERF counters in this mode.
        some_badpred1.fanout(hard<2>{});
#endif
#endif

#ifdef MY_SC
        sc_sum.fanout(hard<3>{});
        // use_sc.fanout(hard<LINEINST*3>{});
        // match1.fanout(hard<LINEINST*2>{});
        // branch_pc.fanout(hard<LINEINST*2>{});
        // threshold.fanout(hard<LINEINST>{});

        // compute intermediate signals once, reuse below
        sc_wrong_arr = arr<val<1>,LINEINST>{[&](u64 offset){
            val<1> sc_pred_bit = sc_pred[offset];
            return sc_pred_bit != branch_taken[offset];
        }};
        sc_wrong_arr.fanout(hard<3>{});
        
        
        
        do_update_arr = arr<val<1>,LINEINST>{[&](u64 offset){
            return sc_wrong_arr[offset] | ~use_sc[offset];
        }};
        do_update_arr.fanout(hard<6>{});
        prov_hit_arr = arr<val<1>,LINEINST>{[&](u64 offset){
            return val<NUMG>{match1[offset]} != val<NUMG>{0};
        }};
        prov_hit_arr.fanout(hard<6>{});
        thre_guard_arr = arr<val<1>,LINEINST>{[&](u64 offset){
            val<TOTAL_THREBITS> abs_sc = select(sc_sum[offset]>val<TOTAL_THREBITS,i64>{0}, val<TOTAL_THREBITS>{sc_sum[offset]}, val<TOTAL_THREBITS>{val<TOTAL_THREBITS,i64>{0}-sc_sum[offset]});
            return abs_sc > (threshold[offset]>>hard<1>{});
            // return val<1>{1};
        }};

        val<1> sc_need_update = arr<val<1>,LINEINST>{[&](u64 offset){
            return is_branch[offset] & do_update_arr[offset] & prov_hit_arr[offset];
        }}.fold_or();

        val<1> extra_cycle = some_badpred1 | mispredict | (disagree_mask != hard<0>{}) | sc_need_update;
#endif

        extra_cycle.fanout(hard<NUMG*2+8>{});
        need_extra_cycle(extra_cycle);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        val<1> p1_update = disagree_mask != hard<0>{};
#ifdef MY_SC
        perf_count_extra_cycle(extra_cycle, some_badpred1, mispredict, p1_update, sc_need_update);
#else
        perf_count_extra_cycle(extra_cycle, some_badpred1, mispredict, p1_update);
#endif
#endif
#endif
#ifdef DEBUG_ENERGY
        energy_mark("UpdatePrep");
#endif

#ifdef USE_ALT
        // update meta counter
        arr<val<1>,LINEINST> altdiff = [&](u64 offset){
            // for each offset, tell if primary and secondary predictions differ
            return (match2[offset] != hard<0>{}) & (pred2[offset] != pred1[offset]);
        };
        arr<val<2,i64>,LINEINST> meta_incr = [&](u64 offset) -> val<2,i64> {
            val<1> update_meta = is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
            val<1> bad_pred2 = (pred2[offset] != branch_taken[offset]);
            return select(update_meta.fo1(),concat(bad_pred2.fo1(),val<1>{1}),val<2>{0});
        };
        for (u64 i=METAPIPE-1; i!=0; i--) {
            meta[i] = meta[i-1];
        }
        auto newmeta = meta[0] + meta_incr.fo1().fold_add();
        newmeta.fanout(hard<3>{});
        using meta_t = valt<decltype(meta[0])>;
        meta[0] = select(newmeta>meta_t::maxval, meta_t{meta_t::maxval}, select(newmeta<meta_t::minval, meta_t{meta_t::minval}, meta_t{newmeta}));
#endif


        // overwrite the tag in the allocated entry (mispredict)
        for (u64 i=0; i<NUMG; i++) {
            execute_if(allocate[i], [&](){gtag[i].write(gindex[i],concat(last_offset,htag[i]));});
        }


        // update the u bits
        arr<val<1>,NUMG> update_u = [&](u64 i){
            return primary[i] & altdiffer[i].fo1();
        };
        // if all post entries have the u bit set, reset their u bits
        val<1> noalloc = (candallocmask == hard<0>{});
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        // noalloc is consumed by uclear path and end-of-cycle perf accounting.
        noalloc.fanout(hard<2>{});
#endif
#endif
        val<NUMG> uclearmask = postmask & noalloc.replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<2>{});

        for (u64 i=0; i<NUMG; i++) {
            execute_if(update_u[i].fo1() | allocate[i] | uclear[i], [&]() {
                val<1> newu = goodpred[i].fo1() & ~allocate[i] & ~uclear[i];
                ubit[i].write(gindex[i],newu.fo1(),extra_cycle);
            });
        }
#ifdef DEBUG_ENERGY
        energy_mark("UpdateAllocUBit");
#endif

        // update P1 prediction if P1 and P2 disagree and the hysteresis bit is weak
        auto p2_split = p2.make_array(val<1>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(p1_weak[offset].fo1(), [&](){
                // update with the P2 prediction, not with the actual branch direction
                table1_pred[offset].write(index1,p2_split[offset].fo1());
            });
        }
        // update P1 hysteresis
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(is_branch[offset],[&](){
                table1_hyst[offset].write(index1,~disagree[offset]);
            });
        }

        // update incorrect bimodal prediction if primary provider and hysteresis is weak
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(b_weak[offset].fo1(), [&](){
                bim[offset].write(bindex,branch_taken[offset]);
            });
        }

        // update bimodal hysteresis if bimodal is primary provider
        for (u64 offset=0; offset<LINEINST; offset++) {
            val<1> bim_primary = match1[offset] >> NUMG;
            execute_if(is_branch[offset] & bim_primary.fo1(), [&](){
                bhyst[offset].write(bindex,~primary_wrong[offset]);
            });
        }
#ifdef DEBUG_ENERGY
        energy_mark("UpdateLocalTables");
#endif

        // update incorrect global prediction if primary provider and the hysteresis is weak;
        // initialize global prediction in the allocated entry
        for (u64 i=0; i<NUMG; i++) {
            execute_if(g_weak[i].fo1() | allocate[i], [&](){
                gpred[i].write(gindex[i],bdir[i]);
            });
        }

        // update global prediction hysteresis if primary provider or allocated entry
#ifdef GATE
        val<1> ghyst_read_en = ~gating;
#else
        val<1> ghyst_read_en = 1;
#endif
        ghyst_read_en.fanout(hard<2*NUMG + 1>{});

        for (u64 i=0; i<NUMG; i++) {
            val<1> do_ghyst_write = primary[i] | allocate[i];
            // if allocated entry, set hysteresis to 0;
            // otherwise, increment hysteresis if correct pred, decrement if incorrect
            val<2> newhyst = select(allocate[i],val<2>{0},update_ctr(readh[i],~badpred1[i]));
#if MY_BP_V1_GHYST_USE_WB_RWRAM
  #if MY_BP_V1_GHYST_HIT_UPDATE
            ghyst[i].write_update(gindex[i], newhyst.fo1(), do_ghyst_write,
                                  primary[i] & ~allocate[i], ~badpred1[i],
                                  gindex[i], ~extra_cycle, extra_cycle);
#else
            ghyst[i].write(gindex[i],newhyst.fo1(),do_ghyst_write,gindex[i],~extra_cycle,extra_cycle);
#endif
#else
            ghyst[i].write(gindex[i],newhyst.fo1(),do_ghyst_write,extra_cycle);
#endif
        }

#ifdef RESET_UBITS
        uctr.fanout(hard<3>{});
        val<NUMG> allocmask1  = collamask1.reverse();
        allocmask1.fanout(hard<2>{});
        val<1> faralloc = (((last_match1>>3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
        val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
        uctrsat.fanout(hard<2>{});
        uctr = select(correct_pred,uctr,select(uctrsat,val<decltype(uctr)::size>{0},update_ctr(uctr,faralloc.fo1())));
        execute_if(uctrsat,[&](){for (auto &uram : ubit) uram.reset();});
#endif
#ifdef DEBUG_ENERGY
        energy_mark("UpdateGlobalTables");
#endif

#ifdef GATE

        gate_hit.fanout(hard<3>{});
        
        val<1> gate_alloc = (disagree_mask == hard<0>{}) & (~mispredict) & (~gating) & (~gate_hit) & g_sat.fold_or();
        val<1> gate_inc = (disagree_mask == hard<0>{}) & (~mispredict) & gate_hit &  g_sat.fold_or();
        val<1> gate_dec = mispredict & gate_hit;
        
        gate_idx.fanout(hard<(1<<LOGGATEBITS)>{});
        line_tag.fanout(hard<(1<<LOGGATEBITS)>{});
        gate_alloc.fanout(hard<(1<<LOGGATEBITS)>{});
        gate_inc.fanout(hard<(1<<LOGGATEBITS)*2>{});
        gate_dec.fanout(hard<(1<<LOGGATEBITS)>{});
        gate_ctr[0].fanout(hard<4>{});
        gate_ctr[1].fanout(hard<2>{});
        // gate_ctr[2].fanout(hard<3>{});
        // gate_inc
        for (u64 i=0; i<(1<<LOGGATEBITS); i++) {
            val<1> tag_match = gate_idx == val<LOGGATEBITS>{i};
            tag_match.fanout(hard<2>{});
            //one pipe delay
            gate_tag[2][i] = gate_tag[1][i];
            gate_ctr[2][i] = gate_ctr[1][i];
            gate_tag[1][i] = gate_tag[0][i];
            gate_ctr[1][i] = gate_ctr[0][i];
            execute_if(tag_match & gate_alloc, [&](){
                gate_tag[0][i] = line_tag;
            });
            execute_if(tag_match & (gate_inc | gate_dec), [&](){
                gate_ctr[0][i] = select(gate_inc, update_ctr(gate_ctr[0][i], val<1>{1}), gate_ctr[0][i]>>1);
            });
        }
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_gate(gating, mispredict, gate_inc, gate_dec);
#endif
#endif
#endif

        // update global history
        val<1> line_end = block_entry >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            next_pc.fanout(hard<2>{});
            global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
            gfolds.update(val<PATHBITS>{next_pc>>2});
        });

#if defined(MY_SC) && defined(SC_FGEHL)
        // FGEHL history: use existing framework signals only (last resolved branch + block next_pc)
        val<64> last_branch_pc = branch_pc[num_branch-1];
        val<1> last_branch_taken = branch_dir[num_branch-1];
        val<1> last_forward = next_pc > last_branch_pc;
        val<FHIST_BITS> next_pc_fold = val<FHIST_BITS>{next_pc >> 2};
        val<FHIST_BITS> branch_pc_fold = val<FHIST_BITS>{last_branch_pc >> 1};
        val<FHIST_BITS> new_fhist = (val<FHIST_BITS>{fhist} << hard<3>{}) ^ next_pc_fold ^ branch_pc_fold;
        execute_if(last_branch_taken & last_forward, [&](){
            fhist = new_fhist;
        });
#endif
#ifdef DEBUG_ENERGY
        energy_mark("UpdateHistory");
#endif


#ifdef MY_SC


        //global_thre: majority vote across all qualifying offsets
        arr<val<1>,LINEINST> thre_update_en = [&](u64 offset){
            return is_branch[offset] & do_update_arr[offset] & thre_guard_arr[offset];
        };
        thre_update_en.fanout(hard<2>{});
        // Meta-style signed accumulation:
        // +1 when SC is wrong, -1 when SC is right, 0 when no update.
        arr<val<2,i64>,LINEINST> gthre_incr = [&](u64 offset) -> val<2,i64> {
            val<1> update_en = thre_update_en[offset];
            val<1> dec_sign = ~sc_wrong_arr[offset]; // 0:+1, 1:-1
            return select(update_en.fo1(), concat(dec_sign.fo1(), val<1>{1}), val<2>{0});
        };
        auto gthre_vote = gthre_incr.fo1().fold_add();
        val<1> global_thre_upd = gthre_vote != val<decltype(gthre_vote)::size, i64>{0};
        val<1> global_thre_inc = gthre_vote > val<decltype(gthre_vote)::size, i64>{0};
        val<1> global_thre_decsign = ~global_thre_inc;
        global_thre_inc.fanout(hard<2>{});
        val<2,i64> gthre_delta = select(global_thre_upd.fo1(), concat(global_thre_decsign.fo1(), val<1>{1}), val<2>{0});
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        if (static_cast<bool>(global_thre_upd)) {
            perf_global_thre_update++;
            if (static_cast<bool>(global_thre_inc))
                perf_global_thre_inc++;
            else
                perf_global_thre_dec++;
        }
#endif
#endif
        pipe_global_thre_delta(gthre_delta, mispredict);

        
#ifdef SC_USE_BIAS
        bias_high_idx.fanout(hard<4>{});
#endif


        //per-offset updates: bias_pc and thre1 (no write conflicts)
        for (u64 offset = 0; offset < LINEINST; offset++) {
#ifdef SC_USE_BIAS
            val<LOGBIAS-LOGLINEINST-2> high_idx = bias_high_idx;
            val<2> write_tage_info = tage_info[offset];
#endif
            val<PRE_PC_THREBITS> old_thre1 = thre1[offset];
            val<PRE_PC_THREBITS> new_thre1 = update_ctr(old_thre1, sc_wrong_arr[offset]);
#if defined(SC_USE_BIAS) || defined(SC_USE_GEHL) || defined(SC_FGEHL)
            [[maybe_unused]] val<1> sc_update_en = is_branch[offset] & do_update_arr[offset] & prov_hit_arr[offset];
#endif
#ifdef SC_USE_BIAS
            write_tage_info.fanout(hard<4>{});
            high_idx.fanout(hard<4>{});
            execute_if(sc_update_en, [&](){
                arr<val<PERCWIDTH,i64>,4> new_bias_vec = arr<val<PERCWIDTH,i64>,4>{[&](u64 bank){
                    val<PERCWIDTH,i64> old_lane = val<PERCWIDTH,i64>{bias_bank_map[bank][offset]};
                    val<PERCWIDTH,i64> new_lane = update_ctr(old_lane, branch_taken[offset]);
                    return select(write_tage_info==val<2>{bank}, new_lane, old_lane);
                }};
                bias_pc[offset].write(high_idx, new_bias_vec);
            });
#endif
#ifdef SC_USE_GEHL
            for (u64 k = 0; k < NUMGEHL; k++) {
                val<LOGGEHL-LOGLINEINST> gehl_write_idx = gehl_idx[k];
                val<PERCWIDTH,i64> old_gehl = gehl_map[k][offset];
                gehl_write_idx.fanout(hard<2>{});
                old_gehl.fanout(hard<2>{});
                execute_if(sc_update_en, [&](){
                    gehl[k][offset].write(gehl_write_idx, update_ctr(old_gehl, branch_taken[offset]));
                });
            }
#endif
#ifdef SC_FGEHL
            val<LOGFGEHL-LOGLINEINST> fgehl_write_idx = fgehl_idx;
            val<PERCWIDTH,i64> old_fgehl = fgehl_map[offset];
            fgehl_write_idx.fanout(hard<2>{});
            old_fgehl.fanout(hard<2>{});
            execute_if(sc_update_en, [&](){
                fgehl[offset].write(fgehl_write_idx, update_ctr(old_fgehl, branch_taken[offset]));
            });
#endif
            thre1[offset] = select(thre_update_en[offset],new_thre1,thre1[offset]);
        }
#ifdef DEBUG_ENERGY
        energy_mark("UpdateSCThre1");
#ifdef SC_USE_BIAS
        energy_mark("UpdateSCBias");
#endif
#ifdef SC_USE_GEHL
        energy_mark("UpdateSCGehl");
#endif
#ifdef SC_FGEHL
        energy_mark("UpdateSCFgehl");
#endif
#endif
#endif

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_end_of_cycle(is_branch, branch_taken, mispredict, allocate, postmask, noalloc);
#endif
#endif


#ifdef DEBUG_ENERGY
        energy_mark("CycleEnd");
#endif
        num_branch = 0; // done
    }
};
