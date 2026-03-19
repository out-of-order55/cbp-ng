// this is a basic TAGE, not necessarily well optimized

#define USE_META
#define RESET_UBITS
// #define BANK
// #define GATE
#define MY_SC
// #define HASH_TAG
#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <sqlite3.h>

using namespace hcm;
#ifdef DEBUG_ENERGY
    struct energy_monitor monitor;
#endif
template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=12, u64 LOGB=12, u64 TAGW=12, u64 GHIST=100, u64 LOGP1=14, u64 GHIST1=6,u64 LOGBANKS = 2,u64 LOGBIAS = 11>
struct my_bp_v1 : predictor {
    // provides 2^(LOGLB-2) predictions per cycle
    // P2 is a TAGE, P1 is a gshare
    static_assert(LOGLB>2);
    static_assert(NUMG>0);

    static constexpr u64 PERCWIDTH = 6;
    static constexpr u64 TOTAL_THREBITS = 9;
    static constexpr u64 GLOBAL_THREBITS = 7;
    static constexpr u64 PRE_PC_THREBITS = 7;
    // static constexpr u64 LOGTHREBITS = 5;

    static constexpr u64 MINHIST = 2;
    static constexpr u64 METABITS = 4;
    static constexpr u64 UCTRBITS = 8;
    static constexpr u64 PATHBITS = 6;

    static constexpr u64 LOGGATEBITS = 10;
    static constexpr u64 TAGGATEBITS = 6;

#ifdef USE_META
    static constexpr u64 METAPIPE = 2;
#endif
    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGB > LOGLINEINST);
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 bindex_bits = LOGB-LOGLINEINST;
    static_assert(TAGW > LOGLINEINST); // the unhashed line offset is part of the tag
#ifdef  HASH_TAG
    static constexpr u64 HTAGBITS = TAGW; // hashed tag bits
#else
    static constexpr u64 HTAGBITS = TAGW-LOGLINEINST; // hashed tag bits
#endif
    geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;
    reg<1> true_block = 1;

    // for P1
    reg<GHIST1> global_history1;
    reg<index1_bits> index1;
    arr<reg<1>,LINEINST> readp1; // prediction bits read from P1 table for each offset
    reg<LINEINST> p1; // P1 predictions

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

    arr<reg<NUMG+1>,LINEINST> match; // all matches for each offset
    arr<reg<NUMG+1>,LINEINST> match1; // longest match for each offset
    arr<reg<NUMG+1>,LINEINST> match2; // second longest match for each offset
    arr<reg<1>,LINEINST> prov_weak;
    arr<reg<1>,LINEINST> prov_mid; 
    arr<reg<1>,LINEINST> prov_sat;  

    arr<reg<1>,LINEINST> pred1; // primary P2 prediction for each offset
    arr<reg<1>,LINEINST> pred2; // alternate P2 prediction for each offset
    arr<reg<1>,LINEINST> use_sc; // alternate P2 prediction for each offset
    reg<LINEINST> p2; // final P2 predictions
    reg<LINEINST> tage_p2; // final P2 predictions
    reg<LINEINST> sc_p2; // final P2 predictions
    arr<reg<1>,LINEINST> sc_pred;
    

#ifdef USE_META
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
    u64 perf_alt_used[NUMG] = {};
    u64 perf_alt_correct[NUMG] = {};

    // Bimodal tracking
    u64 perf_bimodal_used = 0;
    u64 perf_bimodal_correct = 0;

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
#ifdef MY_SC
    u64 perf_extra_cycle_sc_update = 0; // SC needs update
    // SC prediction source counters
    u64 perf_sc_override = 0;        // SC overrode TAGE (p2 != tage_p2)
    u64 perf_sc_override_correct = 0; // SC override was correct
    // threshold update counters
    u64 perf_thre_update = 0;        // thre1 updated
    u64 perf_thre_update_inc = 0;    // thre1 incremented (sc wrong)
    u64 perf_thre_update_dec = 0;    // thre1 decremented (sc correct)
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
        for (u64 j=0; j<NUMG; j++) {
            total_prov += perf_provider_used[j];
            total_prov_ok += perf_provider_correct[j];
            total_alt += perf_alt_used[j];
            total_alt_ok += perf_alt_correct[j];
        }

        u64 total_all = total_prov + total_alt + perf_bimodal_used;

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
        std::cerr << "Threshold updates: total=" << perf_thre_update
                  << "  inc=" << perf_thre_update_inc
                  << "  dec=" << perf_thre_update_dec << "\n";
#endif

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
#endif

#ifdef GATE
    rwram<TAGGATEBITS,(1<<LOGGATEBITS),4> gate {"gate"};
    reg<1> gating;
    reg<TAGGATEBITS> atag[2];
    // reg<LOGGATEBITS> aidx[2];
#endif  
    
#ifdef BANK
    // P1 (gshare)
    rwram<1,(1<<index1_bits),1<<LOGBANKS> table1_pred[LINEINST] {"P1 pred"}; // P1 prediction bit

    // P2 (TAGE)
    rwram<TAGW,(1<<LOGG),1<<LOGBANKS> gtag[NUMG] {"tags"}; // tags
    // rwram<TAGW,(1<<LOGG),4> gtag[NUMG] {"tags"}; // tags
    rwram<1,(1<<LOGG),1<<LOGBANKS> gpred[NUMG] {"gpred"}; // predictions
    rwram<2,(1<<LOGG),1<<LOGBANKS> ghyst[NUMG] {"ghyst"}; // hysteresis

    rwram<1,(1<<LOGG),1<<LOGBANKS> ubit[NUMG] {"u"}; // "useful" bits
    rwram<1,(1<<bindex_bits),1<<LOGBANKS> bim[LINEINST] {"bpred"}; // bimodal prediction bits

    zone UPDATE_ONLY;
    rwram<1,(1<<index1_bits),1<<LOGBANKS> table1_hyst[LINEINST] {"P1 hyst"}; // P1 hysteresis
    rwram<1,(1<<bindex_bits),1<<LOGBANKS> bhyst[LINEINST] {"bhyst"}; // bimodal hysteresis
#else
    // P1 (gshare)
    ram<val<1>,(1<<index1_bits)> table1_pred[LINEINST] {"P1 pred"}; // P1 prediction bit

    // P2 (TAGE)
    ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"tags"}; // tags
    // rwram<TAGW,(1<<LOGG),4> gtag[NUMG] {"tags"}; // tags
    ram<val<1>,(1<<LOGG)> gpred[NUMG] {"gpred"}; // predictions
    rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"}; // hysteresis

    rwram<1,(1<<LOGG),4> ubit[NUMG] {"u"}; // "useful" bits
    rwram<1,(1<<bindex_bits),4> bim[LINEINST] {"bpred"}; // bimodal prediction bits

    zone UPDATE_ONLY;
    ram<val<1>,(1<<index1_bits)> table1_hyst[LINEINST] {"P1 hyst"}; // P1 hysteresis
    ram<val<1>,(1<<bindex_bits)> bhyst[LINEINST] {"bhyst"}; // bimodal hysteresis
#endif
    #ifdef MY_SC
    //bias_pc is concat pc and tage
    //idx = (pc>>(LOGLB+2)) 2bit tage info(prov_weak,prov_taken)
    arr<reg<TOTAL_THREBITS>,LINEINST> threshold;
    ram<val<PERCWIDTH,i64>,(1<<(LOGBIAS-LOGLINEINST-2))> bias_pc[4][LINEINST] {"Bias pc"};
    arr<reg<PERCWIDTH,i64>,4> bias_lmap;
    arr<reg<PERCWIDTH,i64>,LINEINST> bias_map;
    arr<reg<2>,LINEINST> tage_info;
    reg<LOGBIAS-LOGLINEINST-2> bias_high_idx;

    arr<reg<PRE_PC_THREBITS>,LINEINST> thre1;
    reg<GLOBAL_THREBITS> global_thre;
    arr<reg<TOTAL_THREBITS,i64>,LINEINST> sc_sum;
    // intermediate update signals
    arr<reg<1>,LINEINST> sc_wrong_arr;
    arr<reg<1>,LINEINST> do_update_arr;
    arr<reg<1>,LINEINST> prov_hit_arr;
    arr<reg<1>,LINEINST> thre_guard_arr;
    // rwram<PERCWIDTH,(1<<LOGBIAS)/NUMBANKS,NUMBANKS> bias_pc {"Bias pc"};
    #endif

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
        // global_thre = val<GLOBAL_THREBITS>{23};
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
    
        inst_pc.fanout(hard<2>{});
        new_block(inst_pc);
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
        return (block_entry & p1) != hard<0>{};
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        return ((block_entry<<block_size) & p1) != hard<0>{};
    };
#ifdef GATE
    void gate_predict2([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        val<LOGGATEBITS> lineaddr = inst_pc >> LOGLB;
        val<TAGGATEBITS> tag = inst_pc >> (LOGLB+LOGGATEBITS);

        // aidx[0] = lineaddr;
        val<TAGGATEBITS> read_tag = gate.read(lineaddr);
        atag[1] = atag[0];
        atag[0] = read_tag;

        gating = (atag[1] == tag);
    };
#endif
    void tage_pred(val<64> inst_pc){
        // inst_pc.fanout(hard<2>{});
#ifdef GATE
        inst_pc.fanout(hard<2>{});
        gate_predict2(inst_pc);
#endif
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
#ifdef DEBUG_ENERGY
        energy_checkpoint(monitor);
#endif
        for (u64 i=0; i<NUMG; i++) {
            readt[i] = gtag[i].read(gindex[i]);
            readc[i] = gpred[i].read(gindex[i]);
            readh[i] = ghyst[i].read(gindex[i]);
            readu[i] = ubit[i].read(gindex[i]);
        }
#ifdef DEBUG_ENERGY
        energy_checkpoint(monitor);
#endif
#endif
        readt.fanout(hard<LINEINST+1>{});
        readc.fanout(hard<3>{});
#ifdef MY_SC
        readh.fanout(hard<1+4*LINEINST>{});
#endif
        readu.fanout(hard<2>{});
        notumask = ~readu.concat();
        notumask.fanout(hard<2>{});

        // gather prediction bits for each offset
        val<NUMG> gpreds = readc.concat();
        gpreds.fanout(hard<LINEINST>{});
        arr<val<NUMG+1>,LINEINST> preds = [&](u64 offset){return concat(readb[offset],gpreds);};
        preds.fanout(hard<2*LINEINST>{});
#ifdef HASH_TAG
        // generate match mask for each offset
        static_loop<LINEINST>([&]<u64 offset>(){
            arr<val<1>,NUMG> tagcmp = [&](int i){return readt[i] == (htag[i]^val<HTAGBITS>{offset});};
            match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() ); // bimodal is default when no match
        });
#else
        // hashed tags comparisons
        arr<val<1>,NUMG> htagcmp_split = [&](int i){return val<HTAGBITS>{readt[i]} == htag[i];};
        val<NUMG> htagcmp = htagcmp_split.fo1().concat();
        htagcmp.fanout(hard<LINEINST>{});

        // generate match mask for each offset
        static_loop<LINEINST>([&]<u64 offset>(){
            arr<val<1>,NUMG> tagcmp = [&](int i){return val<LOGLINEINST>{readt[i]>>HTAGBITS} == hard<offset>{};};
            match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp); // bimodal is default when no match

        });
#endif
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
        pred1.fanout(hard<2>{});

        // for each offset, find second longest match and select secondary prediction
        for (u64 offset=0; offset<LINEINST; offset++) {
            match2[offset] = (match[offset]^match1[offset]).one_hot();
        }
        match2.fanout(hard<2>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
        }
        pred2.fanout(hard<2>{});


#ifdef USE_META
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
        // Store prediction source for each offset (used in update_cycle)
        for (u64 offset=0; offset<LINEINST; offset++) {
#ifdef USE_META
            val<1> use_alt_o = metasign & newly_alloc[offset] & (match2[offset]!=hard<0>{});
            // bimodal is bit NUMG (highest bit), tables are bits 0..NUMG-1
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
#endif
#endif
    }
#ifdef MY_SC
    void sc_predict(val<64> inst_pc){

        // inst_pc.fanout(hard<5>{});
        match1.fanout(hard<2>{});
        pred1.fanout(hard<2>{});
        prov_weak.fanout(hard<2>{});
        //bias
        val<LOGBIAS-LOGLINEINST-2> bias_pc_high_index =  inst_pc>>(LOGLB+2);
        bias_pc_high_index.fanout(hard<(LINEINST*4+1)>{});
#ifdef DEBUG_ENERGY
        energy_checkpoint(monitor);
#endif
        bias_map = arr<val<PERCWIDTH>,LINEINST>{[&](u64 offset){
            arr<val<PERCWIDTH>,4> bank_out = [&](u64 b){ return bias_pc[b][offset].read(bias_pc_high_index); };
            val<1> prov_hit  = val<NUMG>{match1[offset]}!=val<NUMG>{0};
            val<1> prov_pred = pred1[offset] & prov_hit;
            val<1> is_prov_weak =  prov_weak[offset];
            tage_info[offset] = concat(prov_pred,is_prov_weak);
            val<PERCWIDTH> result = bank_out.select(tage_info[offset]);
            return result;
        }};
#ifdef DEBUG_ENERGY
        energy_checkpoint(monitor);
#endif
        bias_map.fanout(hard<2>{});
        bias_high_idx = bias_pc_high_index;
        // val<LOGTHREBITS> base_idx1 = concat(val<LOGTHREBITS-LOGLINEINST>{inst_pc>>LOGLB},val<LOGLINEINST>{0}) ^ (inst_pc>>(2+2)); 
        // val<LOGTHREBITS> base_idx2 = concat(val<LOGTHREBITS-LOGLINEINST>{inst_pc>>LOGLB},val<LOGLINEINST>{0}) ^ (inst_pc>>(2+5));
        
        // base_idx1.fanout(hard<LINEINST>{});
        // base_idx2.fanout(hard<LINEINST>{});
        threshold = arr<val<TOTAL_THREBITS>,LINEINST>{[&](u64 offset){
            // val<LOGTHREBITS> idx1 = base_idx1 ^ val<LOGTHREBITS>{offset};
            // val<LOGTHREBITS> idx2 = base_idx2 ^ val<LOGTHREBITS>{offset};
            return global_thre + thre1[offset];
        }};
        threshold.fanout(hard<5>{});

        //100 ps 2.5 cycle
        sc_sum = arr<val<TOTAL_THREBITS,i64>,LINEINST>{[&](u64 offset){
            val<PERCWIDTH> bias = bias_map[offset];
            bias.fanout(hard<2>{});
            val<1> sign_bit = val<1>{bias >> hard<PERCWIDTH-1>{}};
            // sign_bit.fanout(hard<TOTAL_THREBITS-PERCWIDTH>{});
            val<TOTAL_THREBITS,i64> res = concat(sign_bit.replicate(hard<TOTAL_THREBITS-PERCWIDTH>{}).concat(), val<PERCWIDTH>{bias});
            return res;
        }};
        
        sc_sum.fanout(hard<5>{});
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
        // sc_p2.fanout(hard<2>{});
        
    }
#endif
    val<1> predict2(val<64> inst_pc)
    {
        tage_pred(inst_pc);
        
    #ifndef MY_SC
        p2 = tage_p2;
    #else
        sc_predict(inst_pc);
        p2 = sc_p2;
    #endif
        // p2 = tage_p2;
        p2.fanout(hard<LINEINST>{});
        val<1> taken = (block_entry & p2) != hard<0>{};
        taken.fanout(hard<2>{});

        reuse_prediction(~val<1>{block_entry>>(LINEINST-1)});
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

    void update_cycle(instruction_info &block_end_info)
    {
        val<1> &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc = block_end_info.next_pc;
        // updates for all conditional branches in the predicted block
        if (num_branch == 0) {
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
            return; // stop here
        }
        mispredict.fanout(hard<NUMG+2>{});
        val<1> correct_pred = ~mispredict;
        correct_pred.fanout(hard<NUMG+2>{});
        index1.fanout(hard<LINEINST*3>{});
        p2.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST*3>{});
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
#ifdef USE_META
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

        val<LINEINST> actualdirs = branch_dir.concat();
        actualdirs.fanout(hard<LINEINST>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
    // Count table reads and hits
    u64 num_branches = 0;
    for (u64 offset=0; offset<LINEINST; offset++) {
        if (static_cast<bool>(is_branch[offset])) num_branches++;
    }

    for (u64 j=0; j<NUMG; j++) {
        perf_table_reads[j] += num_branches;

        // Count hits for this table (any tag match, not just longest)
        for (u64 offset=0; offset<LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;
            val<NUMG> all_hits = val<NUMG+1>{match[offset]} >> 1; // all tag matches, strip bimodal bit
            val<1> hit_j = (all_hits >> j) & val<1>{1};
            if (static_cast<bool>(hit_j)) {
                perf_table_hits[j]++;
            }
        }
    }
#endif
#endif

        arr<val<1>,LINEINST> branch_taken = [&](u64 offset){
            return (actualdirs & update_mask[offset]) != hard<0>{};
        };

        //TODO: fix it
        branch_taken.fanout(hard<8>{});

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

        // select some candidate entries for allocation
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
#ifdef HASH_TAG
        arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == (last_offset^htag[i]);};
#else
        arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == concat(last_offset,htag[i]);};
#endif
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
    // Track allocations per table
    for (u64 j=0; j<NUMG; j++) {
        if (static_cast<bool>(allocate[j])) {
            perf_table_alloc[j]++;
        }
    }
#endif
#endif
#ifdef HASH_TAG
        arr<val<1>,NUMG> bdir = [&](u64 i) {
            arr<val<1>,LINEINST> tag_oh = [&](u64 offset){
                return match[offset].make_array(val<1>{})[i];
            };
            arr<val<LOGLINEINST>,LINEINST> offset_arr = [&](u64 offset){
                return select(tag_oh[offset],val<LOGLINEINST>{offset},val<LOGLINEINST>{0});
            };
            val<LOGLINEINST> tag_offset = offset_arr.fo1().fold_or();
            val<LOGLINEINST> offset = select(allocate[i],last_offset,tag_offset.fo1());
            offset.fanout(hard<LINEINST>{});
            arr<val<1>,LINEINST> match_offset = [&](u64 j){return branch_offset[j] == offset;};
            return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
        };
#else
        // associate a branch direction to each global table
        arr<val<1>,NUMG> bdir = [&](u64 i) {
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            val<LOGLINEINST> offset = select(allocate[i],last_offset,tag_offset.fo1());
            offset.fanout(hard<LINEINST>{});
            arr<val<1>,LINEINST> match_offset = [&](u64 j){return branch_offset[j] == offset;};
            return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
        };
#endif
        bdir.fanout(hard<2>{});

        // tell if global prediction is incorrect
        arr<val<1>,NUMG> badpred1 = [&](u64 i){
            return readc[i] != bdir[i];
        };
        badpred1.fanout(hard<3>{});
#ifdef HASH_TAG
        // associate to each global table a bit telling if local prediction differs from secondary prediction
        arr<val<1>,NUMG> altdiffer = [&](u64 i){
            arr<val<1>,LINEINST> tag_oh = [&](u64 offset){
                return match2[offset].make_array(val<1>{})[i];
            };
            arr<val<LOGLINEINST>,LINEINST> offset_arr = [&](u64 offset){
                return select(tag_oh[offset],val<LOGLINEINST>{offset},val<LOGLINEINST>{0});
            };
            val<LOGLINEINST> tag_offset = offset_arr.fo1().fold_or();
            return readc[i] != pred2.select(tag_offset.fo1());
        };

        // associate to each global table a bit telling if prediction for owning branch is correct
        arr<val<1>,NUMG> goodpred = [&](u64 i){
            arr<val<1>,LINEINST> tag_oh = [&](u64 offset){
                return match[offset].make_array(val<1>{})[i];
            };
            arr<val<LOGLINEINST>,LINEINST> offset_arr = [&](u64 offset){
                return select(tag_oh[offset],val<LOGLINEINST>{offset},val<LOGLINEINST>{0});
            };
            val<LOGLINEINST> tag_offset = offset_arr.fo1().fold_or();
            return (tag_offset.fo1() != last_offset) | correct_pred;
        };
#else
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
#endif
        // do P1 and P2 agree?
        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
        disagree_mask.fanout(hard<2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});

        // read the P1 hysteresis if P1 and P2 disagree
        arr<val<1>,LINEINST> p1_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff disagreement and hysteresis is weak
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(index1); // hyst=0 means weak
            });
        };

        // read the bimodal hysteresis if bimodal caused a misprediction
        arr<val<1>,LINEINST> b_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff cause of misprediction and hysteresis is weak
            val<1> bim_primary = actual_match1[offset] >> NUMG;
            return execute_if(bim_primary.fo1() & primary_wrong[offset], [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
        };

        // determine which primary global predictions are incorrect with a weak hysteresis
        arr<val<1>,NUMG> g_weak = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & badpred1[i] & (readh[i]==hard<0>{});
        };

        // need extra cycle for modifying prediction bits and for TAGE allocation
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
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
        sc_wrong_arr.fanout(hard<2>{});
        
        
        
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

        thre_guard_arr.fanout(hard<2>{});
        val<1> sc_need_update = arr<val<1>,LINEINST>{[&](u64 offset){
            return is_branch[offset] & do_update_arr[offset] & prov_hit_arr[offset];
        }}.fold_or();
        val<1> extra_cycle = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{}) | sc_need_update;
#else
        val<1> extra_cycle = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
#endif
        extra_cycle.fanout(hard<NUMG*2+1>{});
        need_extra_cycle(extra_cycle);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_extra_cycle_total      += static_cast<u64>(extra_cycle);
        perf_extra_cycle_badpred    += static_cast<u64>(some_badpred1);
        perf_extra_cycle_mispredict += static_cast<u64>(mispredict);
        perf_extra_cycle_p1_update  += static_cast<u64>(disagree_mask != hard<0>{});
#ifdef MY_SC
        perf_extra_cycle_sc_update  += static_cast<u64>(sc_need_update);
#endif
#endif
#endif

#ifdef USE_META
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
#ifdef HASH_TAG
        for (u64 i=0; i<NUMG; i++) {
            execute_if(allocate[i], [&](){gtag[i].write(gindex[i],last_offset^htag[i]);});
        }
#else
        // overwrite the tag in the allocated entry (mispredict)
        for (u64 i=0; i<NUMG; i++) {
            execute_if(allocate[i], [&](){gtag[i].write(gindex[i],concat(last_offset,htag[i]));});
        }

#endif
        // update the u bits
        arr<val<1>,NUMG> update_u = [&](u64 i){
            return primary[i] & altdiffer[i].fo1();
        };
        // if all post entries have the u bit set, reset their u bits
        val<1> noalloc = (candallocmask == hard<0>{});
        val<NUMG> uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<2>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
    if (static_cast<bool>(noalloc) && static_cast<bool>(mispredict)) {
        perf_alloc_failures++;
        if (static_cast<bool>(postmask == hard<0>{})) {
            perf_alloc_fail_highest++;
        } else {
            perf_alloc_fail_noubit++;
        }
    }
#endif
#endif
        for (u64 i=0; i<NUMG; i++) {
            execute_if(update_u[i].fo1() | allocate[i] | uclear[i], [&]() {
                val<1> newu = goodpred[i].fo1() & ~allocate[i] & ~uclear[i];
                ubit[i].write(gindex[i],newu.fo1(),extra_cycle);
            });
        }

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
                bim[offset].write(bindex,branch_taken[offset],extra_cycle);
            });
        }
        // update bimodal hysteresis if bimodal is primary provider
        for (u64 offset=0; offset<LINEINST; offset++) {
            val<1> bim_primary = match1[offset] >> NUMG;
            execute_if(is_branch[offset] & bim_primary.fo1(), [&](){
                bhyst[offset].write(bindex,~primary_wrong[offset]);
            });
        }

        // update incorrect global prediction if primary provider and the hysteresis is weak;
        // initialize global prediction in the allocated entry
        for (u64 i=0; i<NUMG; i++) {
            execute_if(g_weak[i].fo1() | allocate[i], [&](){
                gpred[i].write(gindex[i],bdir[i]);
            });
        }
        // update global prediction hysteresis if primary provider or allocated entry
        for (u64 i=0; i<NUMG; i++) {
            execute_if(primary[i] | allocate[i], [&](){
                // if allocated entry, set hysteresis to 0;
                // otherwise, increment hysteresis if correct pred, decrement if incorrect
                val<2> newhyst = select(allocate[i],val<2>{0},update_ctr(readh[i],~badpred1[i]));
                ghyst[i].write(gindex[i],newhyst.fo1(),extra_cycle);
            });
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

        // update global history
        val<1> line_end = block_entry >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            next_pc.fanout(hard<2>{});
            global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
            gfolds.update(val<PATHBITS>{next_pc>>2});
        });


#ifdef MY_SC


        //global_thre: majority vote across all qualifying offsets
        arr<val<1>,LINEINST> thre_update_en = [&](u64 offset){
            return is_branch[offset] & do_update_arr[offset] & thre_guard_arr[offset] & prov_hit_arr[offset];
        };
        thre_update_en.fanout(hard<2>{});
        // arr<val<1>,LINEINST> thre_wrong_bits = [&](u64 offset){
        //     return thre_update_en[offset] & sc_wrong_arr[offset];
        // };
        
        // thre_wrong_bits.fanout(hard<1>{});

        
        bias_high_idx.fanout(hard<LINEINST>{});
        //per-offset updates: bias_pc and thre1 (no write conflicts)
        for (u64 offset = 0; offset < LINEINST; offset++) {
            val<LOGBIAS-LOGLINEINST-2> high_idx = bias_high_idx;
            val<PERCWIDTH,i64> old_bias = bias_map[offset];
            val<2> write_tage_info = tage_info[offset];
            val<PRE_PC_THREBITS> old_thre1 = thre1[offset];
            val<PRE_PC_THREBITS> new_thre1 = update_ctr(old_thre1, sc_wrong_arr[offset]);
            old_bias.fanout(hard<4>{});
            write_tage_info.fanout(hard<4>{});
            high_idx.fanout(hard<4>{});
            for (u64 bank = 0; bank < 4; bank++){
                execute_if((write_tage_info==val<2>{bank}) & prov_hit_arr[offset] & (is_branch[offset] & do_update_arr[offset]), [&](){
                    bias_pc[bank][offset].write(high_idx, update_ctr(old_bias, branch_taken[offset]));
                });
            }
            
            thre1[offset] = select(thre_update_en[offset],new_thre1,thre1[offset]);
        }
#endif

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
    // Track prediction sources and accuracy - only for offsets that have branches
    for (u64 offset=0; offset<LINEINST; offset++) {
        if (!static_cast<bool>(is_branch[offset])) continue;

        perf_predictions++;

        // Use stored source from predict2() — avoids stale meta issue
        u64 src = static_cast<u64>(val<2>{pred_source_stored[offset]});
        // bimodal = bit NUMG, tables = bits 0..NUMG-1
        val<NUMG> prov_mask = val<NUMG+1>{pred_match1_stored[offset]} & hard<(1<<NUMG)-1>{};
        val<NUMG> alt_mask  = val<NUMG+1>{pred_match2_stored[offset]} & hard<(1<<NUMG)-1>{};

        if (src == 2) { // alt
            for (u64 j=0; j<NUMG; j++) {
                if (static_cast<u64>(alt_mask >> j) & 1) { perf_alt_used[j]++; break; }
            }
        } else if (src == 1) { // provider
            for (u64 j=0; j<NUMG; j++) {
                if (static_cast<u64>(prov_mask >> j) & 1) { perf_provider_used[j]++; break; }
            }
        } else { // bimodal
            perf_bimodal_used++;
        }

        // Check if prediction was correct
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
        }

        // Look up PC for this offset by scanning branch slots
        u64 pc_val = 0;
        for (u64 bi = 0; bi < num_branch; bi++) {
            if (static_cast<u64>(branch_offset[bi]) == offset) {
                pc_val = static_cast<u64>(branch_pc[bi]);
                break;
            }
        }

        // match1 is NUMG+1 bits: bits 0..NUMG-1 = TAGE tables, bit NUMG = bimodal
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

        // Determine prediction source from stored value
        u64 pred_source = src; // already computed above from pred_source_stored
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

        // Find alloc event for this branch (only on mispredict, last branch)
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

        // Record exec trace entry via SQLite
        {
            u64 conf_val = 0;
            if (pred_source == 1 || pred_source == 2) {
                u64 tbl = (pred_source == 1) ? pred_table :
                          static_cast<u64>(static_cast<u64>(alt_mask) ? [&]{ for(u64 j=0;j<NUMG;j++) if((static_cast<u64>(alt_mask)>>j)&1) return j; return (u64)0; }() : 0);
                if (tbl < NUMG) conf_val = static_cast<u64>(readh[tbl]);
            }
            // Compute sc_override: p2 differs from tage_p2 for this offset
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

        // Confidence distribution: bucket hysteresis of provider table
        if (pred_source == 1 || pred_source == 2) {
            u64 tbl = (pred_source == 1) ? pred_table :
                      static_cast<u64>(static_cast<u64>(alt_mask) ? [&]{ for(u64 j=0;j<NUMG;j++) if((static_cast<u64>(alt_mask)>>j)&1) return j; return (u64)0; }() : 0);
            if (tbl < NUMG) {
                u64 hval = static_cast<u64>(readh[tbl]);
                if (hval < 4) perf_conf[tbl][hval]++;
            }
        }

        // Misprediction trace (per-PC summary)
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
    // SC source counters: only track when SC overrides TAGE
    for (u64 offset=0; offset<LINEINST; offset++) {
        if (!static_cast<bool>(is_branch[offset])) continue;
        val<1> p2_bit = (p2 >> offset) & val<1>{1};
        val<1> tage_p2_bit = (tage_p2 >> offset) & val<1>{1};
        val<1> actual = branch_taken[offset];
        val<1> sc_override = (p2_bit != tage_p2_bit);
        if (static_cast<bool>(sc_override)) {
            perf_sc_override++;
            perf_sc_override_correct += static_cast<u64>(p2_bit == actual);
        }
        // threshold update counters
        if (static_cast<bool>(is_branch[offset] & do_update_arr[offset] & thre_guard_arr[offset] & prov_hit_arr[offset])) {
            perf_thre_update++;
            if (static_cast<bool>(sc_wrong_arr[offset]))
                perf_thre_update_inc++;
            else
                perf_thre_update_dec++;
        }
    }
#endif
#endif
#endif



        num_branch = 0; // done
    }
};
