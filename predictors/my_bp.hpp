// this is a basic TAGE, not necessarily well optimized

#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <sqlite3.h>

// #define SC
#ifndef SC
#define USE_ALT_PRED
#endif
#define HASH_TAG
#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include "tutorial/energy_monitor.hpp"
using namespace hcm;
#ifdef DEBUG_ENERGY
    struct energy_monitor monitor;
#endif
// #define PERF_COUNTERS
template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=400, u64 LOGP1=14, u64 GHIST1=6
, u64 NUMBANKS=4,u64 NUMWAYS=1,u64 CTRBIT=3,u64 UBIT=2>
struct my_bp : predictor {
    static_assert(LOGLB>2);
    static_assert(NUMG>0);
    static constexpr u64 LOGBIAS = 8;
    static constexpr u64 PERCWIDTH = 6;
    static constexpr u64 MINHIST = 4;
    static constexpr u64 USE_ALT_PRED_BITS = 4;
    static constexpr u64 UCTRBITS = 4;
    static constexpr u64 PATHBITS = 11;

    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGB > LOGLINEINST);
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 bindex_bits = LOGB-LOGLINEINST;

    static_assert(TAGW > LOGLINEINST);
#ifdef HASH_TAG
    static constexpr u64 HTAGBITS = TAGW;
#else
    static constexpr u64 HTAGBITS = TAGW-LOGLINEINST;
#endif
    static constexpr u64 NUMGSETS = (1<<LOGG)/NUMBANKS/NUMWAYS;
    static constexpr u64 NUMBSETS = (1<<(bindex_bits))/NUMBANKS;
    static constexpr u64 NUMP1SETS = (1<<(index1_bits))/NUMBANKS;
    static constexpr u64 NUMGTSETS = (1<<LOGG)/NUMWAYS;

    //rwram will compute the bank and local index for us, but we need to compute the folded global history indexes and tags
    static constexpr u64 BANKBITS       = static_cast<u64>(std::log2(static_cast<double>(NUMBANKS)));
    static constexpr u64 GINDEXBITS     = static_cast<u64>(std::log2(static_cast<double>(NUMGSETS)));
    
    static constexpr u64 BINDEXBITS     = static_cast<u64>(std::log2(static_cast<double>(NUMBSETS)));
    static constexpr u64 P1INDEXBITS    = static_cast<u64>(std::log2(static_cast<double>(NUMP1SETS)));
    // static constexpr u64 BANK = static_cast<u64>(std::log2(static_cast<double>(NUMGSETS)));
    geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;
    reg<1> true_block = 1;

    reg<GHIST1> path_history;
    reg<index1_bits> index1;
    
    arr<reg<1>,LINEINST> readp1;
    reg<LINEINST> p1;

    reg<bindex_bits> bindex;
    arr<reg<LOGG>,NUMG> gindex;
    arr<reg<HTAGBITS>,NUMG> htag;

    arr<reg<1>,LINEINST> readb;
    arr<reg<1>,LINEINST> readb_low;  // Store bimodal hysteresis
    arr<reg<TAGW>,NUMG> readt[NUMWAYS];
    arr<reg<1>,NUMG> readctr_pred[NUMWAYS];
    arr<reg<CTRBIT-1>,NUMG> readctr_cnt[NUMWAYS];
    arr<reg<UBIT>,NUMG> readu[NUMWAYS];

    arr<reg<NUMG>,NUMWAYS> notumask;

    arr<reg<NUMG>,LINEINST> match[NUMWAYS];
    arr<reg<NUMG>,LINEINST> match_provider[NUMWAYS];
    arr<reg<NUMG>,LINEINST> match_alt[NUMWAYS];

    arr<reg<1>,LINEINST> provider[NUMWAYS];
    arr<reg<1>,LINEINST> alt[NUMWAYS];
    reg<LINEINST> p2;
    reg<LINEINST> tage_p2;

    // HCpred: high-confidence prediction (strongest non-weak entry)
    arr<reg<1>,LINEINST> hc_pred_val[NUMWAYS];      // HC predicted direction
    arr<reg<NUMG>,LINEINST> hc_pred_source[NUMWAYS]; // HC source bitmask (0=bimodal, one-hot=TAGE table)
    arr<reg<2>,LINEINST> hc_pred_type[NUMWAYS]; // HC source type: 0=bim, 1=prov, 2=alt


    arr<reg<1>,LINEINST> use_provider[NUMWAYS];
    reg<USE_ALT_PRED_BITS> use_alt_on_na;



    reg<UCTRBITS> uctr;


#ifdef PERF_COUNTERS
    // Performance counters (CHEATING_MODE only)
    u64 perf_p1_predictions = 0;
    u64 perf_p1_correct = 0;
    u64 perf_bim_predictions = 0;
    u64 perf_bim_correct = 0;
    u64 perf_tage_predictions[NUMWAYS][NUMG] = {};
    u64 perf_tage_correct[NUMWAYS][NUMG] = {};        // Correct when used as prediction source
    u64 perf_tage_alloc[NUMWAYS][NUMG] = {};
    u64 perf_tage_provider_used[NUMWAYS][NUMG] = {};  // Count when table is the provider (longest match)
    u64 perf_tage_provider_correct[NUMWAYS][NUMG] = {}; // Correct when used as provider
    u64 perf_tage_alt_used[NUMWAYS][NUMG] = {};       // Count when table is the alt prediction source
    u64 perf_tage_alt_correct[NUMWAYS][NUMG] = {};    // Correct when used as alt
    u64 perf_tage_ctr_updates[NUMWAYS][NUMG] = {};    // Count CTR updates
    u64 perf_tage_useful_updates[NUMWAYS][NUMG] = {}; // Count useful bit updates
    u64 perf_tage_reads[NUMWAYS][NUMG] = {};          // Count table reads
    u64 perf_tage_hits[NUMWAYS][NUMG] = {};           // Count table hits (matches)
    u64 perf_tage_conflicts[NUMWAYS][NUMG] = {};      // Count potential conflicts
    u64 perf_tage_conf[NUMWAYS][NUMG][4] = {};        // CTR confidence dist at hit: 0=low,1=med-lo,2=med-hi,3=high
    u64 perf_tage_ubit_resets_total = 0;
    u64 perf_tage_skip_alloc_total=0;
    u64 perf_tage_alloc_fail_total = 0;  // Count allocation failures (no free entries)
    u64 perf_extra_cycle_tage_update = 0;    // Extra cycles due to TAGE update
    u64 perf_extra_cycle_mispredict = 0;     // Extra cycles due to misprediction
    u64 perf_extra_cycle_bim_update = 0;     // Extra cycles due to bimodal update
    u64 perf_extra_cycle_p1_update = 0;      // Extra cycles due to P1 update
    u64 perf_extra_cycle_total = 0;          // Total extra cycles
    u64 perf_hc_prov_count[NUMWAYS] = {};  // times HCpred from provider per way
    u64 perf_hc_alt_count[NUMWAYS] = {};   // times HCpred from alt per way
    u64 perf_hc_bim_count[NUMWAYS] = {};   // times HCpred from bimodal per way
    u64 perf_hc_prov_correct[NUMWAYS] = {};
    u64 perf_hc_alt_correct[NUMWAYS] = {};
    u64 perf_hc_bim_correct[NUMWAYS] = {};

    u64 perf_useful_alloc = 0;
    u64 perf_useful_inc= 0;
    u64 perf_useful_reset= 0;

    // Update condition counters for pred, hyst, uctr
    u64 perf_gpred_update_prov_weak_wrong[NUMWAYS][NUMG] = {};
    u64 perf_gpred_update_hc_weak_wrong[NUMWAYS][NUMG] = {};
    u64 perf_gpred_update_alt_weak_wrong[NUMWAYS][NUMG] = {};
    u64 perf_gpred_update_alloc[NUMWAYS][NUMG] = {};
    u64 perf_ghyst_update_prov_hit[NUMWAYS][NUMG] = {};
    u64 perf_ghyst_update_hc_used[NUMWAYS][NUMG] = {};
    u64 perf_ghyst_update_alt_hit[NUMWAYS][NUMG] = {};
    u64 perf_ghyst_update_alloc[NUMWAYS][NUMG] = {};
    u64 perf_uctr_update_inc[NUMWAYS][NUMG] = {};
    u64 perf_uctr_update_alloc[NUMWAYS][NUMG] = {};
    u64 perf_uctr_update_reset[NUMWAYS][NUMG] = {};
    u64 perf_bim_pred_update_weak_wrong = 0;
    u64 perf_bim_hyst_update_total = 0;
    u64 perf_p1_pred_update_disagree_weak = 0;
    u64 perf_p1_hyst_update_total = 0;

    // SQLite streaming trace (no in-memory vectors)
    sqlite3      *trace_db   = nullptr;
    sqlite3_stmt *exec_stmt  = nullptr;  // INSERT into exec_trace
    sqlite3_stmt *bim_stmt   = nullptr;  // INSERT into bim_update_trace
    u64           trace_seq  = 0;        // global sequence counter

    void db_exec(const char *sql) {
        char *errmsg = nullptr;
        if (sqlite3_exec(trace_db, sql, nullptr, nullptr, &errmsg) != SQLITE_OK) {
            std::cerr << "SQLite error: " << errmsg << "\n";
            sqlite3_free(errmsg);
        }
    }

    void open_trace_db() {
        if (sqlite3_open("trace.db", &trace_db) != SQLITE_OK) {
            std::cerr << "Cannot open trace.db: " << sqlite3_errmsg(trace_db) << "\n";
            return;
        }
        db_exec("PRAGMA journal_mode=WAL;");
        db_exec("PRAGMA synchronous=NORMAL;");
        db_exec(
            "CREATE TABLE IF NOT EXISTS exec_trace ("
            " seq INTEGER PRIMARY KEY,"
            " cycle INTEGER, pc INTEGER, offset INTEGER,"
            " actual_dir INTEGER, predicted_dir INTEGER, mispredict INTEGER,"
            " pred_source TEXT, pred_table INTEGER, bim_index INTEGER,"
            " hit INTEGER, hit_table INTEGER, hit_gtag INTEGER, hit_gindex INTEGER,"
            " alloc INTEGER, alloc_table INTEGER, alloc_gindex INTEGER, alloc_tag INTEGER"
            ");"
        );
        db_exec(
            "CREATE TABLE IF NOT EXISTS bim_update_trace ("
            " seq INTEGER PRIMARY KEY,"
            " cycle INTEGER, pc INTEGER, offset INTEGER, bim_index INTEGER,"
            " old_pred INTEGER, old_hyst INTEGER, actual_dir INTEGER,"
            " pred_written INTEGER, new_pred INTEGER,"
            " hyst_written INTEGER, new_hyst INTEGER,"
            " no_hit INTEGER, cond_extra_bim INTEGER,"
            " hc_is_bim INTEGER, prov_weak_wrong INTEGER,"
            " prov_pred INTEGER, prov_ctr INTEGER"
            ");"
        );
        db_exec("BEGIN;");

        sqlite3_prepare_v2(trace_db,
            "INSERT INTO exec_trace VALUES"
            "(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18);",
            -1, &exec_stmt, nullptr);
        sqlite3_prepare_v2(trace_db,
            "INSERT INTO bim_update_trace VALUES"
            "(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18);",
            -1, &bim_stmt, nullptr);
    }

    void close_trace_db() {
        if (!trace_db) return;
        db_exec("COMMIT;");
        // Create useful indexes for queries
        db_exec("CREATE INDEX IF NOT EXISTS idx_exec_pc         ON exec_trace(pc);");
        db_exec("CREATE INDEX IF NOT EXISTS idx_exec_mispredict ON exec_trace(mispredict);");
        db_exec("CREATE INDEX IF NOT EXISTS idx_exec_alloc      ON exec_trace(alloc);");
        db_exec("CREATE INDEX IF NOT EXISTS idx_bim_pc          ON bim_update_trace(pc);");
        sqlite3_finalize(exec_stmt);  exec_stmt = nullptr;
        sqlite3_finalize(bim_stmt);   bim_stmt  = nullptr;
        sqlite3_close(trace_db);      trace_db  = nullptr;
        std::cerr << "Trace written to trace.db\n";
    }

    // Insert one exec_trace row directly (no buffering)
    void insert_exec(u64 cycle, u64 pc, u64 offset,
                     u64 actual_dir, u64 predicted_dir, u64 mispredict,
                     const char *pred_source, u64 pred_table, u64 bim_index,
                     u64 hit, u64 hit_table, u64 hit_gtag, u64 hit_gindex,
                     u64 alloc, u64 alloc_table, u64 alloc_gindex, u64 alloc_tag) {
        sqlite3_reset(exec_stmt);
        sqlite3_bind_int64(exec_stmt,  1, static_cast<sqlite3_int64>(trace_seq++));
        sqlite3_bind_int64(exec_stmt,  2, static_cast<sqlite3_int64>(cycle));
        sqlite3_bind_int64(exec_stmt,  3, static_cast<sqlite3_int64>(pc));
        sqlite3_bind_int64(exec_stmt,  4, static_cast<sqlite3_int64>(offset));
        sqlite3_bind_int64(exec_stmt,  5, static_cast<sqlite3_int64>(actual_dir));
        sqlite3_bind_int64(exec_stmt,  6, static_cast<sqlite3_int64>(predicted_dir));
        sqlite3_bind_int64(exec_stmt,  7, static_cast<sqlite3_int64>(mispredict));
        sqlite3_bind_text (exec_stmt,  8, pred_source, -1, SQLITE_STATIC);
        sqlite3_bind_int64(exec_stmt,  9, static_cast<sqlite3_int64>(pred_table));
        sqlite3_bind_int64(exec_stmt, 10, static_cast<sqlite3_int64>(bim_index));
        sqlite3_bind_int64(exec_stmt, 11, static_cast<sqlite3_int64>(hit));
        sqlite3_bind_int64(exec_stmt, 12, static_cast<sqlite3_int64>(hit_table));
        sqlite3_bind_int64(exec_stmt, 13, static_cast<sqlite3_int64>(hit_gtag));
        sqlite3_bind_int64(exec_stmt, 14, static_cast<sqlite3_int64>(hit_gindex));
        sqlite3_bind_int64(exec_stmt, 15, static_cast<sqlite3_int64>(alloc));
        sqlite3_bind_int64(exec_stmt, 16, static_cast<sqlite3_int64>(alloc_table));
        sqlite3_bind_int64(exec_stmt, 17, static_cast<sqlite3_int64>(alloc_gindex));
        sqlite3_bind_int64(exec_stmt, 18, static_cast<sqlite3_int64>(alloc_tag));
        sqlite3_step(exec_stmt);
    }

    // Insert one bim_update_trace row directly (no buffering)
    void insert_bim(u64 cycle, u64 pc, u64 offset, u64 bim_index,
                    u64 old_pred, u64 old_hyst, u64 actual_dir,
                    u64 pred_written, u64 new_pred,
                    u64 hyst_written, u64 new_hyst,
                    u64 no_hit, u64 cond_extra_bim,
                    u64 hc_is_bim, u64 prov_weak_wrong,
                    u64 prov_pred, u64 prov_ctr) {
        sqlite3_reset(bim_stmt);
        sqlite3_bind_int64(bim_stmt,  1, static_cast<sqlite3_int64>(trace_seq++));
        sqlite3_bind_int64(bim_stmt,  2, static_cast<sqlite3_int64>(cycle));
        sqlite3_bind_int64(bim_stmt,  3, static_cast<sqlite3_int64>(pc));
        sqlite3_bind_int64(bim_stmt,  4, static_cast<sqlite3_int64>(offset));
        sqlite3_bind_int64(bim_stmt,  5, static_cast<sqlite3_int64>(bim_index));
        sqlite3_bind_int64(bim_stmt,  6, static_cast<sqlite3_int64>(old_pred));
        sqlite3_bind_int64(bim_stmt,  7, static_cast<sqlite3_int64>(old_hyst));
        sqlite3_bind_int64(bim_stmt,  8, static_cast<sqlite3_int64>(actual_dir));
        sqlite3_bind_int64(bim_stmt,  9, static_cast<sqlite3_int64>(pred_written));
        sqlite3_bind_int64(bim_stmt, 10, static_cast<sqlite3_int64>(new_pred));
        sqlite3_bind_int64(bim_stmt, 11, static_cast<sqlite3_int64>(hyst_written));
        sqlite3_bind_int64(bim_stmt, 12, static_cast<sqlite3_int64>(new_hyst));
        sqlite3_bind_int64(bim_stmt, 13, static_cast<sqlite3_int64>(no_hit));
        sqlite3_bind_int64(bim_stmt, 14, static_cast<sqlite3_int64>(cond_extra_bim));
        sqlite3_bind_int64(bim_stmt, 15, static_cast<sqlite3_int64>(hc_is_bim));
        sqlite3_bind_int64(bim_stmt, 16, static_cast<sqlite3_int64>(prov_weak_wrong));
        sqlite3_bind_int64(bim_stmt, 17, static_cast<sqlite3_int64>(prov_pred));
        sqlite3_bind_int64(bim_stmt, 18, static_cast<sqlite3_int64>(prov_ctr));
        sqlite3_step(bim_stmt);
    }

    void print_perf_counters() {
        std::cerr << "\n╔════════════════════════════════════════════════════════════════╗\n";
        std::cerr << "║              BRANCH PREDICTOR PERFORMANCE COUNTERS              ║\n";
        std::cerr << "╚════════════════════════════════════════════════════════════════╝\n";

        // P1 statistics
        std::cerr << "\n┌─ Level 1 Predictor ─────────────────────────────────────────────┐\n";
        std::cerr << "│ Predictions: " << std::setw(50) << std::left << perf_p1_predictions << "│\n";
        std::cerr << "│ Correct:     " << std::setw(50) << std::left << perf_p1_correct << "│\n";
        if (perf_p1_predictions > 0) {
            double accuracy = (100.0 * perf_p1_correct) / perf_p1_predictions;
            std::cerr << "│ Accuracy:    " << std::setw(50) << std::left << std::fixed << std::setprecision(2) << accuracy << "%│\n";
        }
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Bimodal statistics
        std::cerr << "\n┌─ Bimodal Predictor ─────────────────────────────────────────────┐\n";
        std::cerr << "│ Predictions: " << std::setw(50) << std::left << perf_bim_predictions << "│\n";
        std::cerr << "│ Correct:     " << std::setw(50) << std::left << perf_bim_correct << "│\n";
        if (perf_bim_predictions > 0) {
            double accuracy = (100.0 * perf_bim_correct) / perf_bim_predictions;
            std::cerr << "│ Accuracy:    " << std::setw(50) << std::left << std::fixed << std::setprecision(2) << accuracy << "%│\n";
        }
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Aggregate HC statistics for merging with TAGE
        u64 total_prov_count = 0, total_alt_count = 0, total_bim_count = 0;
        u64 total_prov_correct = 0, total_alt_correct = 0, total_bim_correct = 0;
        for (u64 w=0; w<NUMWAYS; w++) {
            total_prov_count += perf_hc_prov_count[w];
            total_alt_count += perf_hc_alt_count[w];
            total_bim_count += perf_hc_bim_count[w];
            total_prov_correct += perf_hc_prov_correct[w];
            total_alt_correct += perf_hc_alt_correct[w];
            total_bim_correct += perf_hc_bim_correct[w];
        }

        // TAGE statistics per way per table with HCpred source info
        for (u64 w=0; w<NUMWAYS; w++) {
            std::cerr << "\n┌─ TAGE Way " << w << " ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
            std::cerr << "│ Tbl │ HistLen │ Reads │ Hits │ Hit% │ Prov │ PrvAcc% │ Alt │ AltAcc% │ Total │ Use% │ TotAcc% │ Alloc │ Conf0(lo) │ Conf1(mlo) │ Conf2(mhi) │ Conf3(hi) │\n";
            std::cerr << "├─────┼─────────┼───────┼──────┼──────┼──────┼─────────┼─────┼─────────┼───────┼──────┼─────────┼───────┼───────────┼────────────┼────────────┼───────────┤\n";

            for (u64 j=0; j<NUMG; j++) {
                u64 hits      = perf_tage_hits[w][j];
                u64 reads     = perf_tage_reads[w][j];
                u64 prov      = perf_tage_provider_used[w][j];
                u64 prov_ok   = perf_tage_provider_correct[w][j];
                u64 alt_u     = perf_tage_alt_used[w][j];
                u64 alt_ok    = perf_tage_alt_correct[w][j];
                u64 total     = prov + alt_u;
                u64 total_ok  = prov_ok + alt_ok;
                u64 alloc     = perf_tage_alloc[w][j];
                u64 c0 = perf_tage_conf[w][j][0], c1 = perf_tage_conf[w][j][1];
                u64 c2 = perf_tage_conf[w][j][2], c3 = perf_tage_conf[w][j][3];
                u64 conf_total = c0+c1+c2+c3;

                std::cerr << "│ " << std::setw(3) << std::left  << j             << " │ ";
                std::cerr <<        std::setw(7) << std::right << gfolds.HLEN[j] << " │ ";
                std::cerr <<        std::setw(5) << std::right << reads           << " │ ";
                std::cerr <<        std::setw(4) << std::right << hits            << " │ ";

                // Hit%
                if (reads > 0)
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(4) << std::right << (100.0*hits/reads) << "% │ ";
                else
                    std::cerr << " N/A │ ";

                // Provider count + accuracy
                std::cerr << std::setw(4) << std::right << prov << " │ ";
                if (prov > 0)
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(6) << std::right << (100.0*prov_ok/prov) << "% │ ";
                else
                    std::cerr << "   N/A │ ";

                // Alt count + accuracy
                std::cerr << std::setw(3) << std::right << alt_u << " │ ";
                if (alt_u > 0)
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(6) << std::right << (100.0*alt_ok/alt_u) << "% │ ";
                else
                    std::cerr << "   N/A │ ";

                // Total used + use%
                std::cerr << std::setw(5) << std::right << total << " │ ";
                if (hits > 0)
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(4) << std::right << (100.0*total/hits) << "% │ ";
                else
                    std::cerr << " N/A │ ";

                // Total accuracy
                if (total > 0)
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(6) << std::right << (100.0*total_ok/total) << "% │ ";
                else
                    std::cerr << "   N/A │ ";

                // Alloc
                std::cerr << std::setw(5) << std::right << alloc << " │ ";

                // Confidence distribution at hit time: CTR value 0/1/2/3
                if (conf_total > 0) {
                    std::cerr << std::setw(5) << std::right << c0 << "(" << std::fixed << std::setprecision(0) << std::setw(3) << std::right << (100.0*c0/conf_total) << "%) │ ";
                    std::cerr << std::setw(5) << std::right << c1 << "(" << std::fixed << std::setprecision(0) << std::setw(3) << std::right << (100.0*c1/conf_total) << "%)  │ ";
                    std::cerr << std::setw(5) << std::right << c2 << "(" << std::fixed << std::setprecision(0) << std::setw(3) << std::right << (100.0*c2/conf_total) << "%)  │ ";
                    std::cerr << std::setw(5) << std::right << c3 << "(" << std::fixed << std::setprecision(0) << std::setw(3) << std::right << (100.0*c3/conf_total) << "%) │\n";
                } else {
                    std::cerr << "       N/A │         N/A │         N/A │       N/A │\n";
                }
            }
            std::cerr << "└─────┴─────────┴───────┴──────┴──────┴──────┴─────────┴─────┴─────────┴───────┴──────┴─────────┴───────┴───────────┴────────────┴────────────┴───────────┘\n";
        }

        std::cerr << "\n┌─ Skip Allocation Statistics ─────────────────────────────────┐\n";
        std::cerr << "│ Total Skip Alloc: " << std::setw(42) << std::left << perf_tage_skip_alloc_total << "│\n";
        std::cerr << "└───────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ UBIT Reset Statistics ───────────────────────────────────────┐\n";
        std::cerr << "│ Total UBIT Resets: " << std::setw(41) << std::left << perf_tage_ubit_resets_total << "│\n";
        std::cerr << "└───────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ Allocation Failure Statistics ───────────────────────────────┐\n";
        std::cerr << "│ Total Alloc Failures: " << std::setw(38) << std::left
                  << perf_tage_alloc_fail_total << "│\n";
        std::cerr << "└───────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ Extra Cycle Statistics ──────────────────────────────────────┐\n";
        std::cerr << "│ Total Extra Cycles:     " << std::setw(36) << std::left << perf_extra_cycle_total << "│\n";
        std::cerr << "│   TAGE Update:          " << std::setw(36) << std::left << perf_extra_cycle_tage_update << "│\n";
        std::cerr << "│   Misprediction:        " << std::setw(36) << std::left << perf_extra_cycle_mispredict << "│\n";
        std::cerr << "│   Bimodal Update:       " << std::setw(36) << std::left << perf_extra_cycle_bim_update << "│\n";
        std::cerr << "│   P1 Update:            " << std::setw(36) << std::left << perf_extra_cycle_p1_update << "│\n";
        std::cerr << "└───────────────────────────────────────────────────────────────┘\n";
        printf("Useful Allocations: %lu\n", static_cast<unsigned long>(perf_useful_alloc));
        printf("Useful Increments: %lu\n", static_cast<unsigned long>(perf_useful_inc));
        printf("Useful Resets: %lu\n", static_cast<unsigned long>(perf_useful_reset));
        // HCpred Source Statistics - merged with TAGE info
        std::cerr << "\n┌─ Prediction Source Distribution (via HCpred) ─────────────────┐\n";
        std::cerr << "│ Source   │ Count │ Correct │ Accuracy │ % of Total │\n";
        std::cerr << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";

        u64 total_tage_count = total_prov_count + total_alt_count;
        u64 total_tage_correct = total_prov_correct + total_alt_correct;
        u64 total_all = total_tage_count + total_bim_count;

        // TAGE Provider row
        std::cerr << "│ TAGE Prov│ " << std::setw(5) << std::right << total_prov_count << " │ ";
        std::cerr << std::setw(7) << std::right << total_prov_correct << " │ ";
        if (total_prov_count > 0) {
            double acc = (100.0 * total_prov_correct) / total_prov_count;
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right << acc << "% │ ";
        } else {
            std::cerr << "    N/A │ ";
        }
        if (total_all > 0) {
            double pct = (100.0 * total_prov_count) / total_all;
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right << pct << "% │\n";
        } else {
            std::cerr << "    N/A │\n";
        }

        // TAGE Alt row
        std::cerr << "│ TAGE Alt │ " << std::setw(5) << std::right << total_alt_count << " │ ";
        std::cerr << std::setw(7) << std::right << total_alt_correct << " │ ";
        if (total_alt_count > 0) {
            double acc = (100.0 * total_alt_correct) / total_alt_count;
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right << acc << "% │ ";
        } else {
            std::cerr << "    N/A │ ";
        }
        if (total_all > 0) {
            double pct = (100.0 * total_alt_count) / total_all;
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right << pct << "% │\n";
        } else {
            std::cerr << "    N/A │\n";
        }

        // Bimodal row
        std::cerr << "│ Bimodal  │ " << std::setw(5) << std::right << total_bim_count << " │ ";
        std::cerr << std::setw(7) << std::right << total_bim_correct << " │ ";
        if (total_bim_count > 0) {
            double acc = (100.0 * total_bim_correct) / total_bim_count;
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right << acc << "% │ ";
        } else {
            std::cerr << "    N/A │ ";
        }
        if (total_all > 0) {
            double pct = (100.0 * total_bim_count) / total_all;
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right << pct << "% │\n";
        } else {
            std::cerr << "    N/A │\n";
        }

        // Total row
        std::cerr << "├──────────┼───────┼─────────┼──────────┼────────────┤\n";
        std::cerr << "│ Total    │ " << std::setw(5) << std::right << total_all << " │ ";
        u64 total_correct = total_tage_correct + total_bim_correct;
        std::cerr << std::setw(7) << std::right << total_correct << " │ ";
        if (total_all > 0) {
            double acc = (100.0 * total_correct) / total_all;
            std::cerr << std::fixed << std::setprecision(2) << std::setw(8) << std::right << acc << "% │ ";
            std::cerr << std::fixed << std::setprecision(1) << std::setw(9) << std::right << 100.0 << "% │\n";
        } else {
            std::cerr << "    N/A │     N/A │\n";
        }
        std::cerr << "└──────────┴───────┴─────────┴──────────┴────────────┘\n";

        // Update Condition Statistics
        std::cerr << "\n╔════════════════════════════════════════════════════════════════╗\n";
        std::cerr << "║           UPDATE CONDITION BREAKDOWN STATISTICS                 ║\n";
        std::cerr << "╚════════════════════════════════════════════════════════════════╝\n";

        for (u64 w=0; w<NUMWAYS; w++) {
            std::cerr << "\n┌─ TAGE Way " << w << " Update Conditions ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐\n";
            std::cerr << "│ Tbl │ HistLen │ Pred Updates                                    │ Hyst Updates                                    │ Uctr Updates                    │\n";
            std::cerr << "│     │         │ ProvWW │ HcWW │ AltWW │ Alloc │ Total │ HystUpd │ ProvHit │ HcUse │ AltHit │ Alloc │ Total │ Inc │ Alloc │ Reset │ Total │\n";
            std::cerr << "├─────┼─────────┼────────┼──────┼───────┼───────┼───────┼─────────┼─────────┼───────┼────────┼───────┼───────┼─────┼───────┼───────┼───────┤\n";

            for (u64 j=0; j<NUMG; j++) {
                u64 pred_prov_ww = perf_gpred_update_prov_weak_wrong[w][j];
                u64 pred_hc_ww = perf_gpred_update_hc_weak_wrong[w][j];
                u64 pred_alt_ww = perf_gpred_update_alt_weak_wrong[w][j];
                u64 pred_alloc = perf_gpred_update_alloc[w][j];
                u64 pred_total = pred_prov_ww + pred_hc_ww + pred_alt_ww + pred_alloc;

                u64 hyst_prov = perf_ghyst_update_prov_hit[w][j];
                u64 hyst_hc = perf_ghyst_update_hc_used[w][j];
                u64 hyst_alt = perf_ghyst_update_alt_hit[w][j];
                u64 hyst_alloc = perf_ghyst_update_alloc[w][j];
                u64 hyst_total = hyst_prov + hyst_hc + hyst_alt + hyst_alloc;

                u64 uctr_inc = perf_uctr_update_inc[w][j];
                u64 uctr_alloc = perf_uctr_update_alloc[w][j];
                u64 uctr_reset = perf_uctr_update_reset[w][j];
                u64 uctr_total = uctr_inc + uctr_alloc + uctr_reset;

                std::cerr << "│ " << std::setw(3) << std::left << j << " │ ";
                std::cerr << std::setw(7) << std::right << gfolds.HLEN[j] << " │ ";
                std::cerr << std::setw(6) << std::right << pred_prov_ww << " │ ";
                std::cerr << std::setw(4) << std::right << pred_hc_ww << " │ ";
                std::cerr << std::setw(5) << std::right << pred_alt_ww << " │ ";
                std::cerr << std::setw(5) << std::right << pred_alloc << " │ ";
                std::cerr << std::setw(5) << std::right << pred_total << " │ ";
                std::cerr << std::setw(7) << std::right << hyst_total << " │ ";
                std::cerr << std::setw(7) << std::right << hyst_prov << " │ ";
                std::cerr << std::setw(5) << std::right << hyst_hc << " │ ";
                std::cerr << std::setw(6) << std::right << hyst_alt << " │ ";
                std::cerr << std::setw(5) << std::right << hyst_alloc << " │ ";
                std::cerr << std::setw(5) << std::right << hyst_total << " │ ";
                std::cerr << std::setw(3) << std::right << uctr_inc << " │ ";
                std::cerr << std::setw(5) << std::right << uctr_alloc << " │ ";
                std::cerr << std::setw(5) << std::right << uctr_reset << " │ ";
                std::cerr << std::setw(5) << std::right << uctr_total << " │\n";
            }
            std::cerr << "└─────┴─────────┴────────┴──────┴───────┴───────┴───────┴─────────┴─────────┴───────┴────────┴───────┴───────┴─────┴───────┴───────┴───────┘\n";
        }

        std::cerr << "\n┌─ Bimodal Update Conditions ─────────────────────────────────────┐\n";
        std::cerr << "│ Pred Updates (weak & wrong): " << std::setw(31) << std::left << perf_bim_pred_update_weak_wrong << "│\n";
        std::cerr << "│ Hyst Updates (total):        " << std::setw(31) << std::left << perf_bim_hyst_update_total << "│\n";
        std::cerr << "└──────────────────────────────────────────────────────────────────┘\n";

        std::cerr << "\n┌─ P1 Update Conditions ──────────────────────────────────────────┐\n";
        std::cerr << "│ Pred Updates (disagree & weak): " << std::setw(28) << std::left << perf_p1_pred_update_disagree_weak << "│\n";
        std::cerr << "│ Hyst Updates (total):           " << std::setw(28) << std::left << perf_p1_hyst_update_total << "│\n";
        std::cerr << "└──────────────────────────────────────────────────────────────────┘\n";

        // Commit and close the SQLite trace database
        close_trace_db();

        std::cerr << "\n";
    }

#endif



    u64 num_branch = 0;
    u64 block_size = 0;



    arr<reg<LOGLINEINST>,LINEINST> branch_offset;
    arr<reg<64>,LINEINST> branch_pc;
    arr<reg<1>,LINEINST> branch_dir;
    reg<LINEINST> inst_oh;

    rwram<1,(1<<index1_bits),NUMBANKS> table1_pred[LINEINST] {"P1 pred"};

    //205ps
    rwram<TAGW,NUMGTSETS,NUMBANKS> gtag[NUMWAYS][NUMG] {"tags"};
    // ram<val<TAGW>,(1<<LOGG)> gtag[NUMWAYS][NUMG] {"tags"}; // tags
    //183ps
    rwram<1,NUMGTSETS,NUMBANKS> gpred[NUMWAYS][NUMG] {"gpred"};
    //
    rwram<UBIT,NUMGTSETS,NUMBANKS> ubit[NUMWAYS][NUMG] {"uctr"};

    //116ps
    rwram<1,1<<bindex_bits,NUMBANKS> bim_hi[LINEINST] {"bpred"};


    #ifdef SC
    //bias_pc is concat pc and tage
    //idx = (pc>>(LOGLB+2)) 2bit tage info(prov_weak,prov_taken)
    rwram<PERCWIDTH,(1<<LOGBIAS)/NUMBANKS,NUMBANKS> bias_pc {"Bias pc"};
    arr<reg<PERCWIDTH>,4> bias_lmap {"Bias LMAP"};

    rwram<PERCWIDTH,(1<<LOGBIAS)/NUMBANKS,NUMBANKS> bias_pc {"Bias pc"};
    #endif

    zone UPDATE_ONLY;
    rwram<1,1<<index1_bits,NUMBANKS> table1_hyst[LINEINST] {"P1 hyst"};
    rwram<1,1<<bindex_bits,NUMBANKS> bim_low[LINEINST] {"bhyst"};
    rwram<CTRBIT-1,NUMGTSETS,NUMBANKS> ghyst[NUMWAYS][NUMG] {"ghyst"};

    my_bp()
    {
        constexpr u64 lineinst = 1 << LOGLINEINST;
        constexpr u64 total_bits = (1ULL << index1_bits) * lineinst
            + (1ULL << LOGG) * (TAGW + CTRBIT + UBIT) * NUMG * NUMWAYS
            + (1ULL << bindex_bits) * 2 * lineinst;
#ifdef VERBOSE
        std::cerr << "TAGE history lengths: ";
        for (u64 i=0; i<NUMG; i++) std::cerr << gfolds.HLEN[i] << " ";
        std::cerr << std::endl;
        std::cerr << "Total storage: " << total_bits << " bits (" << (total_bits / 8192.0) << " KB)" << std::endl;
#endif
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        open_trace_db();
#endif
#endif
    }


    void new_block(val<64> inst_pc)
    {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        inst_oh = offset.fo1().decode().concat();
        inst_oh.fanout(hard<6*LINEINST+1>{});
        block_size = 1;
    }

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        new_block(inst_pc);
        // val<std::max(P1INDEXBITS+BANKBITS,GHIST1)> lineaddr = inst_pc >> (LOGLB);
        // val<BANKBITS> bankid = inst_pc >> (LOGLB);
        val<index1_bits> raw_index = inst_pc >> (LOGLB);
        raw_index.fanout(hard<2>{});
        if constexpr (GHIST1 <= index1_bits) {
            index1 = raw_index ^ (val<index1_bits>{path_history}<<(index1_bits-GHIST1));
        } else {
            index1 = path_history.make_array(val<index1_bits>{}).append(raw_index).fold_xor();
        }
        index1.fanout(hard<LINEINST+1>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            readp1[offset] = table1_pred[offset].read(index1);
        }
        p1 = readp1.fo1().concat();
        p1.fanout(hard<LINEINST+1>{});
        return (inst_oh & p1) != hard<0>{};
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        return ((inst_oh<<block_size) & p1) != hard<0>{};
    };

    val<1> is_weak(val<1> pred, val<CTRBIT-1> ctr){
        return select(pred, ctr==ctr.minval, ctr==ctr.maxval);
    }
    void tage_predict(val<64> inst_pc)
    {
        val<HTAGBITS>   raw_tag = inst_pc >> (LOGLB+LOGG);
        val<LOGG>       raw_gindex = inst_pc >> LOGLB;

        val<bindex_bits> raw_bindex = inst_pc >> (LOGLB);

        raw_tag.fanout(hard<NUMG>{});
        raw_gindex.fanout(hard<NUMG>{});
        gfolds.fanout(hard<2>{});
        

        bindex = raw_bindex.fo1();
        bindex.fanout(hard<2*LINEINST>{});

        // gfolds 140ps
        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = raw_gindex ^ gfolds.template get<0>(i);
        }
        gindex.fanout(hard<3*NUMWAYS>{});

        for (u64 i=0; i<NUMG; i++) {
            htag[i] = raw_tag ^ gfolds.template get<1>(i);
        }
        // htag.fanout(hard<2>{});

        //150ps
        for (u64 offset=0; offset<LINEINST; offset++) {
            readb[offset] = bim_hi[offset].read(bindex);
            readb_low[offset] = bim_low[offset].read(bindex);
        }
        // readb.fanout(hard<NUMWAYS>{});

        //260ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for(u64 j=0; j<NUMG; j++){
                //146ps
                readt[w][j] = gtag[w][j].read(gindex[j]);
                readu[w][j] = ubit[w][j].read(gindex[j]);
                readctr_pred[w][j] = gpred[w][j].read(gindex[j]);
                readctr_cnt[w][j] = ghyst[w][j].read(gindex[j]);
            }
        }

        for (u64 w=0; w<NUMWAYS; w++) {
            readt[w].fanout(hard<LINEINST>{});
            readu[w].fanout(hard<2>{});
            readctr_cnt[w].fanout(hard<2>{});
            match_provider[w].fanout(hard<5>{});
            readctr_pred[w].fanout(hard<2>{});
        }

        //20ps
        arr<val<NUMG>,NUMWAYS> umask = [&](u64 w){
            arr<val<1>,NUMG> u = [&](u64 j){ return readu[w][j]!=hard<0>{}; };
            return u.concat();
        };

        for(u64 way=0;way<NUMWAYS;way++){
            notumask[way] = ~umask[way].fo1();
        }
        notumask.fanout(hard<LINEINST*NUMWAYS+1>{});

        //70ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
#ifdef HASH_TAG
                arr<val<1>,NUMG> wayhit = [&](u64 j){
                    return readt[w][j] == (htag[j] ^ val<HTAGBITS>{offset});
                };
#else
                arr<val<1>,NUMG> wayhit = [&](u64 j){
                    return (val<HTAGBITS>{readt[w][j]} == htag[j]) & (val<LOGLINEINST>{readt[w][j]>>HTAGBITS} == val<LOGLINEINST>{offset});
                };                
#endif
                match[w][offset] = wayhit.concat();
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) match[w].fanout(hard<3>{});

        //50ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                match_provider[w][offset] = match[w][offset].one_hot();
            }
        }

        arr<val<NUMG>,NUMWAYS> gpreds = [&](u64 w){ return readctr_pred[w].concat(); };
        gpreds.fanout(hard<3*LINEINST>{});
        //30ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                provider[w][offset] = (match_provider[w][offset] & gpreds[w]) != hard<0>{};
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) provider[w].fanout(hard<3>{});

        //20ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                match_alt[w][offset] = (match[w][offset] ^ match_provider[w][offset]).one_hot();
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) match_alt[w].fanout(hard<2>{});
        //20ps
        //bim pred is in alt
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                alt[w][offset] = (match_alt[w][offset] & gpreds[w]) != hard<0>{};
            }
        }
        // for (u64 w=0; w<NUMWAYS; w++) alt[w].fanout(hard<2>{});


        // use_alt_on_na.fanout(hard<2>{});
        arr<val<1>,USE_ALT_PRED_BITS> use_alt_on_na_array = use_alt_on_na.make_array(val<1>{});
        val<1> use_alt_on_na_pos = use_alt_on_na_array[USE_ALT_PRED_BITS-1];
        arr<val<1>,LINEINST> use_alt_on_na_pos_dup = use_alt_on_na_pos.fo1().replicate(hard<LINEINST>{});//reduce fanout
        use_alt_on_na_pos_dup.fanout(hard<NUMWAYS*2>{});
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 inst=0; inst<LINEINST; inst++) {
                arr<val<1>,NUMG> prov_oh = val<NUMG>{match_provider[w][inst]}.make_array(val<1>{});
                arr<val<1>,NUMG> alt_oh = val<NUMG>{match_alt[w][inst]}.make_array(val<1>{});
                prov_oh.fanout(hard<3>{});
                arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                    return select(prov_oh[j],readctr_cnt[w][j],val<CTRBIT-1>{0});
                };
                val<CTRBIT-1> prov_ctr  = prov_ctrs.fo1().fold_or();
                // provider & alt
                val<1> prov_is_weak = (prov_ctr == val<CTRBIT-1>{0}) & 
                (match_provider[w][inst] & (notumask.fold_or()!=val<NUMG>{0}));
                prov_is_weak.fanout(hard<2>{});
#ifndef SC
                use_provider[w][inst] = prov_oh.fold_or() &
                select(alt_oh.fo1().fold_or(),(~prov_is_weak | ~use_alt_on_na_pos_dup[inst]),val<1>{1});
#else
                use_provider[w][inst] = prov_oh.fold_or();
#endif
                // ── HCpred ──────────────────────────────────────────────
                // Provider: not weak and exists → HCpred = provider
                
                val<1> use_provider_hc = prov_oh.fold_or() & (~prov_is_weak | ~use_alt_on_na_pos_dup[inst]);
                use_provider_hc.fanout(hard<3>{});
                // Alt: from all alt hits (not just one_hot), find longest non-weak
                val<NUMG> all_alt_hits = val<NUMG>{match[w][inst]} ^ val<NUMG>{match_provider[w][inst]};
                arr<val<1>,NUMG> hc_alt_cand = [&](u64 j){
                    val<1> hit = all_alt_hits.make_array(val<1>{})[j];
                    val<1> not_weak = ~(readctr_cnt[w][j] == val<CTRBIT-1>{0});
                    return hit.fo1() & not_weak.fo1();
                };
                hc_alt_cand.fanout(hard<2>{});
                // Select the longest (highest-indexed) non-weak alt
                val<NUMG> hc_alt_mask = hc_alt_cand.concat().one_hot();
                hc_alt_mask.fanout(hard<2>{});
                val<1> has_hc_alt = hc_alt_cand.fold_or();
                val<1> hc_alt_pred_val = (hc_alt_mask & gpreds[w]) != hard<0>{};

                val<1> use_hcpred = ~use_provider_hc & has_hc_alt.fo1();
                use_hcpred.fanout(hard<2>{});
                // Final HCpred value and source
                val<1> hc_val = select(use_provider_hc, provider[w][inst],
                                select(use_hcpred,     hc_alt_pred_val.fo1(),readb[inst]));
                val<NUMG> hc_src = select(use_provider_hc, val<NUMG>{match_provider[w][inst]},
                                   select(use_hcpred,      hc_alt_mask,
                                                            val<NUMG>{0}));

                hc_pred_val[w][inst]    = hc_val.fo1();
                hc_pred_source[w][inst] = hc_src.fo1();
#ifdef PERF_COUNTERS
                val<2> hc_type = select(use_provider_hc, val<2>{1},   // provider
                                 select(use_hcpred,       val<2>{2},   // alternate
                                                           val<2>{0})); // bimodal
                hc_pred_type[w][inst] = hc_type;
#endif
            }
        }

        // p2 now uses HCpred: strongest non-weak entry (provider > longest non-weak alt > bimodal)
#ifndef SC
        // timing is bad (If tage is 3 cycle,this may be used)
        tage_p2 = arr<val<1>,LINEINST>{[&](u64 offset){
            arr<val<1>,NUMWAYS> pred_use_tage = [&](u64 w){
                return hc_pred_source[w][offset] != val<NUMG>{0};
            };
            arr<val<1>,NUMWAYS> final_pred = [&](u64 w){
                return select(pred_use_tage[w], val<1>{hc_pred_val[w][offset]}, val<1>{0});
            };
            return select(pred_use_tage.fold_or(), final_pred.fo1().fold_or(), hc_pred_val[0][offset]);
        }}.concat();
#else
        tage_p2 = arr<val<1>,LINEINST>{[&](u64 offset){
            arr<val<1>,NUMWAYS> pred_use_tage = [&](u64 w){
                return match_provider[w][offset] != val<NUMG>{0};
            };
            arr<val<1>,NUMWAYS> final_pred = [&](u64 w){
                return select(pred_use_tage[w], val<1>{provider[w][offset]}, val<1>{0});
            };
            return select(pred_use_tage.fold_or(), final_pred.fo1().fold_or(), readb[offset]);
        }}.concat();

#endif

    }
    val<1> predict2(val<64> inst_pc)
    {

        tage_predict(inst_pc);
#ifndef SC
        p2 = tage_p2;
#endif
        p2.fanout(hard<LINEINST+2>{});
        val<1> taken = (inst_oh & p2) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{inst_oh>>(LINEINST-1)});

        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = ((inst_oh<<block_size) & p2) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{inst_oh>>(LINEINST-1-block_size)});
        block_size++;
        return taken;
    }

    void update_condbr(val<64> branch_pc, val<1> taken, [[maybe_unused]] val<64> next_pc)
    {
        //TODO: fix fanout
        assert(num_branch < LINEINST);
        val<LOGLINEINST> offset = branch_pc >> 2;
        branch_offset[num_branch] = offset;
        this->branch_pc[num_branch] = branch_pc;
        branch_dir[num_branch] = taken.fo1();
        num_branch++;
    }

    void update_cycle(instruction_info &block_end_info) override
    {

        val<1> &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc = block_end_info.next_pc;
        gfolds.fanout(hard<2>{});
        if (num_branch == 0) {
            val<1> line_end = inst_oh >> (LINEINST-block_size);
            val<1> actual_block = ~(true_block & line_end.fo1());
            actual_block.fanout(hard<GHIST+NUMG*2+2>{});
            execute_if(actual_block, [&](){
                next_pc.fanout(hard<2>{});
                path_history = (path_history << 1) ^ val<GHIST1>{next_pc>>2};
                gfolds.update(val<PATHBITS>{next_pc>>2});
                true_block = 1;
            });
            return;
        }

        mispredict.fanout(hard<4>{});

        val<1> correct_pred = ~mispredict;
        gindex.fanout(hard<3*NUMWAYS>{});
        for (u64 w=0; w<NUMWAYS; w++) {
            match_provider[w].fanout(hard<10>{});
            provider[w].fanout(hard<5>{});
            // alt[w].fanout(hard<2>{});
            readctr_cnt[w].fanout(hard<4>{});
        }
        branch_offset.fanout(hard<2*LINEINST+3*NUMWAYS>{});
        branch_dir.fanout(hard<LINEINST+NUMWAYS*(NUMG+1)+1>{});
        
        index1.fanout(hard<LINEINST*3>{});

        bindex.fanout(hard<LINEINST*3>{});
        readb.fanout(hard<NUMWAYS+1>{});


        val<LOGLINEINST> last_offset = branch_offset[num_branch-1];
        last_offset.fanout(hard<2*NUMG*NUMWAYS>{});

        u64 update_valid = (u64(1)<<num_branch)-1;
        arr<val<LINEINST>,LINEINST> update_mask = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){ return branch_offset[i] == offset; };
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});

        arr<val<1>,LINEINST> is_branch = [&](u64 offset){
            return update_mask[offset] != hard<0>{};
        };
        is_branch.fanout(hard<(5+4*NUMG)*NUMWAYS+2>{});

        val<LINEINST> branch_mask = is_branch.concat();

        // ==================== P1 Disagree Check (before extra_cycle) ====================
        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
        disagree_mask.fanout(hard<2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});

        // Read P1 hysteresis if P1 and P2 disagree (must read before extra_cycle)
        arr<val<1>,LINEINST> p1_weak = [&](u64 offset) -> val<1> {
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(index1);
            });
        };

        arr<val<1>,LINEINST> actualdirs = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){ return (branch_offset[i] == offset) & branch_dir[i]; };
            return (match_offset.fo1().concat() & update_valid) != val<LINEINST>{0};
        };
        actualdirs.fanout(hard<2+(4*NUMG+1)*NUMWAYS>{});

        // Aggregate match_provider across branches (like tage.hpp actual_match1)
        arr<val<NUMG>,NUMWAYS> primary_mask = [&](u64 w){
            arr<val<NUMG>,LINEINST> m = [&](u64 offset){
                return val<NUMG>{is_branch[offset].replicate(hard<NUMG>{}).concat()} & 
                val<NUMG>{match_provider[w][offset]};
            };
            return m.fo1().fold_or();
        };
        primary_mask.fanout(hard<4>{});

        /*
        bim update:
        1. pred wrong
        2. pred right and weak
        tage update:
        1. ctr
            1. pred right increase final pred table
            2. pred wrong decrease final pred table
            3. alloc set ctr to 1<<(CTRBIT-1)
        2. tag: only update when alloc
        3. ubit:
            1. provider right ,alt wrong increase provider
            2. alloc,set ubit to 0
            3. reset

        */
        // ==================== Extra Cycle Determination ====================
        // Determine which updates are needed for extra cycle
        // 1. TAGE update: any table was used as HCpred, provider, or extra alt
        // 2. Bimodal update: no hit or HCpred fell back to bimodal
        // 3. P1 update: P1 and P2 disagree
        // 4. Misprediction: need allocation

        // Check if any TAGE table needs update (HCpred used, provider hit, or extra alt)
        val<1> some_tage_update = (primary_mask.fold_or() != val<NUMG>{0});


        // Check if any bimodal needs update
        arr<val<1>,LINEINST> bim_update_arr = [&](u64 offset){
            val<1> branch_valid = is_branch[offset];
            val<1> actual_dir = actualdirs[offset];
            // Condition 1: no TAGE hit at all
            arr<val<1>,NUMWAYS> has_provider_arr = [&](u64 w){
                return match_provider[w][offset] != hard<0>{};
            };
            has_provider_arr.fanout(hard<2>{});
            val<1> no_hit = ~has_provider_arr.fold_or();
            // Condition 2: HCpred fell back to bimodal (hc_pred_type == 0)
            arr<val<1>,NUMWAYS> hc_is_bim_arr = [&](u64 w){
                return hc_pred_source[w][offset] == hard<0>{};
            };
            val<1> hc_is_bim = hc_is_bim_arr.fo1().fold_or();

            // Provider weak & wrong check (for condition 2)
            arr<val<1>,NUMWAYS> prov_weak_wrong_arr = [&](u64 w){
                // use match_provider[w][offset] (per-offset), not primary_mask[w] (aggregated)
                arr<val<1>,NUMG> prov_bits = val<NUMG>{match_provider[w][offset]}.make_array(val<1>{});
                arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                    return select(prov_bits[j], readctr_cnt[w][j], val<CTRBIT-1>{0});
                };
                val<CTRBIT-1> prov_ctr = prov_ctrs.fo1().fold_or();
                val<1> prov_pred = provider[w][offset];
                prov_pred.fanout(hard<2>{});
                val<1> prov_weak = prov_ctr.fo1() == val<CTRBIT-1>{0};
                val<1> prov_wrong = prov_pred != actual_dir.fo1();
                return has_provider_arr[w] & prov_weak.fo1() & prov_wrong.fo1();
            };
            val<1> prov_weak_wrong = prov_weak_wrong_arr.fo1().fold_or();

            val<1> cond_extra_bim = hc_is_bim.fo1() & prov_weak_wrong.fo1();

            // Combined update condition
            return branch_valid.fo1() & (no_hit.fo1() | cond_extra_bim.fo1());
        };

        bim_update_arr.fanout(hard<3>{});
        val<1> some_bim_update = bim_update_arr.fold_or();

        val<1> some_p1_update = (disagree_mask != hard<0>{});
        val<1> extra_cycle = some_tage_update | mispredict | some_bim_update | some_p1_update.fo1();
        extra_cycle.fanout(hard<4*LINEINST+3*NUMWAYS*NUMG>{});



        need_extra_cycle(extra_cycle);

        // ==================== Bimodal Update ====================
        // Bimodal update conditions (HCpred-based):
        // 1. No TAGE hit at all (no provider) → always update bimodal
        // 2. HCpred fell back to bimodal (provider weak & wrong) → update bimodal
        
        for (u64 offset = 0; offset < LINEINST; offset++) {
            val<1> bim_pred = readb[offset];
            val<1> curr_hyst = readb_low[offset];
            val<1> actual_dir = actualdirs[offset];
            val<1> dir_match = bim_pred.fo1()==actual_dir;
            actual_dir.fanout(hard<2>{});
            dir_match.fanout(hard<2>{});
            //TODO：add is branch
            val<1> pred_need_update = (curr_hyst.fo1() == val<1>{0}) & bim_update_arr[offset] & ~dir_match;
            execute_if(bim_update_arr[offset], [&](){
                bim_low[offset].write(bindex, dir_match, extra_cycle);
            });
            execute_if(pred_need_update.fo1(), [&](){
                bim_hi[offset].write(bindex, actual_dir, extra_cycle);
            });
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            if (static_cast<bool>(is_branch[offset])) {
                if (static_cast<bool>(bim_update_arr[offset])) perf_bim_hyst_update_total++;
                if (static_cast<bool>(pred_need_update)) perf_bim_pred_update_weak_wrong++;
            }
#endif
#endif

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            // Record bim update trace for this offset if update triggered
            if (static_cast<bool>(is_branch[offset]) && static_cast<bool>(bim_update_arr[offset])) {
                // find PC for this offset
                u64 pc_for_offset = 0;
                u64 dir_for_offset = 0;
                for (u64 bi = 0; bi < num_branch; bi++) {
                    if (static_cast<u64>(branch_offset[bi]) == offset) {
                        pc_for_offset  = static_cast<u64>(branch_pc[bi]);
                        dir_for_offset = static_cast<u64>(branch_dir[bi]);
                        break;
                    }
                }
                // recompute no_hit and cond_extra_bim for this offset
                u64 no_hit_val = 1;
                for (u64 w = 0; w < NUMWAYS; w++) {
                    if (static_cast<u64>(val<NUMG>{match_provider[w][offset]}) != 0) {
                        no_hit_val = 0; break;
                    }
                }
                // hc_is_bim: any way has hc_pred_type==0
                u64 hc_is_bim_val = 0;
                for (u64 w = 0; w < NUMWAYS; w++) {
                    if (static_cast<u64>(val<2>{hc_pred_type[w][offset]}) == 0) {
                        hc_is_bim_val = 1; break;
                    }
                }
                // prov_weak_wrong: use match_provider[w][offset] (per-offset)
                u64 prov_weak_wrong_val = 0;
                u64 prov_pred_val = 0;
                u64 prov_ctr_val = 0;
                for (u64 w = 0; w < NUMWAYS; w++) {
                    u64 pmask = static_cast<u64>(val<NUMG>{match_provider[w][offset]});
                    if (pmask != 0) {
                        prov_pred_val = static_cast<u64>(provider[w][offset]);
                        for (u64 j = 0; j < NUMG; j++) {
                            if (pmask & (u64(1)<<j)) {
                                prov_ctr_val = static_cast<u64>(readctr_cnt[w][j]);
                                break;
                            }
                        }
                        u64 prov_wrong = (prov_pred_val != dir_for_offset) ? 1 : 0;
                        // is_weak: ctr at weak boundary (0 for taken, max for not-taken)
                        u64 max_ctr = (u64(1) << (CTRBIT-1)) - 1;
                        u64 prov_weak = (prov_pred_val == 1) ? (prov_ctr_val == 0 ? 1 : 0)
                                                              : (prov_ctr_val == max_ctr ? 1 : 0);
                        if (prov_weak && prov_wrong) prov_weak_wrong_val = 1;
                        break;
                    }
                }
                // Compute bim update outcome before inserting
                u64 dir_match_val = (static_cast<u64>(bim_pred) == dir_for_offset) ? 1 : 0;
                u64 pred_need_val = (static_cast<u64>(curr_hyst) == 0) && !dir_match_val ? 1 : 0;
                // Stream directly into SQLite
                insert_bim(
                    static_cast<u64>(panel.cycle), pc_for_offset, offset,
                    static_cast<u64>(bindex),
                    static_cast<u64>(bim_pred), static_cast<u64>(curr_hyst),
                    dir_for_offset,
                    pred_need_val, pred_need_val ? dir_for_offset : static_cast<u64>(bim_pred),
                    1, dir_match_val,
                    no_hit_val, hc_is_bim_val & prov_weak_wrong_val,
                    hc_is_bim_val, prov_weak_wrong_val,
                    prov_pred_val, prov_ctr_val
                );
                
            }
#endif
#endif
        }

        // ==================== TAGE Allocation ====================
        // Compute allocation mask (only on misprediction)
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
        // mispmask.fanout(hard<2>{});

        arr<val<NUMG>,NUMWAYS> last_match = [&](u64 w){
            arr<val<NUMG>,LINEINST> last = [&](u64 offset){
                return select(offset==last_offset,match_provider[w][offset],val<NUMG>{0});
            };
            return last.fold_or();
        };
        
        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match.append(1).concat().one_hot()-1);
        // postmask.fanout(hard<2>{});

        // Check if provider was not used but correct (to skip allocation)
        arr<val<1>,NUMWAYS> skip_alloc_arr = [&](u64 w){
            arr<val<NUMG>,LINEINST> prov_for_last = [&](u64 offset){
                return select(branch_offset[num_branch-1] == offset, val<NUMG>{match_provider[w][offset]}, val<NUMG>{0});
            };
            val<NUMG> last_prov_mask = prov_for_last.fo1().fold_or();
            val<1> has = last_prov_mask.fo1() != hard<0>{};

            arr<val<1>,LINEINST> used_for_last = [&](u64 offset){
                return (branch_offset[num_branch-1] == offset) & use_provider[w][offset].fo1();
            };
            val<1> used = used_for_last.fo1().fold_or();

            arr<val<1>,LINEINST> pred_for_last = [&](u64 offset){
                return (branch_offset[num_branch-1] == offset) & provider[w][offset];
            };
            val<1> pred = pred_for_last.fo1().fold_or();

            val<1> actual_last = val<1>{branch_dir[num_branch-1]};
            val<1> correct = pred.fo1() == actual_last.fo1();

            return has.fo1() & ~used.fo1() & correct.fo1();
        };
        val<1> skip_alloc = skip_alloc_arr.fo1().fold_or();


        val<NUMG> candallocmask = postmask.fo1() & notumask.fold_or() & ~(skip_alloc.replicate(hard<NUMG>{}).concat());
        // candallocmask.fanout(hard<1>{});
        val<1> alloc_fail = mispredict & (candallocmask == hard<0>{});



        // ==================== TAGE UBIT Reset Counter ====================
        // Increment uctr and check if reset is needed
        val<UCTRBITS> uctr_threshold = val<UCTRBITS>{1 << (UCTRBITS - 1)};
        val<1> should_reset = (uctr == uctr_threshold);

        // Update uctr: increment normally, reset to 0 when threshold reached
        uctr = select(should_reset, val<UCTRBITS>{0}, select(alloc_fail,val<UCTRBITS>{uctr + 1}, uctr));
        should_reset.fanout(hard<2*NUMWAYS*NUMG+1>{});


        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask_low = collamask.one_hot();
        collamask_low.fanout(hard<2>{});
        val<NUMG> collamask_high = (collamask ^ collamask_low).one_hot();
        val<NUMG> collamask_final = select(val<2>{static_cast<u64>(std::rand())%2}==hard<0>{}, collamask_high.fo1(), collamask_low);
        arr<val<1>,NUMG> allocate = collamask_final.reverse().make_array(val<1>{});
        // allocate.fanout(hard<1>{});


        u64 way_sel = static_cast<u64>(std::rand()) % NUMWAYS;

        // ==================== TAGE Update ====================
        // CTR update, tag update, and ubit update
        val<1> last_dir = branch_dir[num_branch-1];

        arr<val<NUMG>,NUMWAYS> hc_used = [&](u64 w){
            arr<val<NUMG>,LINEINST> hc_source_mask = [&](u64 offset){
                return select(is_branch[offset], val<NUMG>{hc_pred_source[w][offset]}, val<NUMG>{0});
            };
            return hc_source_mask.fo1().fold_or();
        };

        arr<val<NUMG>,NUMWAYS> provider_hit = [&](u64 w){
            arr<val<NUMG>,LINEINST> provider_mask = [&](u64 offset){
                return select(is_branch[offset], val<NUMG>{match_provider[w][offset]}, val<NUMG>{0});
            };
            return provider_mask.fo1().fold_or();
        };

        arr<val<NUMG>,NUMWAYS> alt_hits = [&](u64 w){
            arr<val<NUMG>,LINEINST> alt_mask = [&](u64 offset){
                return select(is_branch[offset], match_alt[w][offset], val<NUMG>{0});
            };
            return alt_mask.fo1().fold_or();
        };

        arr<val<NUMG>,NUMWAYS> prov_correct = [&](u64 w){
            arr<val<1>,NUMG> prov_correct_tmp = [&](u64 num){
                arr<val<1>,LINEINST> prov_correct_arr = [&](u64 offset){
                    val<1> j_hit = match_provider[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & j_hit & (provider[w][offset] == actualdirs[offset]);
                };
                return prov_correct_arr.fo1().fold_or();
            };
            return prov_correct_tmp.fo1().concat();
        };

        arr<val<NUMG>,NUMWAYS> prov_weak = [&](u64 w){
            arr<val<1>,NUMG> prov_weak_tmp = [&](u64 num){
                arr<val<1>,LINEINST> prov_weak_arr = [&](u64 offset){
                    val<1> j_hit = match_provider[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & j_hit & (readctr_cnt[w][num] == val<CTRBIT-1>{0});
                };
                return prov_weak_arr.fo1().fold_or();
            };
            return prov_weak_tmp.fo1().concat();
        };

        arr<val<NUMG>,NUMWAYS> alt_wrong = [&](u64 w){
            arr<val<1>,NUMG> alt_wrong_tmp = [&](u64 num){
                arr<val<1>,LINEINST> alt_wrong_arr = [&](u64 offset){
                    val<1> j_hit = match_alt[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & select(j_hit, (alt[w][offset] != actualdirs[offset]),readb[offset] != actualdirs[offset]);
                };
                return alt_wrong_arr.fo1().fold_or();
            };
            return alt_wrong_tmp.fo1().concat();
        };

        arr<val<NUMG>,NUMWAYS> base_cond = [&](u64 w){

            arr<val<1>,NUMG> base_cond_tmp = [&](u64 num){
                return prov_weak[w].make_array(val<1>{})[num] & (~prov_correct[w].make_array(val<1>{})[num]);
            };
            return base_cond_tmp.fo1().concat();
        };

        // gpred update cond
        arr<val<NUMG>,NUMWAYS> prov_gpred_need_update = [&](u64 w){
            arr<val<1>,NUMG> prov_gpred_need_update_tmp = [&](u64 num){
                return base_cond[w].make_array(val<1>{})[num] & (provider_hit[w].make_array(val<1>{})[num]);
            };
            return prov_gpred_need_update_tmp.fo1().concat();
        };
#ifdef UPDATEALTONWEAKMISP
        arr<val<NUMG>,NUMWAYS> hc_gpred_need_update = [&](u64 w){
            arr<val<1>,NUMG> hc_gpred_need_update_tmp = [&](u64 num){
                arr<val<1>,LINEINST> hc_mispred_weak_arr = [&](u64 offset){
                    val<1> j_hit = hc_pred_source[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & j_hit & (hc_pred_val[w][offset] != actualdirs[offset]) & (readctr_cnt[w][num] == val<CTRBIT-1>{0});
                };
                return hc_mispred_weak_arr.fo1().fold_or() ;
            };
            return hc_gpred_need_update_tmp.fo1().concat()& base_cond[w];
        };

        arr<val<NUMG>,NUMWAYS> alt_gpred_need_update = [&](u64 w){
            arr<val<1>,NUMG> alt_gpred_need_update_tmp = [&](u64 num){
                arr<val<1>,LINEINST> alt_mispred_weak_arr = [&](u64 offset){
                    val<1> j_hit = match_alt[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & j_hit & (alt[w][offset] != actualdirs[offset]) & (readctr_cnt[w][num] == val<CTRBIT-1>{0});
                };
                return alt_mispred_weak_arr.fo1().fold_or();
            };
            return alt_gpred_need_update_tmp.fo1().concat() & base_cond[w];
        };
#endif

        val<CTRBIT-1> ghyst_max = ghyst_max.maxval;
        // ghyst update cond !(weak&mispred | correct&sat) & hit
        arr<val<NUMG>,NUMWAYS> prov_ghyst_need_update = [&](u64 w){
            return provider_hit[w] & (~prov_gpred_need_update[w]);
        };
#ifdef UPDATEALTONWEAKMISP
        arr<val<NUMG>,NUMWAYS> hc_ghyst_need_update = [&](u64 w){
            return hc_used[w] & (~hc_gpred_need_update[w]);
        };
        arr<val<NUMG>,NUMWAYS> alt_ghyst_need_update = [&](u64 w){
            return alt_hits[w] & (~alt_gpred_need_update[w]);
        };
#endif

        arr<val<NUMG>,NUMWAYS> do_alloc = [&](u64 w){
            return (w == way_sel) ? allocate.concat() : val<NUMG>{0};
        };
        arr<val<NUMG>,NUMWAYS> matched_dir = [&](u64 w){
            arr<val<1>,NUMG> matched_dir_tmp = [&](u64 num){
                arr<val<1>,LINEINST> offset_match = [&](u64 offset){
                    val<1> j_matches = match[w][offset].make_array(val<1>{})[num];
                    return is_branch[offset] & j_matches;
                };
                arr<val<1>,LINEINST> offset_dir = [&](u64 offset){
                    return select(offset_match[offset], actualdirs[offset], val<1>{0});
                };
                return offset_dir.fo1().fold_or();
            };
            return matched_dir_tmp.fo1().concat();
        };


        arr<val<NUMG>,NUMWAYS> inc_useful = [&](u64 w){
            return prov_correct[w] & alt_wrong[w] & provider_hit[w];
        };


        // val<UBIT> new_u[NUMWAYS][NUMG];
        // val<UBIT> final_u[NUMWAYS][NUMG];
        // for (u64 w=0; w<NUMWAYS; w++) {
        //     for (u64 j=0; j<NUMG; j++) {
        //         new_u[w][j] = select(readu[w][j] == hard<new_u.maxval>{}, readu[w][j], val<UBIT>{readu[w][j] + 1});
        //         final_u[w][j] = select(do_alloc[w][j], val<UBIT>{0}, new_u[w][j]);
        //     }
        // }

        // val<UBIT> ubit_value[NUMWAYS][NUMG];
        // for (u64 w=0; w<NUMWAYS; w++) {
        //     for (u64 j=0; j<NUMG; j++) {
        //         ubit_value[w][j] = select(should_reset, val<UBIT>{0}, final_u[w][j]);
        //     }
        // }


        // ==================== SRAM writes  ====================
        for (u64 w=0; w<NUMWAYS; w++) {
#ifdef UPDATEALTONWEAKMISP
            arr<val<1>,NUMG> gpred_update = (prov_gpred_need_update[w] | hc_gpred_need_update[w] | alt_gpred_need_update[w]).make_array(val<1>{});
            arr<val<1>,NUMG> ghyst_update = (prov_ghyst_need_update[w] | hc_ghyst_need_update[w] | alt_ghyst_need_update[w]).make_array(val<1>{});
#else
            arr<val<1>,NUMG> gpred_update = prov_gpred_need_update[w].make_array(val<1>{});
            arr<val<1>,NUMG> ghyst_update = prov_ghyst_need_update[w].make_array(val<1>{});
#endif
            arr<val<1>,NUMG> do_alloc_arr = do_alloc[w].make_array(val<1>{});
            arr<val<1>,NUMG> dir_arr = matched_dir[w].make_array(val<1>{});

            arr<val<1>,NUMG> inc_dir = inc_useful[w].make_array(val<1>{});
            for (u64 j=0; j<NUMG; j++) {
                val<1> final_gdir = select(do_alloc_arr[j],last_dir, dir_arr[j]);
                val<CTRBIT-1> correct_hyst = update_ctr(dir_arr[j],readctr_cnt[w][j]);
                val<CTRBIT-1> final_ghyst = select(do_alloc_arr[j],val<CTRBIT-1>{0}, correct_hyst);
                val<TAGW> final_gtag = htag[j] ^ val<HTAGBITS>{last_offset};
                val<UBIT> update_uctr = update_ctr(inc_dir[j],readu[w][j]);
                val<UBIT> final_uctr = select(should_reset | do_alloc_arr[j] , val<UBIT>{0}, update_uctr);
                execute_if(gpred_update[j] | do_alloc_arr[j], [&](){
                    gpred[w][j].write(gindex[j], final_gdir, extra_cycle);
                });
                execute_if(ghyst_update[j] | do_alloc_arr[j], [&](){
                    ghyst[w][j].write(gindex[j], final_ghyst, extra_cycle);
                });
                execute_if(do_alloc_arr[j], [&](){
                    gtag[w][j].write(gindex[j], final_gtag, extra_cycle);
                });
                execute_if(do_alloc_arr[j] | inc_dir[j] | should_reset, [&](){
                    ubit[w][j].write(gindex[j], final_uctr, extra_cycle);
                });
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
                if (static_cast<bool>(prov_gpred_need_update[w].make_array(val<1>{})[j])) perf_gpred_update_prov_weak_wrong[w][j]++;
#ifdef UPDATEALTONWEAKMISP
                if (static_cast<bool>(hc_gpred_need_update[w].make_array(val<1>{})[j])) perf_gpred_update_hc_weak_wrong[w][j]++;
                if (static_cast<bool>(alt_gpred_need_update[w].make_array(val<1>{})[j])) perf_gpred_update_alt_weak_wrong[w][j]++;
#endif
                if (static_cast<bool>(do_alloc_arr[j])) perf_gpred_update_alloc[w][j]++;
                if (static_cast<bool>(prov_ghyst_need_update[w].make_array(val<1>{})[j])) perf_ghyst_update_prov_hit[w][j]++;
#ifdef UPDATEALTONWEAKMISP
                if (static_cast<bool>(hc_ghyst_need_update[w].make_array(val<1>{})[j])) perf_ghyst_update_hc_used[w][j]++;
                if (static_cast<bool>(alt_ghyst_need_update[w].make_array(val<1>{})[j])) perf_ghyst_update_alt_hit[w][j]++;
#endif
                if (static_cast<bool>(do_alloc_arr[j])) perf_ghyst_update_alloc[w][j]++;
                if (static_cast<bool>(inc_dir[j])) perf_uctr_update_inc[w][j]++;
                if (static_cast<bool>(do_alloc_arr[j])) perf_uctr_update_alloc[w][j]++;
                if (static_cast<bool>(should_reset)) perf_uctr_update_reset[w][j]++;
#endif
#endif
            }
        }


#ifdef USE_ALT_PRED
        // USE_ALT_ON_NA update - aggregate conditions using arrays

        // // Helper lambda to compute common USE_ALT_PRED conditions
        auto compute_use_alt_conditions = [&](u64 w){
            struct Conditions {
                val<1> has_provider;
                val<1> provider_is_weak;
                val<1> provider_pred;
                val<1> hc_pred;
                val<1> actual_dir;
            };
            arr<val<1>,NUMG> primary_bits = primary_mask[w].make_array(val<1>{});
            arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                return select(primary_bits[j], readctr_cnt[w][j], val<CTRBIT-1>{0});
            };
            val<CTRBIT-1> prov_ctr = prov_ctrs.fold_or();

            arr<val<1>,LINEINST> prov_pred = [&](u64 offset){
                return is_branch[offset] & provider[w][offset];
            };
            val<1> provider_pred = prov_pred.fold_or();
            provider_pred.fanout(hard<2>{});
            arr<val<1>,LINEINST> hc_pred_arr = [&](u64 inst){
                return select(match_provider[w][inst]!=val<NUMG>{0},hc_pred_val[w][inst],val<1>{0});
            };
            val<1> hc_pred = hc_pred_arr.fold_or();
            hc_pred.fanout(hard<2>{});

            arr<val<1>,LINEINST> actual_dir_arr = [&](u64 inst){
                return select(match_provider[w][inst]!=val<NUMG>{0},actualdirs[inst],val<1>{0});
            };

            val<1> has_provider = primary_mask[w] != hard<0>{};
            val<1> provider_is_weak = is_weak(provider_pred, prov_ctr);
            // removed unused alt_or_base_match

            return Conditions{has_provider, provider_is_weak, provider_pred, hc_pred, actual_dir_arr.fold_or()};
        };

        arr<val<1>,NUMWAYS> inc_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & (cond.hc_pred != cond.provider_pred) & (cond.hc_pred == cond.actual_dir);
        };

        arr<val<1>,NUMWAYS> dec_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & (cond.hc_pred != cond.provider_pred) & (cond.hc_pred != cond.actual_dir);
        };

        val<1> any_inc_use_alt = inc_use_alt_arr.fo1().fold_or();
        val<1> any_dec_use_alt = dec_use_alt_arr.fo1().fold_or();

        use_alt_on_na.fanout(hard<3>{});
        val<USE_ALT_PRED_BITS,i64> new_use_alt = select(any_inc_use_alt,
            val<USE_ALT_PRED_BITS,i64>{use_alt_on_na + 1},
            select(any_dec_use_alt,
                val<USE_ALT_PRED_BITS,i64>{use_alt_on_na - 1},
                use_alt_on_na));
        using use_alt_t = valt<decltype(use_alt_on_na)>;
        use_alt_on_na = select(new_use_alt > use_alt_t::maxval, use_alt_t{use_alt_t::maxval},
                              select(new_use_alt < use_alt_t::minval, use_alt_t{use_alt_t::minval}, use_alt_t{new_use_alt}));
#endif

        // ==================== P1 Update (after extra_cycle) ====================
        // Update P1 prediction if P1 and P2 disagree and hysteresis is weak
        arr<val<1>,LINEINST> p2_split = p2.make_array(val<1>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(p1_weak[offset].fo1(), [&](){
                table1_pred[offset].write(index1, p2_split[offset],extra_cycle);
            });
        }

        // Update P1 hysteresis for all branches
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(is_branch[offset], [&](){
                table1_hyst[offset].write(index1, ~disagree[offset],extra_cycle);
            });
        }
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        for (u64 offset=0; offset<LINEINST; offset++) {
            if (static_cast<bool>(is_branch[offset])) {
                if (static_cast<bool>(p1_weak[offset])) perf_p1_pred_update_disagree_weak++;
                perf_p1_hyst_update_total++;
            }
        }
#endif
#endif

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
    // Count extra cycle types
    perf_extra_cycle_tage_update += static_cast<u64>(some_tage_update);
    perf_extra_cycle_mispredict += static_cast<u64>(mispredict);
    perf_extra_cycle_bim_update += static_cast<u64>(some_bim_update);
    perf_extra_cycle_p1_update += static_cast<u64>(some_p1_update);
    perf_extra_cycle_total += static_cast<u64>(extra_cycle);

    // Build execution flow trace: one record per branch in this block
    {
        u64 is_misp = static_cast<u64>(mispredict);

        // Find alloc event: which table/index/tag was allocated (way_sel only)
        u64 alloc_found = 0, alloc_table = 0, alloc_gindex = 0, alloc_tag = 0;
        for (u64 j = 0; j < NUMG; j++) {
            val<1> do_alloc = (way_sel < NUMWAYS) ? allocate[j] : val<1>{0};
            if (static_cast<bool>(do_alloc)) {
                alloc_found  = 1;
                alloc_table  = j;
                alloc_gindex = static_cast<u64>(gindex[j]);
#ifdef HASH_TAG
                alloc_tag = static_cast<u64>(htag[j] ^ val<HTAGBITS>{last_offset});
#else
                alloc_tag = static_cast<u64>(concat(val<LOGLINEINST>{last_offset}, htag[j]));
#endif
                break;
            }
        }

        for (u64 i = 0; i < num_branch; i++) {
            u64 pc_val  = static_cast<u64>(branch_pc[i]);
            u64 dir_val = static_cast<u64>(branch_dir[i]);
            u64 idx     = static_cast<u64>(branch_offset[i]);
            u64 pred_val = static_cast<u64>((p2 >> idx) & val<LINEINST>{1});

            // Find provider hit for this branch offset
            u64 hit_found = 0, hit_table = NUMG, hit_gtag = 0, hit_gindex = 0;
            for (u64 w = 0; w < NUMWAYS && !hit_found; w++) {
                u64 pmask = static_cast<u64>(val<NUMG>{match_provider[w][idx]});
                if (pmask != 0) {
                    for (u64 j = 0; j < NUMG; j++) {
                        if ((pmask >> j) & 1) {
                            hit_found  = 1;
                            hit_table  = j;
                            hit_gtag   = static_cast<u64>(readt[w][j]);
                            hit_gindex = static_cast<u64>(gindex[j]);
                            break;
                        }
                    }
                }
            }

            // Determine prediction source from hc_pred_type/hc_pred_source (way 0)
            u64 pred_source = 0; // 0=bimodal
            u64 pred_table  = NUMG;
            {
                u64 hc_type = static_cast<u64>(val<2>{hc_pred_type[0][idx]});
                if (hc_type == 1 || hc_type == 2) {
                    pred_source = hc_type; // 1=provider, 2=alt
                    u64 src_mask = static_cast<u64>(val<NUMG>{hc_pred_source[0][idx]});
                    for (u64 j = 0; j < NUMG; j++) {
                        if ((src_mask >> j) & 1) { pred_table = j; break; }
                    }
                }
            }

            // alloc only applies to the last branch (the mispredicted one)
            u64 is_last = (i == num_branch - 1);
            // Stream directly into SQLite (no in-memory buffer)
            const char *src_str = (pred_source == 0) ? "bimodal"
                                : (pred_source == 1) ? "tage_prov" : "tage_alt";
            insert_exec(
                static_cast<u64>(panel.cycle), pc_val, idx,
                dir_val, pred_val,
                static_cast<u64>(is_misp & is_last),
                src_str, pred_table,
                static_cast<u64>(bindex),
                hit_found, hit_table, hit_gtag, hit_gindex,
                static_cast<u64>(is_misp & is_last & alloc_found),
                alloc_table, alloc_gindex, alloc_tag
            );
        }
    }
        perf_tage_skip_alloc_total += static_cast<u64>(skip_alloc);
        // Increment ubit reset counter
        perf_tage_ubit_resets_total += static_cast<u64>(should_reset);
        perf_tage_alloc_fail_total += static_cast<u64>(alloc_fail);
        for (u64 w = 0; w < NUMWAYS; w++) {
            for (u64 j = 0; j < NUMG; j++) {
                // Count per-offset: provider used / alt used / correct, and confidence distribution at hit
                for (u64 offset = 0; offset < LINEINST; offset++) {
                    val<1> is_br = is_branch[offset];
                    if (!static_cast<bool>(is_br)) continue;

                    val<1> is_hc_src_j = val<NUMG>{hc_pred_source[w][offset]}.make_array(val<1>{})[j];
                    val<2> hc_type = val<2>{hc_pred_type[w][offset]};
                    val<1> j_as_prov_hc = is_hc_src_j & (hc_type == val<2>{1});
                    val<1> j_as_alt_hc  = is_hc_src_j & (hc_type == val<2>{2});

                    val<1> actual_dir = actualdirs[offset];
                    val<1> pred_dir   = readctr_pred[w][j];

                    perf_tage_provider_used[w][j]    += static_cast<u64>(j_as_prov_hc);
                    perf_tage_provider_correct[w][j] += static_cast<u64>(j_as_prov_hc & (pred_dir == actual_dir));
                    perf_tage_alt_used[w][j]         += static_cast<u64>(j_as_alt_hc);
                    perf_tage_alt_correct[w][j]      += static_cast<u64>(j_as_alt_hc  & (pred_dir == actual_dir));

                    val<1> any_hit_j = val<NUMG>{match[w][offset]}.make_array(val<1>{})[j];
                    if (static_cast<bool>(any_hit_j & is_br)) {
                        u64 ctr_val = static_cast<u64>(concat(readctr_pred[w][j], readctr_cnt[w][j]));
                        u64 conf_idx = ctr_val >> (CTRBIT - 2);
                        if (conf_idx < 4) perf_tage_conf[w][j][conf_idx]++;
                    }
                }
            }
        }
        // Count P1 predictions and correct predictions
        perf_p1_predictions += num_branch;

        arr<val<1>,LINEINST> p1_split = p1.make_array(val<1>{});
        arr<val<1>,LINEINST> p1_correct_arr = [&](u64 offset){
            return is_branch[offset] & (p1_split[offset] == actualdirs[offset]);
        };
        val<LINEINST> p1_correct_mask = p1_correct_arr.concat();
        auto p1_correct_count = p1_correct_mask.ones();
        perf_p1_correct += p1_correct_count;

        // Count bimodal predictions and correct predictions
        arr<val<1>,LINEINST> bim_used_arr = [&](u64 offset){
            arr<val<1>,NUMWAYS> no_provider = [&](u64 w){
                return match_provider[w][offset] == hard<0>{};
            };
            return is_branch[offset] & no_provider.fo1().fold_and();
        };
        val<LINEINST> bim_used_mask = bim_used_arr.concat();
        auto bim_used_count = bim_used_mask.ones();
        perf_bim_predictions += bim_used_count;

        arr<val<1>,LINEINST> bim_correct_arr = [&](u64 offset){
            return bim_used_arr[offset] & (readb[offset] == actualdirs[offset]);
        };
        val<LINEINST> bim_correct_mask = bim_correct_arr.concat();
        auto bim_correct_count = bim_correct_mask.ones();
        perf_bim_correct += bim_correct_count;

        // Count TAGE predictions, correct predictions, and allocations
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 j=0; j<NUMG; j++) {
                // Count allocations
                val<1> do_alloc = (w == way_sel) ? allocate[j] : val<1>{0};
                perf_tage_alloc[w][j] += static_cast<u64>(do_alloc);

                // Count reads: all tables are read for each branch
                perf_tage_reads[w][j] += num_branch;

                // Count hits: when table j matched for way w
                arr<val<1>,NUMG> prim_bits = primary_mask[w].make_array(val<1>{});
                val<1> was_hit = prim_bits[j];
                perf_tage_hits[w][j] += static_cast<u64>(was_hit);

                // Count potential conflicts: when table matches but prediction is wrong
                val<1> pred_val = readctr_pred[w][j];
                arr<val<1>,LINEINST> conflict_arr = [&](u64 offset){
                    arr<val<1>,NUMG> match_bits = primary_mask[w].make_array(val<1>{});
                    return is_branch[offset] & match_bits[j] & (pred_val != actualdirs[offset]);
                };
                val<1> was_conflict = conflict_arr.fo1().fold_or();
                perf_tage_conflicts[w][j] += static_cast<u64>(was_conflict);
            }
        }

        // Count HCpred source type and accuracy - per way for TAGE, once for bimodal
        // For each branch and each way, count which source was used
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                val<1> is_br = is_branch[offset];
                val<2> hc_t = hc_pred_type[w][offset];

                // Check which source was used in this way
                val<1> prov_hit = (val<NUMG>{match_provider[w][offset]} != hard<0>{});
                val<1> hc_src_hit = (val<NUMG>{hc_pred_source[w][offset]} != hard<0>{});

                val<1> is_prov = (hc_t == val<2>{1}) & prov_hit;
                val<1> is_alt = (hc_t == val<2>{2}) & hc_src_hit;

                // Check if prediction was correct
                val<1> hc_correct = val<1>{hc_pred_val[w][offset]} == actualdirs[offset];

                // Count provider and alt per way
                perf_hc_prov_count[w]   += static_cast<u64>(is_br & is_prov);
                perf_hc_alt_count[w]    += static_cast<u64>(is_br & is_alt);

                perf_hc_prov_correct[w] += static_cast<u64>(is_br & is_prov & hc_correct);
                perf_hc_alt_correct[w]  += static_cast<u64>(is_br & is_alt & hc_correct);
            }
        }

        // Count bimodal once per branch (not per way, since it's global)
        for (u64 offset=0; offset<LINEINST; offset++) {
            val<1> is_br = is_branch[offset];

            // Aggregate across ways to find if bimodal was used
            arr<val<1>,NUMWAYS> bim_used_arr = [&](u64 w){
                val<2> hc_t = hc_pred_type[w][offset];
                return (hc_t == val<2>{0});
            };
            val<1> any_bim_used = bim_used_arr.fo1().fold_or();

            // Check if prediction was correct
            arr<val<1>,NUMWAYS> correct_arr = [&](u64 w){
                return val<1>{hc_pred_val[w][offset]} == actualdirs[offset];
            };
            val<1> hc_correct = correct_arr.fo1().fold_or();

            // Count bimodal once per branch
            perf_hc_bim_count[0]    += static_cast<u64>(is_br & any_bim_used);
            perf_hc_bim_correct[0]  += static_cast<u64>(is_br & any_bim_used & hc_correct);
        }
#endif
#endif
        // Global history update
        val<1> line_end = inst_oh >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            next_pc.fanout(hard<2>{});
            path_history = (path_history << 1) ^ val<GHIST1>{next_pc>>2};
            gfolds.update(val<PATHBITS>{next_pc>>2});
        });

        num_branch = 0;
    }
    ~my_bp() {

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        print_perf_counters();
#endif
#endif
    }
};