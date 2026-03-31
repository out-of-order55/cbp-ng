// this is a basic TAGE, not necessarily well optimized

#define USE_ALT
#define RESET_UBITS

#if defined(GATE)
#error "GATE configuration has been removed from my_bp_v1."
#endif

#if defined(DEBUG_ENERGY)
#error "DEBUG_ENERGY configuration has been removed from my_bp_v1."
#endif

#if defined(MY_BP_V1_DISABLE_SC) || defined(MY_SC) || defined(SC_CFG_BASE) || defined(SC_CFG_MEDIUM) || defined(SC_CFG_FULL) || defined(SC_USE_BIAS) || defined(SC_FGEHL) || defined(SC_BGEHL) || defined(SC_USE_TAIMLI) || defined(SC_USE_BRIMLI) || defined(SC_USE_FHIST) || defined(SC_USE_BHIST) || defined(SC_GLOBAL_THRE_INIT) || defined(SC_GLOBAL_THRE_PIPE_STAGES)
#error "SC configuration has been removed from my_bp_v1."
#endif

#ifndef MY_BP_V1_PERF_ENABLE_TRACE_DB
#define MY_BP_V1_PERF_ENABLE_TRACE_DB 1
#endif

#ifndef MY_BP_V1_PERF_PRINT_TRACE
#define MY_BP_V1_PERF_PRINT_TRACE 1
#endif

#ifndef MY_BP_V1_P1_MICRO_TAGE
#define MY_BP_V1_P1_MICRO_TAGE 0
#endif

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include <iomanip>
#include <sstream>
#include <array>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdint>
#include <sqlite3.h>
#include "my_bp_v1_perf.hpp"

using namespace hcm;

template<u64 W, u64 M>
struct split_row_ram2 {
    static_assert(std::has_single_bit(M));
    static constexpr u64 A = std::bit_width(M - 1);
    ram<val<W>, M> lane[2];

    split_row_ram2(const char *label = "") : lane{label} {}

    arr<val<W>, 2> read(val<A> addr)
    {
        addr.fanout(hard<2>{});
        return arr<val<W>, 2>{[&](u64 slot) -> val<W> {
            return lane[slot].read(addr);
        }};
    }

    void write(val<A> addr, arr<val<W>, 2> row)
    {
        addr.fanout(hard<2>{});
        row.fanout(hard<2>{});
        lane[0].write(addr, row[0]);
        lane[1].write(addr, row[1]);
    }

    void reset()
    {
        lane[0].reset();
        lane[1].reset();
    }
};

template<u64 LOGLB=6, u64 NUMG=14, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=300, u64 LOGP1=14, u64 GHIST1=6,u64 LOGBANKS = 1,u64 LOGBIAS = 11>
struct my_bp_v1 : predictor {
    // provides 2^(LOGLB-2) predictions per cycle
    // P2 is a TAGE, P1 is a gshare
    static_assert(LOGLB>2);
    static_assert(NUMG>0);


    //TODO:need review
    static constexpr u64 PERCWIDTH = 6;
    static constexpr u64 NUMGEHL = 2;
    static constexpr u64 LOGGEHL = 10;

    static constexpr u64 MINHIST = 2;
    static constexpr u64 METABITS = 4;
    static constexpr u64 UCTRBITS = 8;
    static constexpr u64 PATHBITS = 6;

    static constexpr const char *GHYST_IMPL_NAME = "rwram";

#ifdef USE_ALT
    static constexpr u64 METAPIPE = 2;
#endif
    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGB > LOGLINEINST);
    static_assert(LOGP1 > LOGLINEINST);
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 bindex_bits = LOGB-LOGLINEINST;
    static_assert(TAGW > LOGLINEINST); // the unhashed line offset is part of the tag

    static constexpr u64 HTAGBITS = TAGW-LOGLINEINST; // hashed tag bits
    static_assert(NUMG == 14);
    static constexpr u64 NUMP = NUMG / 2;
    static constexpr u64 NUMG_GROUP = NUMP;
    static constexpr u64 LOGG_SH = LOGG;

#if MY_BP_V1_P1_MICRO_TAGE
    static constexpr u64 MT_NT = 4;
    static constexpr u64 MT_PIPE = 2;
    static_assert(MT_PIPE >= 1);
    static constexpr u64 MT_PIPE_LAST = MT_PIPE - 1;
    static constexpr u64 MT_LOGSETS = 9;
    static constexpr u64 MT_SETS = 1 << MT_LOGSETS;
    static constexpr u64 MT_TAGW = 12;
    static constexpr u64 MT_TAGW0 = MT_TAGW;
    static constexpr u64 MT_TAGW1 = MT_TAGW;
    static constexpr u64 MT_TAGW2 = MT_TAGW;
    static constexpr u64 MT_TAGW3 = MT_TAGW;
    static constexpr u64 MT_HTAGW0 = MT_TAGW0 - LOGLINEINST;
    static constexpr u64 MT_HTAGW1 = MT_TAGW1 - LOGLINEINST;
    static constexpr u64 MT_HTAGW2 = MT_TAGW2 - LOGLINEINST;
    static constexpr u64 MT_HTAGW3 = MT_TAGW3 - LOGLINEINST;
    static constexpr u64 MT_USE_MIN_PID = 2;
    static constexpr u64 MT_USE_MIN_HYST = 2;
#endif


    geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;

    reg<1> true_block = 1;

    // for P1
    reg<GHIST1> global_history1;
    reg<index1_bits> index1;
    reg<index1_bits> p1_idx_used;
    arr<reg<1>,LINEINST> readp1; // prediction bits read from P1 table for each offset
    reg<LINEINST> p1; // P1 predictions

#if MY_BP_V1_P1_MICRO_TAGE
    // P1 micro-TAGE (ahead pipe: cycle n-1 read, cycle n compare/select)
    geometric_folds<MT_NT,5,24,MT_LOGSETS,MT_HTAGW0> mt_gfolds;

    zone P1_MICRO_TAGE_CLUSTER;
    ram<val<MT_TAGW0>,MT_SETS> mt_tag[MT_NT] {"p1_mt_tag"};
    ram<val<1>,MT_SETS> mt_pred[MT_NT] {"p1_mt_pred"};
    rwram<2,MT_SETS,4> mt_hyst[MT_NT] {"p1_mt_hyst"};
    rwram<1,MT_SETS,4> mt_u[MT_NT] {"p1_mt_u"};

    reg<1> mt_pipe_valid[MT_PIPE];
    reg<64> mt_pipe_pc[MT_PIPE];
    reg<MT_TAGW0> mt_pipe_rtag[MT_NT][MT_PIPE];
    reg<1> mt_pipe_rpred[MT_NT][MT_PIPE];
    reg<2> mt_pipe_rhyst[MT_NT][MT_PIPE];
    reg<1> mt_pipe_ru[MT_NT][MT_PIPE];

    arr<reg<1>,LINEINST> mt_meta_hit;
    arr<reg<1>,LINEINST> mt_meta_tag_hit;
    arr<reg<2>,LINEINST> mt_meta_pid;
    arr<reg<1>,LINEINST> mt_meta_provider_pred;
    arr<reg<1>,LINEINST> mt_meta_provider_weak;
    arr<reg<1>,LINEINST> mt_meta_gshare_pred;
    arr<reg<MT_NT>,LINEINST> mt_meta_hitmask;
    arr<reg<MT_LOGSETS>,LINEINST> mt_meta_idx[MT_NT];
    arr<reg<MT_TAGW0>,LINEINST> mt_meta_tag[MT_NT];
    arr<reg<1>,LINEINST> mt_line_pred;
    arr<reg<1>,LINEINST> mt_meta_rpred[MT_NT];
    arr<reg<2>,LINEINST> mt_meta_rhyst[MT_NT];
    arr<reg<1>,LINEINST> mt_meta_ru[MT_NT];
    reg<8> mt_uclr_ctr = val<8>{1 << 7};
#endif





    arr<reg<1>,LINEINST> pred2; // alternate P2 prediction for each offset
    reg<LINEINST> p2; // final P2 predictions
    reg<LINEINST> tage_p2; // final P2 predictions
    arr<reg<NUMG+1>,LINEINST> match; // all matches for each offset
    arr<reg<NUMG+1>,LINEINST> match1; // longest match for each offset
    arr<reg<NUMG+1>,LINEINST> match2; // second longest match for each offset
    // for P2
    reg<bindex_bits> bindex; // bimodal table index
    arr<reg<LOGG>,NUMG> gindex; // global tables indexes
    arr<reg<HTAGBITS>,NUMG> htag; // computed hashed tags
    arr<reg<LOGG_SH>,NUMG_GROUP> group_row_idx; // grouped indexes (adjacent tables share row index)
    arr<reg<TAGW>,2> snap_tag[NUMG_GROUP]; // prediction-time snapshots, used by update RMW
    arr<reg<1>,2> snap_pred[NUMG_GROUP];
    arr<reg<2>,2> snap_hyst[NUMG_GROUP];
    arr<reg<1>,2> snap_u[NUMG_GROUP];
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
    my_bp_v1_perf::PerfState<NUMG> perf_event_state;
    std::array<std::array<u8, (1 << index1_bits)>, LINEINST> perf_shadow_p1_hyst = {};
    std::array<std::array<u8, (1 << bindex_bits)>, LINEINST> perf_shadow_bhyst = {};

    // Store prediction source per offset (set in predict2, read in update_cycle)
    // 0=bimodal, 1=provider, 2=alt
    arr<reg<2>,LINEINST> pred_source_stored;
    arr<reg<NUMG+1>,LINEINST> pred_match1_stored; // stored match1 for hit info
    arr<reg<NUMG+1>,LINEINST> pred_match2_stored; // stored match2 for alt info

    // Counter aliases: keep existing call sites unchanged while state storage lives in perf module.
#define perf_predictions perf_event_state.perf_predictions
#define perf_correct perf_event_state.perf_correct
#define perf_provider_used perf_event_state.perf_provider_used
#define perf_provider_correct perf_event_state.perf_provider_correct
#define perf_provider_wrong perf_event_state.perf_provider_wrong
#define perf_alt_used perf_event_state.perf_alt_used
#define perf_alt_correct perf_event_state.perf_alt_correct
#define perf_alt_wrong perf_event_state.perf_alt_wrong
#define perf_bimodal_used perf_event_state.perf_bimodal_used
#define perf_bimodal_correct perf_event_state.perf_bimodal_correct
#define perf_bimodal_wrong perf_event_state.perf_bimodal_wrong
#define perf_mispred_blame_tage perf_event_state.perf_mispred_blame_tage
#define perf_mispred_blame_p1 perf_event_state.perf_mispred_blame_p1
#define perf_table_reads perf_event_state.perf_table_reads
#define perf_table_hits perf_event_state.perf_table_hits
#define perf_table_alloc perf_event_state.perf_table_alloc
#define perf_alloc_failures perf_event_state.perf_alloc_failures
#define perf_alloc_fail_highest perf_event_state.perf_alloc_fail_highest
#define perf_alloc_fail_noubit perf_event_state.perf_alloc_fail_noubit
#define perf_extra_cycle_total perf_event_state.perf_extra_cycle_total
#define perf_extra_cycle_badpred perf_event_state.perf_extra_cycle_badpred
#define perf_extra_cycle_mispredict perf_event_state.perf_extra_cycle_mispredict
#define perf_extra_cycle_p1_update perf_event_state.perf_extra_cycle_p1_update
#define perf_phys_tag_w perf_event_state.perf_phys_tag_w
#define perf_phys_pred_w perf_event_state.perf_phys_pred_w
#define perf_phys_hyst_w perf_event_state.perf_phys_hyst_w
#define perf_phys_u_w perf_event_state.perf_phys_u_w
#define perf_hyst_skip_p1 perf_event_state.perf_hyst_skip_p1
#define perf_hyst_skip_bim perf_event_state.perf_hyst_skip_bim
#define perf_hyst_skip_tage perf_event_state.perf_hyst_skip_tage
#define perf_weak_wrong_p1 perf_event_state.perf_weak_wrong_p1
#define perf_weak_wrong_bim perf_event_state.perf_weak_wrong_bim
#define perf_weak_wrong_tage perf_event_state.perf_weak_wrong_tage
#define perf_conf perf_event_state.perf_conf
#define perf_mt_p1_slots perf_event_state.perf_mt_p1_slots
#define perf_mt_p1_match_p2 perf_event_state.perf_mt_p1_match_p2
#define perf_mt_p1_disagree_p2 perf_event_state.perf_mt_p1_disagree_p2
#define perf_mt_p1_from_micro perf_event_state.perf_mt_p1_from_micro
#define perf_mt_p1_from_gshare perf_event_state.perf_mt_p1_from_gshare
#define perf_mt_p1_micro_match_p2 perf_event_state.perf_mt_p1_micro_match_p2
#define perf_mt_p1_micro_disagree_p2 perf_event_state.perf_mt_p1_micro_disagree_p2
#define perf_mt_p1_gshare_match_p2 perf_event_state.perf_mt_p1_gshare_match_p2
#define perf_mt_p1_gshare_disagree_p2 perf_event_state.perf_mt_p1_gshare_disagree_p2
#define perf_mt_provider_hit perf_event_state.perf_mt_provider_hit
#define perf_mt_provider_hit_table perf_event_state.perf_mt_provider_hit_table
#define perf_mt_provider_use_table perf_event_state.perf_mt_provider_use_table
#define perf_mt_provider_correct_table perf_event_state.perf_mt_provider_correct_table
#define perf_mt_provider_wrong_table perf_event_state.perf_mt_provider_wrong_table
#define perf_mt_table_reads perf_event_state.perf_mt_table_reads
#define perf_mt_table_hits perf_event_state.perf_mt_table_hits
#define perf_mt_table_alloc perf_event_state.perf_mt_table_alloc
#define perf_mt_conf perf_event_state.perf_mt_conf
#define perf_mt_provider_wrong perf_event_state.perf_mt_provider_wrong
#define perf_mt_uclr_ctr_gate perf_event_state.perf_mt_uclr_ctr_gate
#define perf_mt_uclr_ctr_inc perf_event_state.perf_mt_uclr_ctr_inc
#define perf_mt_uclr_ctr_dec perf_event_state.perf_mt_uclr_ctr_dec
#define perf_mt_uclr_sat perf_event_state.perf_mt_uclr_sat
#define perf_mt_u_reset perf_event_state.perf_mt_u_reset
#define perf_mt_alloc_pick_table perf_event_state.perf_mt_alloc_pick_table
#define perf_mt_alloc_req_table perf_event_state.perf_mt_alloc_req_table
#define perf_mt_pred_req_table perf_event_state.perf_mt_pred_req_table
#define perf_mt_hyst_req_table perf_event_state.perf_mt_hyst_req_table
#define perf_mt_u_req_table perf_event_state.perf_mt_u_req_table
#define perf_mt_uclear_req_table perf_event_state.perf_mt_uclear_req_table
#define perf_mt_tag_w_table perf_event_state.perf_mt_tag_w_table
#define perf_mt_pred_w_table perf_event_state.perf_mt_pred_w_table
#define perf_mt_hyst_w_table perf_event_state.perf_mt_hyst_w_table
#define perf_mt_u_w_table perf_event_state.perf_mt_u_w_table

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
            "conf INTEGER"
            ");",
            nullptr, nullptr, nullptr);
        sqlite3_exec(trace_db, "BEGIN;", nullptr, nullptr, nullptr);
        sqlite3_prepare_v2(trace_db,
            "INSERT INTO trace VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
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
                      u64 conf) {
#if !MY_BP_V1_PERF_ENABLE_TRACE_DB
        (void)cycle; (void)pc; (void)offset;
        (void)actual_dir; (void)predicted_dir; (void)mispredict;
        (void)pred_source; (void)pred_table; (void)bim_index;
        (void)hit; (void)hit_table; (void)hit_gtag; (void)hit_gindex;
        (void)alloc; (void)alloc_table; (void)alloc_gindex; (void)alloc_tag;
        (void)conf;
        return;
#endif
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
        sqlite3_step(trace_stmt);
        if (trace_seq % 100000 == 0) {
            sqlite3_exec(trace_db, "COMMIT;", nullptr, nullptr, nullptr);
            sqlite3_exec(trace_db, "BEGIN;",  nullptr, nullptr, nullptr);
        }
    }

    void print_perf_counters() {
        u64 ghyst_reads = 0;
        u64 ghyst_writes = 0;
        u64 ghyst_stale_reads = 0;
        u64 ghyst_drop_writes = 0;
        u64 ghyst_write_hit_events = 0;
        for (u64 j = 0; j < NUMP; j++) {
            ghyst_reads += ghyst_p[j].perf_total_reads();
            ghyst_writes += ghyst_p[j].perf_total_writes();
            ghyst_stale_reads += ghyst_p[j].perf_total_stale_reads();
            ghyst_drop_writes += ghyst_p[j].perf_total_drop_writes();
            ghyst_write_hit_events += ghyst_p[j].perf_total_write_updates();
        }

        std::array<std::uint64_t, NUMG> hist_len = {};
        for (u64 j = 0; j < NUMG; j++) {
            hist_len[j] = gfolds.HLEN[j];
        }
        my_bp_v1_perf::PerfReportContext<NUMG> report_ctx{
            .hist_len = hist_len,
            .ghyst_impl_name = GHYST_IMPL_NAME,
            .ghyst_reads = ghyst_reads,
            .ghyst_writes = ghyst_writes,
            .ghyst_stale_reads = ghyst_stale_reads,
            .ghyst_drop_writes = ghyst_drop_writes,
            .ghyst_write_hit_events = ghyst_write_hit_events,
        };
        (void)my_bp_v1_perf::print_main_perf_report(std::cerr, perf_event_state, report_ctx);

        // Finalize SQLite trace (optional)
#if MY_BP_V1_PERF_ENABLE_TRACE_DB
        close_trace_db();
#if MY_BP_V1_PERF_PRINT_TRACE
        std::cerr << "│ Full exec trace written to: trace_v1.db                         │\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";
#endif
#endif

        // Write per-PC mispred summary + top-20 to stderr (optional)
#if MY_BP_V1_PERF_PRINT_TRACE
        {
            std::vector<std::pair<u64,MispredRecord>> sorted_db(mispred_db.begin(), mispred_db.end());
            std::sort(sorted_db.begin(), sorted_db.end(),
                [](const auto &a, const auto &b){ return a.second.count > b.second.count; });

#if MY_BP_V1_PERF_ENABLE_TRACE_DB
            std::cerr << "\n┌─ Top 20 Most Mispredicted PCs (full trace -> trace_v1.db) ���────────────────────────────────┐\n";
#else
            std::cerr << "\n┌─ Top 20 Most Mispredicted PCs ───────────────────────────────────────────────────────────────┐\n";
#endif
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
#endif
    }

    void perf_count_p1_vs_p2(
        arr<val<1>,LINEINST> is_branch,
        arr<val<1>,LINEINST> p1_line,
        arr<val<1>,LINEINST> p2_line)
    {
        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

            bool p1_match_p2 = static_cast<bool>(p1_line[offset] == p2_line[offset]);
            perf_mt_p1_slots++;
            if (p1_match_p2) perf_mt_p1_match_p2++;
            else perf_mt_p1_disagree_p2++;

#if MY_BP_V1_P1_MICRO_TAGE
            bool from_micro = static_cast<bool>(val<1>{mt_meta_hit[offset]});
            bool tag_hit = static_cast<bool>(val<1>{mt_meta_tag_hit[offset]});
            for (u64 t = 0; t < MT_NT; t++) {
                perf_mt_table_reads[t]++;
            }
            if (tag_hit) {
                u64 hitmask = static_cast<u64>(val<MT_NT>{mt_meta_hitmask[offset]});
                for (u64 t = 0; t < MT_NT; t++) {
                    if ((hitmask >> t) & 1ULL) {
                        perf_mt_table_hits[t]++;
                    }
                }
            }
            if (from_micro) {
                perf_mt_p1_from_micro++;
                perf_mt_provider_hit++;
                if (p1_match_p2) perf_mt_p1_micro_match_p2++;
                else perf_mt_p1_micro_disagree_p2++;

                u64 pid = static_cast<u64>(val<2>{mt_meta_pid[offset]});
                if (pid < MT_NT) {
                    perf_mt_provider_hit_table[pid]++;
                    perf_mt_provider_use_table[pid]++;
                    if (p1_match_p2) perf_mt_provider_correct_table[pid]++;
                    else perf_mt_provider_wrong_table[pid]++;

                    u64 conf = static_cast<u64>(val<2>{mt_meta_rhyst[pid][offset]});
                    if (conf < 4) perf_mt_conf[pid][conf]++;
                }
            } else {
                perf_mt_p1_from_gshare++;
                if (p1_match_p2) perf_mt_p1_gshare_match_p2++;
                else perf_mt_p1_gshare_disagree_p2++;
            }
#else
            perf_mt_p1_from_gshare++;
            if (p1_match_p2) perf_mt_p1_gshare_match_p2++;
            else perf_mt_p1_gshare_disagree_p2++;
#endif
        }
    }

#if MY_BP_V1_P1_MICRO_TAGE
    void perf_count_mt_update_core(
        arr<val<1>,LINEINST> mt_do_update,
        arr<val<1>,LINEINST> mt_provider_used,
        arr<val<1>,LINEINST> mt_provider_wrong,
        arr<val<MT_NT>,LINEINST> mt_pick_bits,
        arr<val<1>,LINEINST> mt_ctr_gate,
        arr<val<1>,LINEINST> mt_gshare_wrong,
        val<1> mt_uclr_sat)
    {
        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(mt_do_update[offset])) continue;
            if (static_cast<bool>(mt_provider_used[offset] & mt_provider_wrong[offset])) perf_mt_provider_wrong++;

            if (static_cast<bool>(mt_ctr_gate[offset])) {
                perf_mt_uclr_ctr_gate++;
                if (static_cast<bool>(mt_gshare_wrong[offset])) perf_mt_uclr_ctr_dec++;
                else perf_mt_uclr_ctr_inc++;
            }

            u64 pick_bits = static_cast<u64>(mt_pick_bits[offset]);
            for (u64 t = 0; t < MT_NT; t++) {
                if ((pick_bits >> t) & 1ULL) {
                    perf_mt_alloc_pick_table[t]++;
                    break;
                }
            }
        }

        if (static_cast<bool>(mt_uclr_sat)) {
            perf_mt_uclr_sat++;
            perf_mt_u_reset++;
        }
    }

    template<u64 I>
    void perf_count_mt_table_update(
        arr<val<1>,LINEINST> mt_alloc_req,
        arr<val<1>,LINEINST> mt_pred_req,
        arr<val<1>,LINEINST> mt_hyst_req,
        arr<val<1>,LINEINST> mt_u_req,
        arr<val<1>,LINEINST> mt_uclear_req,
        val<1> mt_tag_we,
        val<1> mt_pred_we,
        val<1> mt_hyst_we,
        val<1> mt_u_we)
    {
        static_assert(I < MT_NT);
        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (static_cast<bool>(mt_alloc_req[offset])) {
                perf_mt_alloc_req_table[I]++;
                perf_mt_table_alloc[I]++;
            }
            if (static_cast<bool>(mt_pred_req[offset])) perf_mt_pred_req_table[I]++;
            if (static_cast<bool>(mt_hyst_req[offset])) perf_mt_hyst_req_table[I]++;
            if (static_cast<bool>(mt_u_req[offset])) perf_mt_u_req_table[I]++;
            if (static_cast<bool>(mt_uclear_req[offset])) perf_mt_uclear_req_table[I]++;
        }
        if (static_cast<bool>(mt_tag_we)) perf_mt_tag_w_table[I]++;
        if (static_cast<bool>(mt_pred_we)) perf_mt_pred_w_table[I]++;
        if (static_cast<bool>(mt_hyst_we)) perf_mt_hyst_w_table[I]++;
        if (static_cast<bool>(mt_u_we)) perf_mt_u_w_table[I]++;
    }
#endif

#ifdef CHEATING_MODE
    template<typename TArray>
    void perf_count_phys_write(TArray &counter, u64 phys, val<1> we)
    {
        if (static_cast<bool>(we)) {
            counter[phys]++;
        }
    }

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

    void perf_count_extra_cycle(val<1> extra_cycle, val<1> some_badpred1, val<1> mispredict, val<1> p1_update) {
        perf_extra_cycle_total      += static_cast<u64>(extra_cycle);
        perf_extra_cycle_badpred    += static_cast<u64>(some_badpred1);
        perf_extra_cycle_mispredict += static_cast<u64>(mispredict);
        perf_extra_cycle_p1_update  += static_cast<u64>(p1_update);
    }

    void perf_count_hyst_skip_events(
        arr<val<1>,LINEINST> &is_branch,
        arr<val<1>,LINEINST> &disagree,
        arr<val<1>,LINEINST> &primary_wrong,
        arr<val<1>,LINEINST> &bhyst_bim_primary,
        arr<val<1>,NUMG> &primary,
        arr<val<1>,NUMG> &allocate,
        arr<val<1>,NUMG> &badpred1,
        arr<val<1>,NUMG> &g_sat)
    {
        const u64 p1_idx = static_cast<u64>(p1_idx_used);
        const u64 bim_idx = static_cast<u64>(bindex);

        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

            const bool p1_correct = !static_cast<bool>(disagree[offset]);
            if (p1_correct && perf_shadow_p1_hyst[offset][p1_idx] != 0) {
                perf_hyst_skip_p1++;
            }

            const bool bim_primary = static_cast<bool>(bhyst_bim_primary[offset]);
            const bool bim_correct = !static_cast<bool>(primary_wrong[offset]);
            if (bim_primary && bim_correct && perf_shadow_bhyst[offset][bim_idx] != 0) {
                perf_hyst_skip_bim++;
            }
        }

        for (u64 i = 0; i < NUMG; i++) {
            const bool tage_skip =
                static_cast<bool>(primary[i]) &&
                !static_cast<bool>(allocate[i]) &&
                !static_cast<bool>(badpred1[i]) &&
                static_cast<bool>(g_sat[i]);
            if (tage_skip) {
                perf_hyst_skip_tage++;
            }
        }
    }

    void perf_count_weak_wrong_events(
        arr<val<1>,LINEINST> is_branch,
        arr<val<1>,LINEINST> p1_weak,
        arr<val<1>,LINEINST> b_weak,
        arr<val<1>,NUMG> g_weak)
    {
        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

            if (static_cast<bool>(p1_weak[offset])) {
                perf_weak_wrong_p1++;
            }

            if (static_cast<bool>(b_weak[offset])) {
                perf_weak_wrong_bim++;
            }
        }

        for (u64 i = 0; i < NUMG; i++) {
            if (static_cast<bool>(g_weak[i])) {
                perf_weak_wrong_tage++;
            }
        }
    }


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
            u64 selected_table = NUMG;
            bool has_selected_table = false;

            if (src == 2) {
                for (u64 j=0; j<NUMG; j++) {
                    if (static_cast<u64>(alt_mask >> j) & 1) {
                        perf_alt_used[j]++;
                        selected_table = j;
                        has_selected_table = true;
                        break;
                    }
                }
            } else if (src == 1) {
                for (u64 j=0; j<NUMG; j++) {
                    if (static_cast<u64>(prov_mask >> j) & 1) {
                        perf_provider_used[j]++;
                        selected_table = j;
                        has_selected_table = true;
                        break;
                    }
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

                perf_mispred_blame_tage++;
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
            u64 pred_table  = has_selected_table ? selected_table : NUMG;

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
                insert_trace(
                    static_cast<u64>(panel.cycle), pc_val, offset,
                    static_cast<u64>(actual), static_cast<u64>(predicted),
                    is_misp & is_last,
                    pred_source, pred_table, static_cast<u64>(bindex),
                    hit_found, hit_table, hit_gtag, hit_gindex,
                    is_misp & is_last & alloc_found,
                    alloc_table, alloc_gindex_val, alloc_tag_val,
                    conf_val);
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

    }
#endif
#endif

    
    // Cluster P1 + bpred tables.
    zone P1_CLUSTER;
    // P1 (gshare)
    ram<val<1>,(1<<index1_bits)> table1_pred[LINEINST] {"P1 pred"}; // P1 prediction bit
    ram<val<1>,(1<<bindex_bits)> bim[LINEINST] {"bpred"}; // bimodal prediction bits
    // Put TAGE tables in a dedicated cluster.
    zone TAGE_CLUSTER;
    // P2 (TAGE)
    split_row_ram2<TAGW,(1<<LOGG_SH)> gtag_p[NUMP] {"tags"}; // tags (2-slot row, split 2xRAM)
    ram<arr<val<1>,2>,(1<<LOGG_SH)> gpred_p[NUMP] {"gpred"}; // predictions (2-slot row)
    rwram<4,(1<<LOGG_SH),4> ghyst_p[NUMP] {"ghyst"}; // packed hysteresis row: {slot1,slot0}
    rwram<2,(1<<LOGG_SH),4> ubit_p[NUMP] {"u"}; // packed u-bit row: {slot1,slot0}
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
#if MY_BP_V1_PERF_ENABLE_TRACE_DB
        open_trace_db();
#endif
#endif
    }


    ~my_bp_v1() {
#ifdef PERF_COUNTERS
        print_perf_counters();
#endif
    }

    void new_block(val<64> inst_pc)
    {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        block_entry = offset.fo1().decode().concat();
        block_entry.fanout(hard<3*LINEINST>{});
        block_size = 1;
    }

    val<index1_bits> p1_compute_index(val<64> inst_pc)
    {
        val<std::max(index1_bits,GHIST1)> lineaddr = inst_pc >> LOGLB;
        if constexpr (GHIST1 <= index1_bits) {
            return lineaddr ^ (val<index1_bits>{global_history1} << (index1_bits-GHIST1));
        } else {
            return global_history1.make_array(val<index1_bits>{}).append(lineaddr).fold_xor();
        }
    }

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<3>{});
        new_block(inst_pc);

        val<64-LOGLB> cur_line_pc = inst_pc >> LOGLB;
        // cur_line_pc.fanout(hard<4>{});
        val<index1_bits> cur_idx = p1_compute_index(inst_pc);
        cur_idx.fanout(hard<LINEINST + 1>{});
        // index1 = cur_idx;
        for (u64 offset=0; offset<LINEINST; offset++) {
            readp1[offset] = table1_pred[offset].read(cur_idx);
        }
        val<LINEINST> p1_vec = readp1.concat();
        p1_idx_used = cur_idx;
#if MY_BP_V1_P1_MICRO_TAGE
        if constexpr (MT_PIPE > 1) {
            for (u64 p = MT_PIPE - 1; p != 0; p--) {
                mt_pipe_valid[p] = mt_pipe_valid[p - 1];
                mt_pipe_pc[p] = mt_pipe_pc[p - 1];
            }
            static_loop<MT_NT>([&]<u64 J>(){
                for (u64 p = MT_PIPE - 1; p != 0; p--) {
                    mt_pipe_rtag[J][p] = mt_pipe_rtag[J][p - 1];
                    mt_pipe_rpred[J][p] = mt_pipe_rpred[J][p - 1];
                    mt_pipe_rhyst[J][p] = mt_pipe_rhyst[J][p - 1];
                    mt_pipe_ru[J][p] = mt_pipe_ru[J][p - 1];
                }
            });
        }
        val<LINEINST> gshare_line = p1_vec;
        auto gshare_split = gshare_line.make_array(val<1>{});
        val<64> line_base_pc = val<64>{cur_line_pc} << LOGLB;
        val<1> pipe_line_valid = val<1>{mt_pipe_valid[MT_PIPE_LAST]};
        line_base_pc.fanout(hard<LINEINST>{});
        pipe_line_valid.fanout(hard<LINEINST>{});
        static_loop<LINEINST>([&]<u64 I>(){
            val<LOGLINEINST> mt_offset = hard<I>{};
            val<64> mt_inst_pc = line_base_pc | (val<64>{mt_offset} << 2);
            val<std::max(index1_bits,LOGG)> mt_lineaddr = mt_inst_pc >> LOGLB;
            mt_lineaddr.fanout(hard<8>{});
            val<MT_LOGSETS> mt_pc9 = mt_lineaddr;
            val<MT_HTAGW0> mt_pctag = val<MT_HTAGW0>{mt_lineaddr}.reverse();
            arr<val<MT_LOGSETS>,MT_NT> mt_idx = [&](u64 j) -> val<MT_LOGSETS> {
                return mt_pc9 ^ mt_gfolds.template get<0>(j);
            };
            arr<val<MT_TAGW0>,MT_NT> mt_tag_vec = [&](u64 j) -> val<MT_TAGW0> {
                return concat(mt_offset, mt_pctag ^ mt_gfolds.template get<1>(j));
            };
            arr<val<1>,MT_NT> mt_hit = [&](u64 j) -> val<1> {
                return pipe_line_valid & (val<MT_TAGW0>{mt_pipe_rtag[j][MT_PIPE_LAST]} == mt_tag_vec[j]);
            };
            val<MT_NT> mt_hitmask = mt_hit.concat();
            mt_hitmask.fanout(hard<2>{});
            val<1> mt_any_hit = mt_hitmask != hard<0>{};
            val<MT_NT> mt_pri_hitmask = mt_hitmask.reverse().one_hot().reverse();
            arr<val<1>,MT_NT> mt_pri_hit = mt_pri_hitmask.make_array(val<1>{});
            val<2> mt_pid = arr<val<2>,MT_NT>{[&](u64 j) -> val<2> {
                return select(mt_pri_hit[j], val<2>{j}, val<2>{0});
            }}.fold_or();
            val<2> mt_hyst_pick = arr<val<2>,MT_NT>{[&](u64 j) -> val<2> {
                return select(mt_pri_hit[j], val<2>{mt_pipe_rhyst[j][MT_PIPE_LAST]}, val<2>{0});
            }}.fold_or();
            val<1> mt_pred_pick = arr<val<1>,MT_NT>{[&](u64 j) -> val<1> {
                return select(mt_pri_hit[j], val<1>{mt_pipe_rpred[j][MT_PIPE_LAST]}, val<1>{0});
            }}.fold_or();
            val<2> mt_sel_hyst = select(mt_any_hit, mt_hyst_pick, val<2>{mt_pipe_rhyst[0][MT_PIPE_LAST]});
            val<1> mt_pred = select(mt_any_hit, mt_pred_pick, val<1>{mt_pipe_rpred[0][MT_PIPE_LAST]});
            val<1> mt_weak = mt_sel_hyst == hard<0>{};
            val<1> mt_allow_pid = mt_pid >= hard<MT_USE_MIN_PID>{};
            val<1> mt_allow_hyst = mt_sel_hyst >= hard<MT_USE_MIN_HYST>{};
            (void)mt_allow_hyst;
            val<1> mt_use_micro = mt_any_hit & mt_allow_pid;

            mt_meta_tag_hit[I] = mt_any_hit;
            mt_meta_hit[I] = mt_use_micro;
            mt_meta_pid[I] = mt_pid;
            mt_meta_provider_pred[I] = mt_pred;
            mt_meta_provider_weak[I] = mt_weak;
            mt_meta_gshare_pred[I] = gshare_split[I];
            mt_meta_hitmask[I] = mt_hitmask;
            static_loop<MT_NT>([&]<u64 J>(){
                mt_meta_idx[J][I] = mt_idx[J];
                mt_meta_tag[J][I] = mt_tag_vec[J];
                mt_meta_rpred[J][I] = val<1>{mt_pipe_rpred[J][MT_PIPE_LAST]};
                mt_meta_rhyst[J][I] = val<2>{mt_pipe_rhyst[J][MT_PIPE_LAST]};
                mt_meta_ru[J][I] = val<1>{mt_pipe_ru[J][MT_PIPE_LAST]};
            });
            mt_line_pred[I] = select(mt_use_micro, mt_pred, gshare_split[I]);
        });
        p1 = mt_line_pred.concat();

        val<std::max(index1_bits,LOGG)> mt_issue_lineaddr = inst_pc >> LOGLB;
        mt_issue_lineaddr.fanout(hard<8>{});
        val<MT_LOGSETS> mt_issue_pc9 = mt_issue_lineaddr;
        mt_pipe_pc[0] = inst_pc;
        static_loop<MT_NT>([&]<u64 J>(){
            val<MT_LOGSETS> mt_issue_idx = mt_issue_pc9 ^ mt_gfolds.template get<0>(J);
            mt_pipe_rtag[J][0] = mt_tag[J].read(mt_issue_idx);
            mt_pipe_rpred[J][0] = mt_pred[J].read(mt_issue_idx);
            mt_pipe_rhyst[J][0] = mt_hyst[J].read(mt_issue_idx);
            mt_pipe_ru[J][0] = mt_u[J].read(mt_issue_idx);
        });
#else
        p1 = p1_vec;
#endif
        p1.fanout(hard<LINEINST>{});
        val<1> p1_taken = (block_entry & p1) != hard<0>{};
        return p1_taken;
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        val<1> p1_taken = ((block_entry<<block_size) & p1) != hard<0>{};
        return p1_taken;
    };
    void read_tag_pred_physical()
    {
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<LOGG_SH> row = val<LOGG_SH>{group_row_idx[G]};
            row.fanout(hard<4>{});

            arr<val<TAGW>,2> tag = gtag_p[G].read(row);
            arr<val<1>,2> pred = gpred_p[G].read(row);
            val<4> hyst_raw = ghyst_p[G].read(row);
            val<2> u_raw = ubit_p[G].read(row);

            hyst_raw.fanout(hard<2>{});
            u_raw.fanout(hard<2>{});
            arr<val<2>,2> hyst = {val<2>{hyst_raw}, val<2>{hyst_raw >> 2}};
            arr<val<1>,2> u = {val<1>{u_raw}, val<1>{u_raw >> 1}};

            static_loop<2>([&]<u64 SLOT>(){
                snap_tag[G][SLOT] = tag[SLOT];
                snap_pred[G][SLOT] = pred[SLOT];
                snap_hyst[G][SLOT] = hyst[SLOT];
                snap_u[G][SLOT] = u[SLOT];
            });
        });

        static_loop<NUMG_GROUP>([&]<u64 G>(){
            snap_tag[G].fanout(hard<2>{});
            snap_pred[G].fanout(hard<2>{});
            snap_hyst[G].fanout(hard<2>{});
            snap_u[G].fanout(hard<2>{});
            readt[2*G] = snap_tag[G][0];
            readt[2*G+1] = snap_tag[G][1];
            readc[2*G] = snap_pred[G][0];
            readc[2*G+1] = snap_pred[G][1];
            readh[2*G] = snap_hyst[G][0];
            readh[2*G+1] = snap_hyst[G][1];
            readu[2*G] = snap_u[G][0];
            readu[2*G+1] = snap_u[G][1];
        });
    }

    void write_tag_physical(
        arr<val<1>,NUMG_GROUP> &tag_we,
        arr<val<TAGW>,NUMG_GROUP> &tag_row0,
        arr<val<TAGW>,NUMG_GROUP> &tag_row1)
    {
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = tag_we[G];
            arr<val<TAGW>,2> row = {tag_row0[G], tag_row1[G]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_tag_w, G, we);
#endif
#endif
            execute_if(we, [&](){ gtag_p[G].write(group_row_idx[G], row); });
        });
    }

    void write_pred_physical(
        arr<val<1>,NUMG_GROUP> &pred_we,
        arr<val<1>,NUMG_GROUP> &pred_row0,
        arr<val<1>,NUMG_GROUP> &pred_row1)
    {
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = pred_we[G];
            arr<val<1>,2> row = {pred_row0[G], pred_row1[G]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_pred_w, G, we);
#endif
#endif
            execute_if(we, [&](){ gpred_p[G].write(group_row_idx[G], row); });
        });
    }

    void write_hyst_physical(
        arr<val<1>,NUMG_GROUP> &hyst_we,
        arr<val<2>,NUMG_GROUP> &hyst_row0,
        arr<val<2>,NUMG_GROUP> &hyst_row1,
        val<1> noconflict)
    {
        noconflict.fanout(hard<NUMG_GROUP>{});
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = hyst_we[G];
            val<4> row = concat(hyst_row1[G], hyst_row0[G]);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_hyst_w, G, we);
#endif
#endif
            ghyst_p[G].write(group_row_idx[G], row, we, noconflict);
        });
    }

    void write_ubit_physical(
        arr<val<1>,NUMG_GROUP> &u_we,
        arr<val<1>,NUMG_GROUP> &u_row0,
        arr<val<1>,NUMG_GROUP> &u_row1,
        val<1> noconflict)
    {
        noconflict.fanout(hard<NUMG_GROUP>{});
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = u_we[G];
            val<2> row = concat(u_row1[G], u_row0[G]);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_u_w, G, we);
#endif
#endif
            ubit_p[G].write(group_row_idx[G], row, we, noconflict);
        });
    }

    void tage_pred(val<64> inst_pc){
        val<std::max(bindex_bits,LOGG)> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<1+NUMG*2>{});
        gfolds.fanout(hard<2>{});
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
        gindex.fanout(hard<8>{});
        for (u64 g=0; g<NUMG_GROUP; g++) {
            group_row_idx[g] = val<LOGG_SH>{gindex[2*g]};
        }
        group_row_idx.fanout(hard<2>{});
        read_tag_pred_physical();


        readt.fanout(hard<LINEINST+1>{});
        readc.fanout(hard<3>{});
        readu.fanout(hard<2>{});
        notumask = ~readu.concat();
        notumask.fanout(hard<2>{});

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
        match1.fanout(hard<3>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
        }
        pred1.fanout(hard<3>{});

        // for each offset, find second longest match and select secondary prediction
        for (u64 offset=0; offset<LINEINST; offset++) {
            match2[offset] = (match[offset]^match1[offset]).one_hot();
        }
        match2.fanout(hard<3>{});
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
    }
    val<1> predict2(val<64> inst_pc)
    {
        tage_pred(inst_pc);
        p2 = tage_p2;
        p2.fanout(hard<LINEINST>{});
        val<1> taken = (block_entry & p2) != hard<0>{};
        // taken.fanout(hard<2>{});

        reuse_prediction(~val<1>{block_entry>>(LINEINST-1)});
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = ((block_entry<<block_size) & p2) != hard<0>{};
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
#if MY_BP_V1_P1_MICRO_TAGE
                mt_gfolds.update(val<PATHBITS>{next_pc>>2});
#endif
                true_block = 1;
            });
            return; // stop here
        }
        mispredict.fanout(hard<3>{});
        val<1> correct_pred = ~mispredict;
        correct_pred.fanout(hard<NUMG+2>{});
        p1_idx_used.fanout(hard<3*LINEINST>{});
        p2.fanout(hard<2>{});
        p1.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST+2>{});
        // gindex.fanout(hard<8>{});
        htag.fanout(hard<2>{});
        // readb.fanout(hard<2>{});
        readt.fanout(hard<4>{});
        readc.fanout(hard<2>{});
        readh.fanout(hard<3>{});
        //SC fix
        match1.fanout(hard<2>{});
        group_row_idx.fanout(hard<4>{});
        // match2.fanout(hard<2>{});
        pred1.fanout(hard<2>{});
        pred2.fanout(hard<2+NUMG>{});
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<2>{});
        gfolds.fanout(hard<2>{});
#ifdef USE_ALT
        meta.fanout(hard<2>{});
#endif
        val<LOGLINEINST> last_offset = branch_offset[num_branch-1];
        last_offset.fanout(hard<5*NUMG>{});

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
        primary.fanout(hard<4>{});

        arr<val<1>,LINEINST> primary_wrong = [&](u64 offset){
            return pred1[offset] != branch_taken[offset];
        };
        primary_wrong.fanout(hard<2>{});
        arr<val<1>,LINEINST> bhyst_bim_primary = [&](u64 offset){ return actual_match1[offset] >> NUMG; };
        bhyst_bim_primary.fanout(hard<3>{});

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
        badpred1.fanout(hard<4>{});

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
        disagree_mask.fanout(hard<2+2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});
        // Read P1 hysteresis before extra_cycle scheduling so read/write can be separated by extra cycle.
        arr<val<1>,LINEINST> p1_weak = [&](u64 offset) -> val<1> {
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(p1_idx_used);
            });
        };
        p1_weak.fanout(hard<4>{});

        // read the bimodal hysteresis if bimodal caused a misprediction
        arr<val<1>,LINEINST> b_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff cause of misprediction and hysteresis is weak
            return execute_if(bhyst_bim_primary[offset] & primary_wrong[offset], [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
        };
        b_weak.fanout(hard<4>{});

        // determine which primary global predictions are incorrect with a weak hysteresis
        arr<val<1>,NUMG> g_weak = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & badpred1[i] & (readh[i]==hard<0>{});
        };
        g_weak.fanout(hard<5>{});

        arr<val<1>,NUMG> g_sat = [&](u64 i) {
            return readh[i]==hard<3>{};
        };

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_hyst_skip_events(
            is_branch,
            disagree,
            primary_wrong,
            bhyst_bim_primary,
            primary,
            allocate,
            badpred1,
            g_sat);
#endif
#endif

        // need extra cycle for modifying prediction bits and for TAGE allocation
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        // some_badpred1 is consumed by extra_cycle and PERF counters in this mode.
        some_badpred1.fanout(hard<2>{});
#endif
#endif

        val<1> extra_cycle_base = some_badpred1 | mispredict | (disagree_mask != hard<0>{});
        val<1> extra_cycle = extra_cycle_base;

        extra_cycle.fanout(hard<NUMG*2+8>{});
        need_extra_cycle(extra_cycle);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        val<1> p1_update = disagree_mask != hard<0>{};
        perf_count_extra_cycle(extra_cycle, some_badpred1, mispredict, p1_update);
#endif
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
        arr<val<1>,NUMG_GROUP> tag_we = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            u64 i1 = i0 + 1;
            return allocate[i0] | allocate[i1];
        }};
        arr<val<TAGW>,NUMG_GROUP> tag_row0 = arr<val<TAGW>,NUMG_GROUP>{[&](u64 g) -> val<TAGW> {
            u64 i0 = 2 * g;
            return select(allocate[i0], concat(last_offset, htag[i0]), val<TAGW>{snap_tag[g][0]});
        }};
        arr<val<TAGW>,NUMG_GROUP> tag_row1 = arr<val<TAGW>,NUMG_GROUP>{[&](u64 g) -> val<TAGW> {
            u64 i1 = 2 * g + 1;
            return select(allocate[i1], concat(last_offset, htag[i1]), val<TAGW>{snap_tag[g][1]});
        }};

        write_tag_physical(tag_we, tag_row0, tag_row1);


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
        update_u.fanout(hard<4>{});
        val<NUMG> uclearmask = postmask & noalloc.replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<6>{});
        allocate.fanout(hard<6>{});
        goodpred.fanout(hard<2>{});
        arr<val<1>,NUMG_GROUP> u_we = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            u64 i1 = i0 + 1;
            return update_u[i0] | allocate[i0] | uclear[i0] | update_u[i1] | allocate[i1] | uclear[i1];
        }};
        arr<val<1>,NUMG_GROUP> u_row0 = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            val<1> we0 = update_u[i0] | allocate[i0] | uclear[i0];
            val<1> newu0 = goodpred[i0] & ~allocate[i0] & ~uclear[i0];
            return select(we0, newu0, val<1>{snap_u[g][0]});
        }};
        arr<val<1>,NUMG_GROUP> u_row1 = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i1 = 2 * g + 1;
            val<1> we1 = update_u[i1] | allocate[i1] | uclear[i1];
            val<1> newu1 = goodpred[i1] & ~allocate[i1] & ~uclear[i1];
            return select(we1, newu1, val<1>{snap_u[g][1]});
        }};
        // u_we.fanout(hard<2>{});
        // u_row0.fanout(hard<2>{});
        // u_row1.fanout(hard<2>{});
        write_ubit_physical(u_we, u_row0, u_row1, extra_cycle);

        auto p1_split = p1.make_array(val<1>{});
        p1_split.fanout(hard<4>{});
        auto p2_split = p2.make_array(val<1>{});
        p2_split.fanout(hard<4>{});
#ifdef PERF_COUNTERS
        perf_count_p1_vs_p2(is_branch, p1_split, p2_split);
        perf_count_weak_wrong_events(is_branch, p1_weak, b_weak, g_weak);
#endif

        // P1 update policy:
        // - update prediction when P1 hysteresis is weak and P1 source is gshare;
        // - update hysteresis for every resolved branch.
        arr<val<1>,LINEINST> p1_hyst_write_req = [&](u64 offset) -> val<1> {
            return is_branch[offset] & ~p1_weak[offset];
        };
        p1_hyst_write_req.fanout(hard<3>{});
        arr<val<1>,LINEINST> p1_pred_write_req = [&](u64 offset) -> val<1> {
#if MY_BP_V1_P1_MICRO_TAGE
            val<1> p1_from_gshare = ~val<1>{mt_meta_hit[offset]};
#else
            val<1> p1_from_gshare = val<1>{1};
#endif
            return is_branch[offset] & ~p1_hyst_write_req[offset] & p1_from_gshare;
        };
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(p1_hyst_write_req[offset], [&](){
                table1_hyst[offset].write(p1_idx_used, ~disagree[offset]);
            });
            execute_if(p1_pred_write_req[offset], [&](){
                table1_pred[offset].write(p1_idx_used, p2_split[offset]);
            });
        }
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        {
            const u64 p1_idx = static_cast<u64>(p1_idx_used);
            for (u64 offset=0; offset<LINEINST; offset++) {
                if (!static_cast<bool>(p1_hyst_write_req[offset])) continue;
                perf_shadow_p1_hyst[offset][p1_idx] = static_cast<u8>(static_cast<bool>(~disagree[offset]));
            }
        }
#endif
#endif

#if MY_BP_V1_P1_MICRO_TAGE
        // Micro-TAGE update on all resolved branch offsets in this block.
        arr<val<1>,LINEINST> mt_do_update = is_branch;
        arr<val<1>,LINEINST> mt_provider_hit = mt_meta_tag_hit;
        arr<val<1>,LINEINST> mt_provider_used = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return val<1>{mt_meta_hit[offset]};
        }};
        arr<val<2>,LINEINST> mt_provider_id = mt_meta_pid;
        arr<val<1>,LINEINST> mt_provider_pred = mt_meta_provider_pred;
        arr<val<1>,LINEINST> mt_provider_weak = mt_meta_provider_weak;
        arr<val<1>,LINEINST> mt_gshare_pred = mt_meta_gshare_pred;
        arr<val<1>,LINEINST> mt_p1_pred = p1_split;
        arr<val<1>,LINEINST> mt_p2_pred = p2_split;
        arr<val<1>,LINEINST> mt_target_dir = mt_p2_pred;

        mt_do_update.fanout(hard<12>{});
        mt_provider_hit.fanout(hard<12>{});
        mt_provider_used.fanout(hard<5>{});
        mt_provider_id.fanout(hard<8>{});
        mt_provider_pred.fanout(hard<8>{});
        mt_provider_weak.fanout(hard<6>{});
        mt_gshare_pred.fanout(hard<6>{});
        mt_p1_pred.fanout(hard<6>{});
        mt_p2_pred.fanout(hard<6>{});
        mt_target_dir.fanout(hard<8>{});

        arr<val<1>,LINEINST> mt_provider_wrong = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return mt_do_update[offset] & (mt_p1_pred[offset] != mt_p2_pred[offset]);
        }};
        arr<val<1>,LINEINST> mt_new_u = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return mt_provider_pred[offset] == mt_p2_pred[offset];
        }};
        arr<val<1>,LINEINST> mt_update_u = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return (mt_provider_pred[offset] != mt_p2_pred[offset]) | (mt_provider_pred[offset] != mt_gshare_pred[offset]);
        }};

        mt_provider_wrong.fanout(hard<11>{});
        mt_new_u.fanout(hard<6>{});
        mt_update_u.fanout(hard<6>{});

        val<1> mt_clear_pipe_valid = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return mt_provider_used[offset] & mt_provider_wrong[offset];
        }}.fold_or();
        mt_pipe_valid[0] = ~mt_clear_pipe_valid;

        // gshare-referenced global u-clear counter update.
        arr<val<1>,LINEINST> mt_ctr_gate = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return mt_do_update[offset] & mt_provider_hit[offset] &
                   (mt_provider_pred[offset] != mt_gshare_pred[offset]) &
                   mt_provider_weak[offset];
        }};
        arr<val<1>,LINEINST> mt_gshare_wrong = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
            return mt_gshare_pred[offset] != mt_target_dir[offset];
        }};
        arr<val<2,i64>,LINEINST> mt_uclr_incr = arr<val<2,i64>,LINEINST>{[&](u64 offset) -> val<2,i64> {
            val<1> dec_sign = mt_gshare_wrong[offset];
            return select(mt_ctr_gate[offset], concat(dec_sign.fo1(), val<1>{1}), val<2>{0});
        }};
        auto mt_uclr_vote = mt_uclr_incr.fo1().fold_add();
        auto mt_uclr_zero = val<decltype(mt_uclr_vote)::size, i64>{0};
        val<1> mt_uclr_upd = mt_uclr_vote != mt_uclr_zero;
        val<1> mt_uclr_inc = mt_uclr_vote > mt_uclr_zero;
        val<8> mt_uclr_base = val<8>{mt_uclr_ctr};
        val<8> mt_uclr_updated = select(mt_uclr_upd, update_ctr(mt_uclr_base, mt_uclr_inc), mt_uclr_base);
        val<1> mt_uclr_sat =
            (mt_uclr_updated == hard<decltype(mt_uclr_ctr)::maxval>{}) |
            (mt_uclr_updated == hard<decltype(mt_uclr_ctr)::minval>{});
        mt_uclr_sat.fanout(hard<MT_NT + 2>{});

        // Allocate on provider_wrong; no-provider starts from lowest table.
        static_assert(MT_NT == 4);
        arr<val<MT_NT>,LINEINST> mt_pick_bits = arr<val<MT_NT>,LINEINST>{[&](u64 offset) -> val<MT_NT> {
            arr<val<1>,MT_NT> mt_allow_vec = arr<val<1>,MT_NT>{[&](u64 t) -> val<1> {
                return select(mt_provider_hit[offset], val<1>{val<2>{t} > mt_provider_id[offset]}, val<1>{1});
            }};
            arr<val<1>,MT_NT> mt_uold_vec = arr<val<1>,MT_NT>{[&](u64 t) -> val<1> {
                return val<1>{mt_meta_ru[t][offset]};
            }};
            val<MT_NT> mt_candidate_bits =
                mt_provider_wrong[offset].replicate(hard<MT_NT>{}).concat() &
                mt_allow_vec.concat() & ~mt_uold_vec.concat();
            return mt_candidate_bits.reverse().one_hot().reverse();
        }};
        mt_pick_bits.fanout(hard<MT_NT>{});
#ifdef PERF_COUNTERS
        perf_count_mt_update_core(
            mt_do_update,
            mt_provider_used,
            mt_provider_wrong,
            mt_pick_bits,
            mt_ctr_gate,
            mt_gshare_wrong,
            mt_uclr_sat);
#endif
        static_loop<MT_NT>([&]<u64 I>(){
            arr<val<1>,LINEINST> mt_sel_i = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_do_update[offset] & mt_provider_hit[offset] & (mt_provider_id[offset] == val<2>{I});
            }};
            arr<val<1>,LINEINST> mt_pick_i = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return (mt_pick_bits[offset] >> hard<I>{}) & val<1>{1};
            }};
            mt_sel_i.fanout(hard<6>{});
            mt_pick_i.fanout(hard<6>{});
            arr<val<1>,LINEINST> mt_alloc_req = mt_pick_i;
            arr<val<1>,LINEINST> mt_pred_update_req = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_sel_i[offset] & mt_provider_wrong[offset] & mt_provider_weak[offset];
            }};
            arr<val<1>,LINEINST> mt_hyst_update_req = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_sel_i[offset] | mt_alloc_req[offset];
            }};
            arr<val<1>,LINEINST> mt_update_u_req = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_sel_i[offset] & mt_update_u[offset];
            }};
            arr<val<1>,LINEINST> mt_uclear_req = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_sel_i[offset] & mt_uclr_sat;
            }};

            // Per-table single-write arbitration across all branch offsets.
            arr<val<1>,LINEINST> mt_tag_req = mt_alloc_req;
            val<LINEINST> mt_tag_pick_mask = mt_tag_req.concat().one_hot();
            arr<val<1>,LINEINST> mt_tag_pick = mt_tag_pick_mask.make_array(val<1>{});
            val<1> mt_tag_we = mt_tag_pick_mask != hard<0>{};
            val<MT_LOGSETS> mt_tag_idx = arr<val<MT_LOGSETS>,LINEINST>{[&](u64 offset) -> val<MT_LOGSETS> {
                return select(mt_tag_pick[offset], val<MT_LOGSETS>{mt_meta_idx[I][offset]}, val<MT_LOGSETS>{0});
            }}.fold_or();
            val<MT_TAGW0> mt_tag_data = arr<val<MT_TAGW0>,LINEINST>{[&](u64 offset) -> val<MT_TAGW0> {
                return select(mt_tag_pick[offset], val<MT_TAGW0>{mt_meta_tag[I][offset]}, val<MT_TAGW0>{0});
            }}.fold_or();
            execute_if(mt_tag_we, [&](){
                mt_tag[I].write(mt_tag_idx, mt_tag_data);
            });

            arr<val<1>,LINEINST> mt_pred_req = mt_pred_update_req;
            val<LINEINST> mt_pred_pick_mask = mt_pred_req.concat().one_hot();
            arr<val<1>,LINEINST> mt_pred_pick = mt_pred_pick_mask.make_array(val<1>{});
            val<1> mt_pred_we = mt_pred_pick_mask != hard<0>{};
            val<MT_LOGSETS> mt_pred_idx = arr<val<MT_LOGSETS>,LINEINST>{[&](u64 offset) -> val<MT_LOGSETS> {
                return select(mt_pred_pick[offset], val<MT_LOGSETS>{mt_meta_idx[I][offset]}, val<MT_LOGSETS>{0});
            }}.fold_or();
            val<1> mt_pred_data = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return select(mt_pred_pick[offset], mt_target_dir[offset], val<1>{0});
            }}.fold_or();
            execute_if(mt_pred_we, [&](){
                mt_pred[I].write(mt_pred_idx, mt_pred_data);
            });

            arr<val<1>,LINEINST> mt_hyst_req = mt_hyst_update_req;
            val<LINEINST> mt_hyst_pick_mask = mt_hyst_req.concat().one_hot();
            arr<val<1>,LINEINST> mt_hyst_pick = mt_hyst_pick_mask.make_array(val<1>{});
            val<1> mt_hyst_we = mt_hyst_pick_mask != hard<0>{};
            val<MT_LOGSETS> mt_hyst_idx = arr<val<MT_LOGSETS>,LINEINST>{[&](u64 offset) -> val<MT_LOGSETS> {
                return select(mt_hyst_pick[offset], val<MT_LOGSETS>{mt_meta_idx[I][offset]}, val<MT_LOGSETS>{0});
            }}.fold_or();
            val<2> mt_hyst_data = arr<val<2>,LINEINST>{[&](u64 offset) -> val<2> {
                val<2> old_h = val<2>{mt_meta_rhyst[I][offset]};
                val<1> provider_correct = ~(mt_provider_pred[offset] != mt_target_dir[offset]);
                val<2> new_h = update_ctr(old_h, provider_correct);
                val<2> write_h = select(mt_alloc_req[offset], val<2>{0}, new_h);
                return select(mt_hyst_pick[offset], write_h, val<2>{0});
            }}.fold_or();
            execute_if(mt_hyst_we, [&](){
                mt_hyst[I].write(mt_hyst_idx, mt_hyst_data, val<1>{1}, extra_cycle);
            });

            arr<val<1>,LINEINST> mt_u_req = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                return mt_update_u_req[offset] | mt_alloc_req[offset] | mt_uclear_req[offset];
            }};
            val<LINEINST> mt_u_pick_mask = mt_u_req.concat().one_hot();
            arr<val<1>,LINEINST> mt_u_pick = mt_u_pick_mask.make_array(val<1>{});
            val<1> mt_u_we = mt_u_pick_mask != hard<0>{};
            val<MT_LOGSETS> mt_u_idx = arr<val<MT_LOGSETS>,LINEINST>{[&](u64 offset) -> val<MT_LOGSETS> {
                return select(mt_u_pick[offset], val<MT_LOGSETS>{mt_meta_idx[I][offset]}, val<MT_LOGSETS>{0});
            }}.fold_or();
            val<1> mt_u_data = arr<val<1>,LINEINST>{[&](u64 offset) -> val<1> {
                val<1> clear_u = mt_alloc_req[offset] | mt_uclear_req[offset];
                val<1> write_u = select(clear_u, val<1>{0}, mt_new_u[offset]);
                return select(mt_u_pick[offset], write_u, val<1>{0});
            }}.fold_or();
#ifdef PERF_COUNTERS
            perf_count_mt_table_update<I>(
                mt_alloc_req,
                mt_pred_update_req,
                mt_hyst_update_req,
                mt_u_req,
                mt_uclear_req,
                mt_tag_we,
                mt_pred_we,
                mt_hyst_we,
                mt_u_we);
#endif
            execute_if(mt_u_we, [&](){
                mt_u[I].write(mt_u_idx, mt_u_data, val<1>{1}, extra_cycle);
            });
        });

        mt_uclr_ctr = select(mt_uclr_sat, val<8>{1 << 7}, mt_uclr_updated);
        execute_if(mt_uclr_sat, [&](){
            static_loop<MT_NT>([&]<u64 I>(){
                mt_u[I].reset();
            });
        });
#endif

        // update incorrect bimodal prediction if primary provider and hysteresis is weak
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(b_weak[offset], [&](){
                bim[offset].write(bindex,branch_taken[offset]);
            });
        }

        // update bimodal hysteresis if bimodal is primary provider
        arr<val<1>,LINEINST> bhyst_write_req = [&](u64 offset) -> val<1> {
            return is_branch[offset] & bhyst_bim_primary[offset] & (~b_weak[offset]);
        };
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(bhyst_write_req[offset], [&](){
                bhyst[offset].write(bindex,~primary_wrong[offset]);
            });
        }
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        {
            const u64 bim_idx = static_cast<u64>(bindex);
            for (u64 offset=0; offset<LINEINST; offset++) {
                if (!static_cast<bool>(bhyst_write_req[offset])) continue;
                perf_shadow_bhyst[offset][bim_idx] = static_cast<u8>(static_cast<bool>(~primary_wrong[offset]));
            }
        }
#endif
#endif

        // update incorrect global prediction if primary provider and the hysteresis is weak;
        // initialize global prediction in the allocated entry
        arr<val<1>,NUMG_GROUP> pred_we = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            u64 i1 = i0 + 1;
            return g_weak[i0] | allocate[i0] | g_weak[i1] | allocate[i1];
        }};
        arr<val<1>,NUMG_GROUP> pred_row0 = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            val<1> we0 = g_weak[i0] | allocate[i0];
            return select(we0, bdir[i0], val<1>{snap_pred[g][0]});
        }};
        arr<val<1>,NUMG_GROUP> pred_row1 = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i1 = 2 * g + 1;
            val<1> we1 = g_weak[i1] | allocate[i1];
            return select(we1, bdir[i1], val<1>{snap_pred[g][1]});
        }};
        // pred_we.fanout(hard<2>{});
        // pred_row0.fanout(hard<2>{});
        // pred_row1.fanout(hard<2>{});
        write_pred_physical(pred_we, pred_row0, pred_row1);

        // update global prediction hysteresis if primary provider or allocated entry
        arr<val<1>,NUMG> g_hyst_write_lane = [&](u64 i) -> val<1> {
            return (primary[i] & select(g_weak[i], val<1>{0}, val<1>{1})) | allocate[i];
        };
        arr<val<1>,NUMG_GROUP> hyst_we = arr<val<1>,NUMG_GROUP>{[&](u64 g) -> val<1> {
            u64 i0 = 2 * g;
            u64 i1 = i0 + 1;
            return (primary[i0] & (~g_weak[i0]) &
                    (~((g_sat[i0])&(~badpred1[i0])))) |
                   allocate[i0] |
                   (primary[i1] & (~g_weak[i1]) &
                    (~((g_sat[i1])&(~badpred1[i1])))) |
                   allocate[i1];
        }};
        arr<val<2>,NUMG_GROUP> hyst_row0 = arr<val<2>,NUMG_GROUP>{[&](u64 g) -> val<2> {
            u64 i0 = 2 * g;
            val<1> we0 = g_hyst_write_lane[i0];
            val<2> newh0 = select(allocate[i0], val<2>{0}, update_ctr(readh[i0], ~badpred1[i0]));
            return select(we0, newh0, val<2>{snap_hyst[g][0]});
        }};
        arr<val<2>,NUMG_GROUP> hyst_row1 = arr<val<2>,NUMG_GROUP>{[&](u64 g) -> val<2> {
            u64 i1 = 2 * g + 1;
            val<1> we1 = g_hyst_write_lane[i1];
            val<2> newh1 = select(allocate[i1], val<2>{0}, update_ctr(readh[i1], ~badpred1[i1]));
            return select(we1, newh1, val<2>{snap_hyst[g][1]});
        }};
        // hyst_we.fanout(hard<2>{});
        // hyst_row0.fanout(hard<2>{});
        // hyst_row1.fanout(hard<2>{});
        write_hyst_physical(hyst_we, hyst_row0, hyst_row1, extra_cycle);

#ifdef RESET_UBITS
        uctr.fanout(hard<3>{});
        val<NUMG> allocmask1  = collamask1.reverse();
        allocmask1.fanout(hard<2>{});
        val<1> faralloc = (((last_match1>>3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
        val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
        uctrsat.fanout(hard<2>{});
        uctr = select(correct_pred,uctr,select(uctrsat,val<decltype(uctr)::size>{0},update_ctr(uctr,faralloc.fo1())));
        execute_if(uctrsat,[&](){for (auto &uram : ubit_p) uram.reset();});
#endif


        // update global history
        val<1> line_end = block_entry >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            next_pc.fanout(hard<2>{});
            global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
            gfolds.update(val<PATHBITS>{next_pc>>2});
#if MY_BP_V1_P1_MICRO_TAGE
            mt_gfolds.update(val<PATHBITS>{next_pc>>2});
#endif
        });

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_count_end_of_cycle(is_branch, branch_taken, mispredict, allocate, postmask, noalloc);
#endif
#endif

        num_branch = 0; // done
    }
};

#ifdef PERF_COUNTERS
#undef perf_predictions
#undef perf_correct
#undef perf_provider_used
#undef perf_provider_correct
#undef perf_provider_wrong
#undef perf_alt_used
#undef perf_alt_correct
#undef perf_alt_wrong
#undef perf_bimodal_used
#undef perf_bimodal_correct
#undef perf_bimodal_wrong
#undef perf_mispred_blame_tage
#undef perf_mispred_blame_p1
#undef perf_table_reads
#undef perf_table_hits
#undef perf_table_alloc
#undef perf_alloc_failures
#undef perf_alloc_fail_highest
#undef perf_alloc_fail_noubit
#undef perf_extra_cycle_total
#undef perf_extra_cycle_badpred
#undef perf_extra_cycle_mispredict
#undef perf_extra_cycle_p1_update
#undef perf_phys_tag_w
#undef perf_phys_pred_w
#undef perf_phys_hyst_w
#undef perf_phys_u_w
#undef perf_hyst_skip_p1
#undef perf_hyst_skip_bim
#undef perf_hyst_skip_tage
#undef perf_weak_wrong_p1
#undef perf_weak_wrong_bim
#undef perf_weak_wrong_tage
#undef perf_conf
#endif
