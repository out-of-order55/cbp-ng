// this is a basic TAGE, not necessarily well optimized

#define USE_ALT

#ifndef MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
#define MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE 1
#endif

#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE != 0 && MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE != 1
#error "MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE must be 0 or 1"
#endif

#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
#define MY_BP_V1_AHEAD_PIPE_DEFAULT_NUMG 14
#else
#define MY_BP_V1_AHEAD_PIPE_DEFAULT_NUMG 8
#endif

#define RESET_UBITS

#ifndef MY_BP_V1_PERF_ENABLE_TRACE_DB
#define MY_BP_V1_PERF_ENABLE_TRACE_DB 1
#endif

#ifndef MY_BP_V1_PERF_PRINT_TRACE
#define MY_BP_V1_PERF_PRINT_TRACE 1
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
struct split_row_ram2_ahead {
    static_assert(std::has_single_bit(M));
    static constexpr u64 A = std::bit_width(M - 1);
    ram<val<W>, M> lane[2];

    split_row_ram2_ahead(const char *label = "") : lane{label} {}

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

template<u64 LOGLB=6, u64 NUMG=MY_BP_V1_AHEAD_PIPE_DEFAULT_NUMG, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=300, u64 LOGP1=13, u64 GHIST1=6,u64 LOGBANKS = 1,u64 LOGBIAS = 11>
struct my_bp_v1_ahead_pipe : predictor {
    // provides 2^(LOGLB-2) predictions per cycle
    // P1/P2 share ahead-pipe TAGE; P1 falls back to bimodal when pipe data is not ready.
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
    static_assert((NUMG % 2) == 0);
    static constexpr u64 NUMP = NUMG / 2;
    static constexpr u64 NUMG_GROUP = NUMP;
    static constexpr u64 LOGG_SH = LOGG;
    static constexpr u64 P2_LINE_BITS = (bindex_bits > LOGG) ? bindex_bits : LOGG;



    geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;

    reg<1> true_block = 1;

    // for P1
    reg<LINEINST> p1; // P1 predictions






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
    reg<1> p2_pipe_valid = 0;
    reg<1> p2_consume_valid = 0;
    reg<bindex_bits> p2_pipe_bindex = 0;
    arr<reg<LOGG>,NUMG> p2_pipe_gindex;
    arr<reg<HTAGBITS>,NUMG> p2_pipe_htag;
    arr<reg<TAGW>,NUMG> p2_pipe_readt;
    arr<reg<1>,NUMG> p2_pipe_readc;
    arr<reg<2>,NUMG> p2_pipe_readh;
    arr<reg<1>,NUMG> p2_pipe_readu;

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
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
        for (u64 j = 0; j < NUMP; j++) {
            ghyst_reads += ghyst_p[j].perf_total_reads();
            ghyst_writes += ghyst_p[j].perf_total_writes();
            ghyst_stale_reads += ghyst_p[j].perf_total_stale_reads();
            ghyst_drop_writes += ghyst_p[j].perf_total_drop_writes();
            ghyst_write_hit_events += ghyst_p[j].perf_total_write_updates();
        }
#else
        for (u64 j = 0; j < NUMG; j++) {
            ghyst_reads += ghyst[j].perf_total_reads();
            ghyst_writes += ghyst[j].perf_total_writes();
            ghyst_stale_reads += ghyst[j].perf_total_stale_reads();
            ghyst_drop_writes += ghyst[j].perf_total_drop_writes();
            ghyst_write_hit_events += ghyst[j].perf_total_write_updates();
        }
#endif

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
        }
    }

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
        arr<val<1>,LINEINST> &primary_wrong,
        arr<val<1>,LINEINST> &bhyst_bim_primary,
        arr<val<1>,NUMG> &primary,
        arr<val<1>,NUMG> &allocate,
        arr<val<1>,NUMG> &badpred1,
        arr<val<1>,NUMG> &g_sat)
    {
        const u64 bim_idx = static_cast<u64>(bindex);

        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

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
        arr<val<1>,LINEINST> b_weak,
        arr<val<1>,NUMG> g_weak)
    {
        for (u64 offset = 0; offset < LINEINST; offset++) {
            if (!static_cast<bool>(is_branch[offset])) continue;

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

    // Keep prediction-path tables in one zone to minimize cross-zone wiring.
    ram<val<1>,(1<<bindex_bits)> bim[LINEINST] {"bpred"}; // bimodal prediction bits
    // P2 (TAGE)
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
    split_row_ram2_ahead<TAGW,(1<<LOGG_SH)> gtag_p[NUMP] {"tags"}; // tags (2-slot row, split 2xRAM)
    ram<arr<val<1>,2>,(1<<LOGG_SH)> gpred_p[NUMP] {"gpred"}; // predictions (2-slot row)
    rwram<4,(1<<LOGG_SH),4> ghyst_p[NUMP] {"ghyst"}; // packed hysteresis row: {slot1,slot0}
    rwram<2,(1<<LOGG_SH),4> ubit_p[NUMP] {"u"}; // packed u-bit row: {slot1,slot0}
#else
    ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"tags"};
    ram<val<1>,(1<<LOGG)> gpred[NUMG] {"gpred"};
    rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"};
    rwram<1,(1<<LOGG),4> ubit[NUMG] {"u"};
#endif
    zone UPDATE_ONLY;
    ram<val<1>,(1<<bindex_bits)> bhyst[LINEINST] {"bhyst"}; // bimodal hysteresis

// #endif
    my_bp_v1_ahead_pipe()
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


    ~my_bp_v1_ahead_pipe() {
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

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<3>{});
        new_block(inst_pc);

        // Move ahead-pipe TAGE to P1 path. P1/P2 consume the same prediction vector.
        tage_pred_ahead_pipe(inst_pc);
        p1 = tage_p2;
        p1.fanout(hard<2*LINEINST>{});
        reuse_prediction(~val<1>{block_entry>>(LINEINST-1)});
        val<1> p1_taken = (block_entry & p1) != hard<0>{};
        return p1_taken;
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        val<1> p1_taken = ((block_entry<<block_size) & p1) != hard<0>{};
        reuse_prediction(~val<1>{block_entry>>(LINEINST-1-block_size)});
        return p1_taken;
    };
    void read_tag_pred_physical()
    {
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
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
#else
        static_loop<NUMG>([&]<u64 I>(){
            val<LOGG> row = val<LOGG>{gindex[I]};
            row.fanout(hard<4>{});
            readt[I] = gtag[I].read(row);
            readc[I] = gpred[I].read(row);
            readh[I] = ghyst[I].read(row);
            readu[I] = ubit[I].read(row);
        });

        static_loop<NUMG_GROUP>([&]<u64 G>(){
            snap_tag[G][0] = readt[2*G];
            snap_tag[G][1] = readt[2*G+1];
            snap_pred[G][0] = readc[2*G];
            snap_pred[G][1] = readc[2*G+1];
            snap_hyst[G][0] = readh[2*G];
            snap_hyst[G][1] = readh[2*G+1];
            snap_u[G][0] = readu[2*G];
            snap_u[G][1] = readu[2*G+1];
        });
#endif
    }

    void write_tag_physical(
        arr<val<1>,NUMG_GROUP> &tag_we,
        arr<val<TAGW>,NUMG_GROUP> &tag_row0,
        arr<val<TAGW>,NUMG_GROUP> &tag_row1)
    {
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
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
#else
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = tag_we[G];
            constexpr u64 I0 = 2 * G;
            constexpr u64 I1 = I0 + 1;
            val<LOGG> row0 = val<LOGG>{gindex[I0]};
            val<LOGG> row1 = val<LOGG>{gindex[I1]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_tag_w, G, we);
#endif
#endif
            execute_if(we, [&](){
                gtag[I0].write(row0, tag_row0[G]);
                gtag[I1].write(row1, tag_row1[G]);
            });
        });
#endif
    }

    void write_pred_physical(
        arr<val<1>,NUMG_GROUP> &pred_we,
        arr<val<1>,NUMG_GROUP> &pred_row0,
        arr<val<1>,NUMG_GROUP> &pred_row1)
    {
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
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
#else
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = pred_we[G];
            constexpr u64 I0 = 2 * G;
            constexpr u64 I1 = I0 + 1;
            val<LOGG> row0 = val<LOGG>{gindex[I0]};
            val<LOGG> row1 = val<LOGG>{gindex[I1]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_pred_w, G, we);
#endif
#endif
            execute_if(we, [&](){
                gpred[I0].write(row0, pred_row0[G]);
                gpred[I1].write(row1, pred_row1[G]);
            });
        });
#endif
    }

    void write_hyst_physical(
        arr<val<1>,NUMG_GROUP> &hyst_we,
        arr<val<2>,NUMG_GROUP> &hyst_row0,
        arr<val<2>,NUMG_GROUP> &hyst_row1,
        val<1> noconflict)
    {
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
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
#else
        noconflict.fanout(hard<NUMG>{});
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = hyst_we[G];
            constexpr u64 I0 = 2 * G;
            constexpr u64 I1 = I0 + 1;
            val<LOGG> row0 = val<LOGG>{gindex[I0]};
            val<LOGG> row1 = val<LOGG>{gindex[I1]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_hyst_w, G, we);
#endif
#endif
            execute_if(we, [&](){
                ghyst[I0].write(row0, hyst_row0[G], val<1>{1}, noconflict);
                ghyst[I1].write(row1, hyst_row1[G], val<1>{1}, noconflict);
            });
        });
#endif
    }

    void write_ubit_physical(
        arr<val<1>,NUMG_GROUP> &u_we,
        arr<val<1>,NUMG_GROUP> &u_row0,
        arr<val<1>,NUMG_GROUP> &u_row1,
        val<1> noconflict)
    {
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
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
#else
        noconflict.fanout(hard<NUMG>{});
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<1> we = u_we[G];
            constexpr u64 I0 = 2 * G;
            constexpr u64 I1 = I0 + 1;
            val<LOGG> row0 = val<LOGG>{gindex[I0]};
            val<LOGG> row1 = val<LOGG>{gindex[I1]};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
            perf_count_phys_write(perf_phys_u_w, G, we);
#endif
#endif
            execute_if(we, [&](){
                ubit[I0].write(row0, u_row0[G], val<1>{1}, noconflict);
                ubit[I1].write(row1, u_row1[G], val<1>{1}, noconflict);
            });
        });
#endif
    }

    void tage_finalize_prediction()
    {
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
        for (u64 offset=0; offset<LINEINST; offset++) {
            match1[offset] = match[offset].one_hot();
        }
        match1.fanout(hard<3>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
        }
        pred1.fanout(hard<3>{});

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

    void tage_pred_ahead_pipe(val<64> inst_pc)
    {
        // Consume stage: only use previous-cycle pipe registers.
        bindex = val<bindex_bits>{p2_pipe_bindex};
        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = val<LOGG>{p2_pipe_gindex[i]};
            htag[i] = val<HTAGBITS>{p2_pipe_htag[i]};
            readt[i] = val<TAGW>{p2_pipe_readt[i]};
            readc[i] = val<1>{p2_pipe_readc[i]};
            readh[i] = val<2>{p2_pipe_readh[i]};
            readu[i] = val<1>{p2_pipe_readu[i]};
        }
        for (u64 offset=0; offset<LINEINST; offset++) {
            readb[offset] = bim[offset].read(bindex);
        }
        for (u64 g=0; g<NUMG_GROUP; g++) {
            group_row_idx[g] = val<LOGG_SH>{gindex[2*g]};
            snap_tag[g][0] = readt[2*g];
            snap_tag[g][1] = readt[2*g+1];
            snap_pred[g][0] = readc[2*g];
            snap_pred[g][1] = readc[2*g+1];
            snap_hyst[g][0] = readh[2*g];
            snap_hyst[g][1] = readh[2*g+1];
            snap_u[g][0] = readu[2*g];
            snap_u[g][1] = readu[2*g+1];
        }
        tage_finalize_prediction();
        val<1> pipe_ready = val<1>{p2_pipe_valid};
        p2_consume_valid = pipe_ready;

        // Issue stage: read current line and fill pipe registers for next cycle.
        val<P2_LINE_BITS> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<8 + NUMG * 4>{});
        gfolds.fanout(hard<2>{});
        val<bindex_bits> bindex_now = lineaddr;
        bindex_now.fanout(hard<LINEINST + 2>{});
        p2_pipe_bindex = bindex_now;
        for (u64 i=0; i<NUMG; i++) {
            val<LOGG> idx = lineaddr ^ gfolds.template get<0>(i);
            p2_pipe_gindex[i] = idx;
            p2_pipe_htag[i] = val<HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
        }
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
        static_loop<NUMG_GROUP>([&]<u64 G>(){
            val<LOGG_SH> row = val<LOGG_SH>{lineaddr ^ gfolds.template get<0>(2*G)};
            row.fanout(hard<4>{});
            arr<val<TAGW>,2> tag = gtag_p[G].read(row);
            arr<val<1>,2> pred = gpred_p[G].read(row);
            val<4> hyst_raw = ghyst_p[G].read(row);
            val<2> u_raw = ubit_p[G].read(row);
            p2_pipe_readt[2*G] = tag[0];
            p2_pipe_readt[2*G+1] = tag[1];
            p2_pipe_readc[2*G] = pred[0];
            p2_pipe_readc[2*G+1] = pred[1];
            p2_pipe_readh[2*G] = val<2>{hyst_raw};
            p2_pipe_readh[2*G+1] = val<2>{hyst_raw >> 2};
            p2_pipe_readu[2*G] = val<1>{u_raw};
            p2_pipe_readu[2*G+1] = val<1>{u_raw >> 1};
        });
#else
        static_loop<NUMG>([&]<u64 I>(){
            val<LOGG> row = val<LOGG>{p2_pipe_gindex[I]};
            row.fanout(hard<4>{});
            p2_pipe_readt[I] = gtag[I].read(row);
            p2_pipe_readc[I] = gpred[I].read(row);
            p2_pipe_readh[I] = ghyst[I].read(row);
            p2_pipe_readu[I] = ubit[I].read(row);
        });
#endif
        p2_pipe_valid = val<1>{1};
    }
    val<1> predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = (block_entry & p1) != hard<0>{};
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = ((block_entry<<block_size) & p1) != hard<0>{};
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
                gfolds.update(val<PATHBITS>{next_pc>>2});
                true_block = 1;
            });
            return; // stop here
        }
        mispredict.fanout(hard<3>{});
        val<1> correct_pred = ~mispredict;
        correct_pred.fanout(hard<NUMG+2>{});
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

        val<1> extra_cycle_base = some_badpred1 | mispredict;
        val<1> extra_cycle = extra_cycle_base;

        extra_cycle.fanout(hard<NUMG*2+8>{});
        need_extra_cycle(extra_cycle);
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        val<1> p1_update = hard<0>{};
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

#ifdef PERF_COUNTERS
        perf_count_weak_wrong_events(is_branch, b_weak, g_weak);
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
#if MY_BP_V1_AHEAD_PIPE_ADJACENTTABLE
        execute_if(uctrsat,[&](){for (auto &uram : ubit_p) uram.reset();});
#else
        execute_if(uctrsat,[&](){for (auto &uram : ubit) uram.reset();});
#endif
#endif


        // update global history
        val<1> line_end = block_entry >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            gfolds.update(val<PATHBITS>{next_pc>>2});
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
