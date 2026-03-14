// this is a basic TAGE, not necessarily well optimized

#define USE_META
#define RESET_UBITS

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

using namespace hcm;

template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=400, u64 LOGP1=14, u64 GHIST1=6>
struct tage : predictor {
    // provides 2^(LOGLB-2) predictions per cycle
    // P2 is a TAGE, P1 is a gshare
    static_assert(LOGLB>2);
    static_assert(NUMG>0);
    static constexpr u64 MINHIST = 2;
    static constexpr u64 METABITS = 4;
    static constexpr u64 UCTRBITS = 8;
    static constexpr u64 PATHBITS = 6;
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
    static constexpr u64 HTAGBITS = TAGW-LOGLINEINST; // hashed tag bits

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

    arr<reg<1>,LINEINST> pred1; // primary P2 prediction for each offset
    arr<reg<1>,LINEINST> pred2; // alternate P2 prediction for each offset
    reg<LINEINST> p2; // final P2 predictions

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

    // Branch execution flow trace (one entry per branch execution)
    struct ExecRecord {
        u64 cycle = 0;
        u64 pc = 0;
        u64 offset = 0;
        u64 actual_dir = 0;
        u64 predicted_dir = 0;
        u64 mispredict = 0;
        // prediction source: 0=bimodal, 1=tage_provider, 2=tage_alt
        u64 pred_source = 0;
        u64 pred_table = 0; // which TAGE table (NUMG=bimodal)
        u64 bim_index = 0;
        u64 hit = 0;
        u64 hit_table = 0;
        u64 hit_gtag = 0;
        u64 hit_gindex = 0;
        u64 alloc = 0;
        u64 alloc_table = 0;
        u64 alloc_gindex = 0;
        u64 alloc_tag = 0;
    };
    std::vector<ExecRecord> exec_trace;

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

        // Per-table statistics
        std::cerr << "\n┌─ TAGE Table Statistics ─────────────────────────────────────────────────────────────────────────────────┐\n";
        std::cerr << "│ Tbl │ HistLen │ Reads │ Hits │ Hit%  │ Prov │ PrvAcc% │ Alt │ AltAcc% │ Total │ TotAcc% │ Alloc │\n";
        std::cerr << "├─────┼─────────┼───────┼──────┼───────┼──────┼─────────┼─────┼─────────┼───────┼─────────┼───────┤\n";

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

            std::cerr << "│ " << std::setw(3) << std::left << j << " │ ";
            std::cerr << std::setw(7) << std::right << gfolds.HLEN[j] << " │ ";
            std::cerr << std::setw(5) << std::right << reads << " │ ";
            std::cerr << std::setw(4) << std::right << hits << " │ ";

            if (reads > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(5) << std::right
                          << (100.0*hits/reads) << "% │ ";
            else
                std::cerr << "  N/A │ ";

            std::cerr << std::setw(4) << std::right << prov << " │ ";
            if (prov > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right
                          << (100.0*prov_ok/prov) << "% │ ";
            else
                std::cerr << "    N/A │ ";

            std::cerr << std::setw(3) << std::right << alt << " │ ";
            if (alt > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right
                          << (100.0*alt_ok/alt) << "% │ ";
            else
                std::cerr << "    N/A │ ";

            std::cerr << std::setw(5) << std::right << total << " │ ";
            if (total > 0)
                std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right
                          << (100.0*total_ok/total) << "% │ ";
            else
                std::cerr << "    N/A │ ";

            std::cerr << std::setw(5) << std::right << alloc << " │\n";
        }
        std::cerr << "└─────┴─────────┴───────┴──────┴───────┴──────┴─────────┴─────┴─────────┴───────┴─────────┴───────┘\n";

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

        // Allocation failures
        std::cerr << "\n┌─ Allocation Statistics ─────────────────────────────────────────┐\n";
        std::cerr << "│ Allocation Failures: " << std::setw(41) << std::left
                  << perf_alloc_failures << "│\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Verification
        std::cerr << "\n┌─ Verification ──────────────────────────────────────────────────┐\n";
        std::cerr << "│ Source total matches predictions: ";
        if (total_all == perf_predictions)
            std::cerr << "✓ PASS                      │\n";
        else
            std::cerr << "✗ FAIL (" << total_all << " vs " << perf_predictions << ")    │\n";
        std::cerr << "└─────────────────────────────────────────────────────────────────┘\n";

        // Write full exec trace to file
        {
            std::ofstream ef("mispred_trace_tage.csv");
            ef << "seq,cycle,pc,offset,actual_dir,predicted_dir,mispredict,pred_source,pred_table,bim_index,hit,hit_table,hit_gtag,hit_gindex,alloc,alloc_table,alloc_gindex,alloc_tag\n";
            for (u64 i = 0; i < exec_trace.size(); i++) {
                auto &r = exec_trace[i];
                std::string src_str = (r.pred_source == 0) ? "bimodal" :
                                      (r.pred_source == 1) ? "tage_prov" : "tage_alt";
                ef << i << ","
                   << r.cycle << ","
                   << "0x" << std::hex << r.pc << std::dec << ","
                   << r.offset << ","
                   << r.actual_dir << ","
                   << r.predicted_dir << ","
                   << r.mispredict << ","
                   << src_str << ","
                   << (r.pred_source != 0 ? std::to_string(r.pred_table) : "N/A") << ","
                   << r.bim_index << ","
                   << r.hit << ","
                   << (r.hit ? std::to_string(r.hit_table) : "N/A") << ","
                   << (r.hit ? ("0x" + [&]{ std::ostringstream s; s << std::hex << r.hit_gtag; return s.str(); }()) : "N/A") << ","
                   << (r.hit ? std::to_string(r.hit_gindex) : "N/A") << ","
                   << r.alloc << ","
                   << (r.alloc ? std::to_string(r.alloc_table) : "N/A") << ","
                   << (r.alloc ? std::to_string(r.alloc_gindex) : "N/A") << ","
                   << (r.alloc ? ("0x" + [&]{ std::ostringstream s; s << std::hex << r.alloc_tag; return s.str(); }()) : "N/A") << "\n";
            }
        }

        // Write per-PC mispred summary + top-20 to stderr
        {
            std::vector<std::pair<u64,MispredRecord>> sorted_db(mispred_db.begin(), mispred_db.end());
            std::sort(sorted_db.begin(), sorted_db.end(),
                [](const auto &a, const auto &b){ return a.second.count > b.second.count; });

            std::cerr << "\n┌─ Top 20 Most Mispredicted PCs (full trace -> mispred_trace_tage.csv) ──────────────────────┐\n";
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

    // P1 (gshare)
    ram<val<1>,(1<<index1_bits)> table1_pred[LINEINST] {"P1 pred"}; // P1 prediction bit

    // P2 (TAGE)
    ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"tags"}; // tags
    ram<val<1>,(1<<LOGG)> gpred[NUMG] {"gpred"}; // predictions
    rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"}; // hysteresis
    rwram<1,(1<<LOGG),4> ubit[NUMG] {"u"}; // "useful" bits
    ram<val<1>,(1<<bindex_bits)> bim[LINEINST] {"bpred"}; // bimodal prediction bits

    zone UPDATE_ONLY;
    ram<val<1>,(1<<index1_bits)> table1_hyst[LINEINST] {"P1 hyst"}; // P1 hysteresis
    ram<val<1>,(1<<bindex_bits)> bhyst[LINEINST] {"bhyst"}; // bimodal hysteresis

    tage()
    {
#ifdef TAGE_VERBOSE
        std::cerr << "TAGE history lengths: ";
        for (u64 i=0; i<NUMG; i++) std::cerr << gfolds.HLEN[i] << " ";
        std::cerr << std::endl;
        if (LOGG == HTAGBITS) {
            std::cerr << "WARNING: the tag function and index function are not different enough\n";
        }
#endif
    }

#ifdef PERF_COUNTERS
    ~tage() {
        print_perf_counters();
    }
#endif

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

    val<1> predict2(val<64> inst_pc)
    {
        val<std::max(bindex_bits,LOGG)> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<1+NUMG*2>{});
        gfolds.fanout(hard<2>{});

        // compute indexes
        bindex = lineaddr;
        bindex.fanout(hard<LINEINST>{});
        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = lineaddr ^ gfolds.template get<0>(i);
        }
        gindex.fanout(hard<4>{});

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
        for (u64 i=0; i<NUMG; i++) {
            readt[i] = gtag[i].read(gindex[i]);
            readc[i] = gpred[i].read(gindex[i]);
            readh[i] = ghyst[i].read(gindex[i]);
            readu[i] = ubit[i].read(gindex[i]);
        }
        readt.fanout(hard<LINEINST+1>{});
        readc.fanout(hard<3>{});
        readh.fanout(hard<2>{});
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
        p2 = arr<val<1>,LINEINST> {[&](u64 offset){
            return select(altsel[offset].fo1(),pred2[offset],pred1[offset]);
        }}.concat();
#else
        p2 = pred1.concat();
#endif
        p2.fanout(hard<LINEINST>{});
        val<1> taken = (block_entry & p2) != hard<0>{};
        taken.fanout(hard<2>{});
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
        is_branch.fanout(hard<6>{});

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
            // bits 0..NUMG-1 = tables, bit NUMG = bimodal
            val<NUMG> all_hits = val<NUMG+1>{pred_match1_stored[offset]} & hard<(1<<NUMG)-1>{};
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
        branch_taken.fanout(hard<3>{});

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
    // Track allocations per table
    for (u64 j=0; j<NUMG; j++) {
        if (static_cast<bool>(allocate[j])) {
            perf_table_alloc[j]++;
        }
    }
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
        val<1> extra_cycle = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
        extra_cycle.fanout(hard<NUMG*2+1>{});
        need_extra_cycle(extra_cycle);

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
        val<NUMG> uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<2>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
    if (static_cast<bool>(noalloc) && static_cast<bool>(mispredict)) {
        perf_alloc_failures++;
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

        // Record exec trace entry
        ExecRecord er;
        er.cycle        = static_cast<u64>(panel.cycle);
        er.pc           = pc_val;
        er.offset       = offset;
        er.actual_dir   = static_cast<u64>(actual);
        er.predicted_dir = static_cast<u64>(predicted);
        er.mispredict   = is_misp & is_last;
        er.pred_source  = pred_source;
        er.pred_table   = pred_table;
        er.bim_index    = static_cast<u64>(bindex);
        er.hit          = hit_found;
        er.hit_table    = hit_table;
        er.hit_gtag     = hit_gtag;
        er.hit_gindex   = hit_gindex;
        er.alloc        = is_misp & is_last & alloc_found;
        er.alloc_table  = alloc_table;
        er.alloc_gindex = alloc_gindex_val;
        er.alloc_tag    = alloc_tag_val;
        exec_trace.push_back(er);

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
#endif
#endif

        num_branch = 0; // done
    }
};
