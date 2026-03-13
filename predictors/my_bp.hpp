// this is a basic TAGE, not necessarily well optimized

#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#define USE_ALT_PRED
#define RESET_UBITS
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
template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGW=12, u64 GHIST=100, u64 LOGP1=14, u64 GHIST1=6
, u64 NUMBANKS=4,u64 NUMWAYS=1,u64 CTRBIT=3,u64 UBIT=2>
struct my_bp : predictor {
    static_assert(LOGLB>2);
    static_assert(NUMG>0);
    static constexpr u64 MINHIST = 2;
    static constexpr u64 USE_ALT_PRED_BITS = 4;
    static constexpr u64 UCTRBITS = 8;
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

    // HCpred: high-confidence prediction (strongest non-weak entry)
    arr<reg<1>,LINEINST> hc_pred_val[NUMWAYS];      // HC predicted direction
    arr<reg<NUMG>,LINEINST> hc_pred_source[NUMWAYS]; // HC source bitmask (0=bimodal, one-hot=TAGE table)
    arr<reg<2>,LINEINST> hc_pred_type[NUMWAYS]; // HC source type: 0=bim, 1=prov, 2=alt


    arr<reg<1>,LINEINST> use_provider[NUMWAYS];
    reg<USE_ALT_PRED_BITS> use_alt_on_na;


#ifdef RESET_UBITS
    reg<UCTRBITS> uctr;
#endif

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
        printf("Useful Allocations: %lu\n", (unsigned long)perf_useful_alloc);
        printf("Useful Increments: %lu\n", (unsigned long)perf_useful_inc);
        printf("Useful Resets: %lu\n", (unsigned long)perf_useful_reset);
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

        std::cerr << "\n";
    }
#endif

#ifdef CHEATING_MODE
    std::unordered_map<u64, u64> mispred_pc_count;
    u64 total_mispred_branches = 0;

    void record_mispred_branch_pc(u64 pc)
    {
        mispred_pc_count[pc]++;
        total_mispred_branches++;
    }

    void dump_mispred_branch_db_csv(const std::string &path = "results/mispred_branches.csv")
    {
        std::vector<std::pair<u64, u64>> rows;
        rows.reserve(mispred_pc_count.size());
        for (const auto &entry : mispred_pc_count) rows.push_back(entry);

        std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
            if (a.second != b.second) return a.second > b.second;
            return a.first < b.first;
        });

        std::ofstream out(path);
        if (!out.is_open()) {
            std::cerr << "[mispred-db] failed to open " << path << " for writing\n";
            return;
        }

        out << "pc_hex,pc_dec,mispred_count\n";
        for (const auto &row : rows) {
            std::ios old_state(nullptr);
            old_state.copyfmt(out);
            out << "0x" << std::hex << row.first;
            out.copyfmt(old_state);
            out << "," << row.first << "," << row.second << "\n";
        }
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
    rwram<CTRBIT,NUMGTSETS,NUMBANKS> gctr[NUMWAYS][NUMG] {"ctr"};
    //
    rwram<UBIT,NUMGTSETS,NUMBANKS> ubit[NUMWAYS][NUMG] {"uctr"};

    //116ps
    rwram<1,1<<bindex_bits,NUMBANKS> bim_hi[LINEINST] {"bpred"};


    zone UPDATE_ONLY;
    rwram<1,1<<index1_bits,NUMBANKS> table1_hyst[LINEINST] {"P1 hyst"};
    rwram<1,1<<bindex_bits,NUMBANKS> bim_low[LINEINST] {"bhyst"};

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
    }


    void new_block(val<64> inst_pc)
    {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        inst_oh = offset.fo1().decode().concat();
        inst_oh.fanout(hard<6*LINEINST>{});
        block_size = 1;
    }

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        new_block(inst_pc);
        // val<std::max(P1INDEXBITS+BANKBITS,GHIST1)> lineaddr = inst_pc >> (LOGLB);
        // val<BANKBITS> bankid = inst_pc >> (LOGLB);
        val<index1_bits> raw_index = inst_pc >> (LOGLB);

        
        if constexpr (GHIST1 <= index1_bits) {
            index1 = raw_index ^ (val<index1_bits>{path_history}<<(index1_bits-GHIST1));
        } else {
            index1 = path_history.make_array(val<index1_bits>{}).append(raw_index).fold_xor();
        }
        index1.fanout(hard<LINEINST>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            readp1[offset] = table1_pred[offset].read(index1);
        }
        readp1.fanout(hard<2>{});
        p1 = readp1.concat();
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

    val<1> predict2(val<64> inst_pc)
    {
        val<HTAGBITS>   raw_tag = inst_pc >> (LOGLB);
        val<LOGG>       raw_gindex = inst_pc >> (LOGLB);

        val<bindex_bits> raw_bindex = inst_pc >> (LOGLB);

        raw_tag.fanout(hard<NUMG>{});
        raw_gindex.fanout(hard<NUMG>{});
        gfolds.fanout(hard<2>{});
        

        bindex = raw_bindex.fo1();
        bindex.fanout(hard<2*LINEINST>{});

        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = raw_gindex ^ gfolds.template get<0>(i);
        }
        gindex.fanout(hard<3*NUMWAYS>{});

        for (u64 i=0; i<NUMG; i++) {
            htag[i] = raw_tag.reverse() ^ gfolds.template get<1>(i);
        }
        // htag.fanout(hard<2>{});

        for (u64 offset=0; offset<LINEINST; offset++) {
            readb[offset] = bim_hi[offset].read(bindex);
            readb_low[offset] = bim_low[offset].read(bindex);
        }
        // readb.fanout(hard<NUMWAYS>{});

        for (u64 w=0; w<NUMWAYS; w++) {
            for(u64 j=0; j<NUMG; j++){
                //146ps
                arr<val<1>,CTRBIT> ctr_val = gctr[w][j].read(gindex[j]).make_array(val<1>{});
                readt[w][j] = gtag[w][j].read(gindex[j]);
                readu[w][j] = ubit[w][j].read(gindex[j]);
                readctr_pred[w][j] = ctr_val[CTRBIT-1];
                readctr_cnt[w][j] = val<CTRBIT-1>{ctr_val.concat()};
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

        //TODO: fix this
        for(u64 way=0;way<NUMWAYS;way++){
            notumask[way] = ~umask[way].fo1();
        }
        // notumask.fanout(hard<2>{});

        //12ps 2000fJ
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

        //20ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                match_provider[w][offset] = match[w][offset].one_hot();
            }
        }

        arr<val<NUMG+1>,NUMWAYS> gpreds = [&](u64 w){ return readctr_pred[w].concat(); };
        gpreds.fanout(hard<2>{});
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
        for (u64 w=0; w<NUMWAYS; w++) alt[w].fanout(hard<2>{});


        use_alt_on_na.fanout(hard<2>{});
        arr<val<1>,USE_ALT_PRED_BITS> use_alt_on_na_array = use_alt_on_na.make_array(val<1>{});
        val<1> use_alt_on_na_pos = use_alt_on_na_array[USE_ALT_PRED_BITS-1];
        arr<val<1>,LINEINST> use_alt_on_na_pos_dup = use_alt_on_na_pos.fo1().replicate(hard<LINEINST>{});//reduce fanout
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 inst=0; inst<LINEINST; inst++) {
                arr<val<1>,NUMG> prov_oh = val<NUMG>{match_provider[w][inst]}.make_array(val<1>{});
                arr<val<1>,NUMG> alt_oh = val<NUMG>{match_alt[w][inst]}.make_array(val<1>{});
                prov_oh.fanout(hard<3>{});
                arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                    return select(prov_oh[j],readctr_cnt[w][j],val<CTRBIT-1>{0});
                };
                val<CTRBIT-1> prov_ctr  = prov_ctrs.fo1().fold_or();
                prov_ctr.fanout(hard<2>{});
                // provider & alt
                val<1> prov_is_weak = is_weak(provider[w][inst], prov_ctr);
                prov_is_weak.fanout(hard<2>{});
#ifdef USE_ALT_PRED
                use_provider[w][inst] = prov_oh.fold_or() &
                select(alt_oh.fo1().fold_or(),(~prov_is_weak | ~use_alt_on_na_pos_dup[inst]),val<1>{1});
#else
                use_provider[w][inst] = prov_oh.fold_or();
#endif
                // ── HCpred ──────────────────────────────────────────────
                // Provider: not weak and exists → HCpred = provider
                
                val<1> use_provider_hc = prov_oh.fold_or() & ~prov_is_weak;
                use_provider_hc.fanout(hard<3>{});
                // Alt: from all alt hits (not just one_hot), find longest non-weak
                val<NUMG> all_alt_hits = val<NUMG>{match[w][inst]} ^ val<NUMG>{match_provider[w][inst]};
                all_alt_hits.fanout(hard<NUMG>{});
                arr<val<1>,NUMG> hc_alt_cand = [&](u64 j){
                    val<1> hit = all_alt_hits.make_array(val<1>{})[j];
                    val<1> not_weak = ~is_weak(readctr_pred[w][j], readctr_cnt[w][j]);
                    return hit.fo1() & not_weak.fo1();
                };
                hc_alt_cand.fanout(hard<2>{});
                // Select the longest (highest-indexed) non-weak alt
                val<NUMG> hc_alt_mask = hc_alt_cand.concat().one_hot();
                hc_alt_mask.fanout(hard<2>{});
                val<1> has_hc_alt = hc_alt_cand.fold_or();
                val<1> hc_alt_pred_val = (hc_alt_mask & gpreds[w]) != hard<0>{};

                val<1> use_alt_hc = ~use_provider_hc & has_hc_alt.fo1();
                use_alt_hc.fanout(hard<2>{});
#ifdef USE_ALT_PRED
                // Final HCpred value and source
                val<1> hc_val = select(use_provider_hc, provider[w][inst],
                                select(use_alt_hc,     hc_alt_pred_val.fo1(),
                                                       select(prov_oh.fold_or(),provider[w][inst],readb[inst])));
                val<NUMG> hc_src = select(use_provider_hc, val<NUMG>{match_provider[w][inst]},
                                   select(use_alt_hc,      hc_alt_mask,
                                                            select(prov_oh.fold_or(),val<NUMG>{match_provider[w][inst]},val<NUMG>{0})));
                val<2> hc_type = select(use_provider_hc, val<2>{1},   // provider
                                 select(use_alt_hc,       val<2>{2},   // alternate
                                                           val<2>{0})); // bimodal
#else
                       // Final HCpred value and source
                val<1> hc_val = select(prov_oh.fold_or(), provider[w][inst],
                                readb[inst]);
                val<NUMG> hc_src = select(prov_oh.fold_or(), val<NUMG>{match_provider[w][inst]},
                                   val<NUMG>{0});         
                val<2> hc_type = select(prov_oh.fold_or(), val<2>{1},   // provider
                                                           val<2>{0}); // bimodal
#endif
                hc_pred_val[w][inst]    = hc_val.fo1();
                hc_pred_source[w][inst] = hc_src.fo1();

                hc_pred_type[w][inst] = hc_type;
            }
        }
        // hc_pred_val is used once in p2 computation - no fanout needed (FO=1)

        // p2 now uses HCpred: strongest non-weak entry (provider > longest non-weak alt > bimodal)
        p2 = arr<val<1>,LINEINST>{[&](u64 offset){
            arr<val<1>,NUMWAYS> final_pred = [&](u64 w){
                return val<1>{hc_pred_val[w][offset]};
            };
            return final_pred.fo1().fold_or();
        }}.concat();



        p2.fanout(hard<LINEINST>{});
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
            match_provider[w].fanout(hard<8>{});
            provider[w].fanout(hard<4>{});
            alt[w].fanout(hard<2>{});
            readctr_cnt[w].fanout(hard<5>{});
        }
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<LINEINST+NUMWAYS>{});
        gfolds.fanout(hard<2>{});
        index1.fanout(hard<LINEINST*3>{});
        readp1.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST*3>{});
        readb.fanout(hard<4>{});


        val<LOGLINEINST> last_offset = branch_offset[num_branch-1];
        last_offset.fanout(hard<4*NUMG+2>{});

        u64 update_valid = (u64(1)<<num_branch)-1;
        arr<val<LINEINST>,LINEINST> update_mask = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){ return branch_offset[i] == offset; };
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});

        arr<val<1>,LINEINST> is_branch = [&](u64 offset){
            return update_mask[offset] != hard<0>{};
        };
        is_branch.fanout(hard<20>{});

        val<LINEINST> branch_mask = is_branch.concat();

        // ==================== P1 Disagree Check (before extra_cycle) ====================
        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask;
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
            return match_offset.fo1().fold_or();
        };
        actualdirs.fanout(hard<7>{});

        // Aggregate match_provider across branches (like tage.hpp actual_match1)
        arr<val<NUMG>,NUMWAYS> primary_mask = [&](u64 w){
            arr<val<NUMG>,LINEINST> m = [&](u64 offset){
                return val<NUMG>{is_branch[offset].replicate(hard<NUMG>{}).concat()} & 
                val<NUMG>{match_provider[w][offset]};
            };
            return m.fo1().fold_or();
        };
        primary_mask.fanout(hard<8>{});

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
            arr<val<1>,NUMWAYS> no_provider_arr = [&](u64 w){
                return match_provider[w][offset] == hard<0>{};
            };
            val<1> no_hit = no_provider_arr.fo1().fold_and();
            // Condition 2: HCpred fell back to bimodal (hc_pred_type == 0)
            arr<val<1>,NUMWAYS> hc_is_bim_arr = [&](u64 w){
                return val<2>{hc_pred_type[w][offset]} == hard<0>{};
            };
            val<1> hc_is_bim = hc_is_bim_arr.fo1().fold_or();

            // Provider weak & wrong check (for condition 2)
            arr<val<1>,NUMWAYS> prov_weak_wrong_arr = [&](u64 w){
                arr<val<1>,NUMG> prov_bits = primary_mask[w].make_array(val<1>{});
                arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                    return select(prov_bits[j], readctr_cnt[w][j], val<CTRBIT-1>{0});
                };
                val<CTRBIT-1> prov_ctr = prov_ctrs.fold_or();
                val<1> prov_pred = provider[w][offset];
                val<1> prov_weak = is_weak(prov_pred, prov_ctr);
                val<1> prov_wrong = prov_pred != actual_dir;
                val<1> has_prov = primary_mask[w] != hard<0>{};
                return has_prov & prov_weak & prov_wrong;
            };
            val<1> prov_weak_wrong = prov_weak_wrong_arr.fo1().fold_or();

            val<1> cond_extra_bim = hc_is_bim & prov_weak_wrong;

            // Combined update condition
            return branch_valid & (no_hit | cond_extra_bim);
        };

        bim_update_arr.fanout(hard<3>{});
        val<1> some_bim_update = bim_update_arr.fold_or();

        val<1> some_p1_update = (disagree_mask != hard<0>{});
        val<1> extra_cycle = some_tage_update | mispredict | some_bim_update | some_p1_update;
        extra_cycle.fanout(hard<7>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        // Count extra cycle types
        perf_extra_cycle_tage_update += static_cast<u64>(some_tage_update);
        perf_extra_cycle_mispredict += static_cast<u64>(mispredict);
        perf_extra_cycle_bim_update += static_cast<u64>(some_bim_update);
        perf_extra_cycle_p1_update += static_cast<u64>(some_p1_update);
        perf_extra_cycle_total += static_cast<u64>(extra_cycle);
#endif
#endif

        need_extra_cycle(extra_cycle);

        // ==================== Bimodal Update ====================
        // Bimodal update conditions (HCpred-based):
        // 1. No TAGE hit at all (no provider) → always update bimodal
        // 2. HCpred fell back to bimodal (provider weak & wrong) → update bimodal
        for (u64 offset = 0; offset < LINEINST; offset++) {
            val<1> bim_pred = readb[offset];
            val<1> curr_hyst = readb_low[offset];
            val<1> actual_dir = actualdirs[offset];
            val<1> dir_match = bim_pred==actual_dir;
            val<1> pred_need_update = curr_hyst == val<1>{0} & bim_update_arr[offset] & ~dir_match;
            execute_if(bim_update_arr[offset], [&](){      
                bim_low[offset].write(bindex, dir_match, extra_cycle);
            });
            execute_if(pred_need_update, [&](){
                bim_hi[offset].write(bindex, actual_dir, extra_cycle);
            });
        }

        // ==================== TAGE Allocation ====================
        // Compute allocation mask (only on misprediction)
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
        mispmask.fanout(hard<2>{});

        arr<val<NUMG>,NUMWAYS> last_match = [&](u64 w){
            arr<val<NUMG>,LINEINST> last = [&](u64 offset){
                return select(offset==last_offset,match_provider[w][offset],val<NUMG>{0});
            };
            return last.fold_or();
        };
        
        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match.append(1).concat().one_hot()-1);
        postmask.fanout(hard<2>{});

        // Check if provider was not used but correct (to skip allocation)
        arr<val<1>,NUMWAYS> skip_alloc_arr = [&](u64 w){
            arr<val<NUMG>,LINEINST> prov_for_last = [&](u64 offset){
                return select(branch_offset[num_branch-1] == offset, val<NUMG>{match_provider[w][offset]}, val<NUMG>{0});
            };
            val<NUMG> last_prov_mask = prov_for_last.fo1().fold_or();
            val<1> has = last_prov_mask != hard<0>{};

            arr<val<1>,LINEINST> used_for_last = [&](u64 offset){
                return (branch_offset[num_branch-1] == offset) & use_provider[w][offset];
            };
            val<1> used = used_for_last.fo1().fold_or();

            arr<val<1>,LINEINST> pred_for_last = [&](u64 offset){
                return (branch_offset[num_branch-1] == offset) & provider[w][offset];
            };
            val<1> pred = pred_for_last.fo1().fold_or();

            val<1> actual_last = val<1>{branch_dir[num_branch-1]};
            val<1> correct = pred == actual_last;

            return has & ~used & correct;
        };
        val<1> skip_alloc = skip_alloc_arr.fo1().fold_or();


        val<NUMG> candallocmask = postmask & notumask.fo1().fold_or() & ~(skip_alloc.replicate(hard<NUMG>{}).concat());
        // candallocmask.fanout(hard<1>{});
        val<1> alloc_fail = mispredict & (candallocmask == hard<0>{});



#ifdef RESET_UBITS
        // ==================== TAGE UBIT Reset Counter ====================
        // Increment uctr and check if reset is needed
        val<UCTRBITS> uctr_threshold = val<UCTRBITS>{1 << (UCTRBITS - 1)};
        val<1> should_reset = (uctr == uctr_threshold);

        // Update uctr: increment normally, reset to 0 when threshold reached
        uctr = select(should_reset, val<UCTRBITS>{0}, select(alloc_fail,val<UCTRBITS>{uctr + 1}, uctr));
        should_reset.fanout(hard<NUMWAYS*NUMG>{});

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_tage_skip_alloc_total += static_cast<u64>(skip_alloc);
        // Increment ubit reset counter
        perf_tage_ubit_resets_total += static_cast<u64>(should_reset);
        perf_tage_alloc_fail_total += static_cast<u64>(alloc_fail);
#endif
#endif
#endif


        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask_low = collamask.one_hot();
        collamask_low.fanout(hard<2>{});
        val<NUMG> collamask_high = (collamask ^ collamask_low).one_hot();
        val<NUMG> collamask_final = select(val<2>{static_cast<u64>(std::rand())%2}==hard<0>{}, collamask_high.fo1(), collamask_low);
        arr<val<1>,NUMG> allocate = collamask_final.reverse().make_array(val<1>{});
        // allocate.fanout(hard<1>{});
#ifdef CHEATING_MODE
        if (num_branch > 0 && static_cast<u64>(mispredict)) {
            u64 last_branch_pc = static_cast<u64>(val<64>{branch_pc[num_branch-1]});
            record_mispred_branch_pc(last_branch_pc);
            if(last_branch_pc == 15588712){
                std::cout << "act dir " << std::bitset<1>(static_cast<u64>(branch_dir[num_branch-1])) 
                          << " ,cand " << std::bitset<NUMG>(static_cast<u64>(candallocmask)) 
                          << " prov hit " << std::bitset<NUMG>(static_cast<u64>(last_match.fold_or())) 
                          << " alloc " << std::bitset<NUMG>(static_cast<u64>(collamask_final.reverse())) << std::endl;
            }
        }
#endif

        u64 way_sel = static_cast<u64>(std::rand()) % NUMWAYS;

        // ==================== TAGE Update ====================
        // CTR update, tag update, and ubit update
        // All reads before extra_cycle, all writes after extra_cycle

        //one pred max 2 way
        for (u64 w=0; w<NUMWAYS; w++) {
            // ==================== Derive HCpred-based update masks ====================
            // For each table j, determine if it should be updated based on:
            // 1. It was the HCpred source (strongest non-weak entry)
            // 2. It was the provider (longest match, always update)
            // 3. It was an alt when provider was weak & wrong

            // Aggregate HCpred source across all branches
            arr<val<NUMG>,LINEINST> hc_source_mask = [&](u64 offset){
                return select(is_branch[offset], val<NUMG>{hc_pred_source[w][offset]}, val<NUMG>{0});
            };
            val<NUMG> hc_used = hc_source_mask.fo1().fold_or();

            // Aggregate provider (always update if hit)
            arr<val<NUMG>,LINEINST> provider_mask = [&](u64 offset){
                return select(is_branch[offset], val<NUMG>{match_provider[w][offset]}, val<NUMG>{0});
            };
            val<NUMG> provider_hit = provider_mask.fo1().fold_or();

            // Compute provider weak & wrong condition for extra alt updates
            arr<val<1>,LINEINST> prov_correct_arr = [&](u64 offset){
                return is_branch[offset] & (provider[w][offset]==actualdirs[offset]);
            };
            val<1> provider_correct = prov_correct_arr.fold_or();

            // Provider weak check (reuse from existing logic)
            arr<val<1>,NUMG> primary_bits = primary_mask[w].make_array(val<1>{});
            arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                return select(primary_bits[j], readctr_cnt[w][j], val<CTRBIT-1>{0});
            };
            val<CTRBIT-1> prov_ctr = prov_ctrs.fold_or();
            arr<val<1>,LINEINST> prov_pred_arr = [&](u64 offset){
                return is_branch[offset] & provider[w][offset];
            };
            val<1> provider_pred = prov_pred_arr.fold_or();
            val<1> provider_is_weak = is_weak(provider_pred, prov_ctr);
            val<1> provider_wrong = ~provider_correct;

            // All alt hits (not just match_alt which is one-hot)
            arr<val<NUMG>,LINEINST> alt_mask = [&](u64 offset){
                
                return select(is_branch[offset], match_alt[w][offset], val<NUMG>{0});
            };
            val<NUMG> alt_hits = alt_mask.fo1().fold_or();

            // Extra alt update condition: provider weak & wrong
            val<1> base_cond = provider_is_weak & provider_wrong;
            // val<NUMG> extra_alt_update = select(extra_alt_cond, alt_hits, val<NUMG>{0});

            arr<val<1>,NUMG> hc_used_bits = hc_used.make_array(val<1>{});
            arr<val<1>,NUMG> prov_hit_bits = provider_hit.make_array(val<1>{});
            arr<val<1>,NUMG> extra_alt_bits = alt_hits.make_array(val<1>{});

            // Compute provider_correct and alt_wrong for ubit update
            arr<val<1>,LINEINST> alt_wrong_arr = [&](u64 offset){
                return is_branch[offset] & (select(alt_hits!=val<NUMG>{0},alt[w][offset]!=actualdirs[offset], readb[offset] != actualdirs[offset]));
            };
            val<1> alt_wrong = alt_wrong_arr.fold_or();

            for (u64 j=0; j<NUMG; j++) {
                // Update conditions:
                // - HCpred used this table
                // - Provider hit (always update longest match)
                // - Extra alt update (provider weak & wrong & this is an alt)
#ifdef UPDATEALTONWEAKMISP
                val<1> was_used = ((hc_used_bits[j] | extra_alt_bits[j]) & base_cond) | prov_hit_bits[j];
#else
                val<1> was_used = prov_hit_bits[j];
#endif
                val<1> do_alloc = (w == way_sel) ? allocate[j] : val<1>{0};

                // Compute actual direction for this table j
                // For allocation: use last_offset branch direction
                // For existing entry: use the branch direction from the offset that matched this table
// #ifdef HASH_TAG
                // With HASH_TAG, we need to find which offset matched this table
                arr<val<1>,LINEINST> offset_match = [&](u64 offset){
                    val<1> j_matches = match[w][offset].make_array(val<1>{})[j];
                    return is_branch[offset] & j_matches;
                };
                arr<val<1>,LINEINST> offset_dir = [&](u64 offset){
                    return select(offset_match[offset], actualdirs[offset], val<1>{0});
                };
                val<1> matched_dir = offset_dir.fo1().fold_or();
// #else
//                 // // Without HASH_TAG, extract offset from tag
//                 val<LOGLINEINST> tag_offset = readt[w][j] >> HTAGBITS;
//                 arr<val<1>,LINEINST> offset_match = [&](u64 offset){
//                     val<1> j_matches = match[w][offset].make_array(val<1>{})[j];
//                     return is_branch[offset] & j_matches;
//                 };
//                 val<1> matched_dir = (offset_match.fo1().concat() & actualdirs.concat()) != hard<0>{};
// #endif
                val<1> last_dir = branch_dir[num_branch-1];
                val<1> actual_val = select(do_alloc, last_dir, matched_dir);

                // Use stored CTR and UBIT values from predict phase
                val<CTRBIT> current_ctr = concat(readctr_pred[w][j], readctr_cnt[w][j]);
                val<UBIT> cur_u = readu[w][j];

                // CTR update: increment if correct, decrement if wrong, init on alloc
                val<CTRBIT> init_ctr = concat(actual_val,val<CTRBIT-1>{0});
                val<CTRBIT> updated_ctr = update_ctr(current_ctr, actual_val);
                val<CTRBIT> new_ctr = select(do_alloc, init_ctr, updated_ctr);

                // UBIT: alloc -> 0, provider right & alt wrong -> +1
                val<1> has_provider = primary_bits[j];
                val<1> inc_useful = has_provider & provider_correct & alt_wrong;

                val<UBIT> new_u = select(cur_u == hard<cur_u.maxval>{}, cur_u, val<UBIT>{cur_u + 1});
                val<UBIT> final_u = select(do_alloc, val<UBIT>{0}, new_u);

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
                // Count per-offset: provider used / alt used / correct, and confidence distribution at hit
                // We iterate over all LINEINST slots so we can use per-offset match/direction data.
                for (u64 offset = 0; offset < LINEINST; offset++) {
                    val<1> is_br = is_branch[offset];
                    if (!static_cast<bool>(is_br)) continue;

                    // Check if table j is the provider (longest match) for this offset
                    // val<1> is_prov_j = val<NUMG>{match_provider[w][offset]}.make_array(val<1>{})[j];
                    // Check if table j is the HC prediction source (hc_pred_source encodes the one-hot table)
                    val<1> is_hc_src_j = val<NUMG>{hc_pred_source[w][offset]}.make_array(val<1>{})[j];
                    // hc_pred_type==1: provider used as HC; hc_pred_type==2: alt used as HC
                    val<2> hc_type = val<2>{hc_pred_type[w][offset]};
                    // Table j is provider-HC: is_hc_src_j && hc_type==1
                    val<1> j_as_prov_hc = is_hc_src_j & (hc_type == val<2>{1});
                    // Table j is alt-HC: is_hc_src_j && hc_type==2
                    val<1> j_as_alt_hc  = is_hc_src_j & (hc_type == val<2>{2});

                    val<1> actual_dir = actualdirs[offset];
                    val<1> pred_dir   = readctr_pred[w][j]; // CTR MSB = prediction

                    perf_tage_provider_used[w][j]    += static_cast<u64>(j_as_prov_hc);
                    perf_tage_provider_correct[w][j] += static_cast<u64>(j_as_prov_hc & (pred_dir == actual_dir));
                    perf_tage_alt_used[w][j]         += static_cast<u64>(j_as_alt_hc);
                    perf_tage_alt_correct[w][j]      += static_cast<u64>(j_as_alt_hc  & (pred_dir == actual_dir));

                    // Confidence distribution: record CTR value when table j hits for this offset
                    // val<1> hits_j = val<NUMG>{match_provider[w][offset]}.make_array(val<1>{})[j]
                    //               | (val<NUMG>{match[w][offset]} ^ val<NUMG>{match_provider[w][offset]}).make_array(val<1>{})[j];
                    // (match[w][offset] covers all hits; match_provider is just the longest-match one)
                    // Simpler: any bit j set in match[w][offset] means table j hit for this offset
                    val<1> any_hit_j = val<NUMG>{match[w][offset]}.make_array(val<1>{})[j];
                    if (static_cast<bool>(any_hit_j & is_br)) {
                        // Full CTR = concat(pred, cnt) — value in [0, 2^CTRBIT)
                        u64 ctr_val = static_cast<u64>(concat(readctr_pred[w][j], readctr_cnt[w][j]));
                        // For CTRBIT=3: values 0-7. Map to 4 confidence levels (0=low, 3=high)
                        // Confidence = ctr_val >> (CTRBIT - 2)  maps [0,1]->0, [2,3]->1, [4,5]->2, [6,7]->3
                        u64 conf_idx = ctr_val >> (CTRBIT - 2);
                        if (conf_idx < 4) perf_tage_conf[w][j][conf_idx]++;
                    }
                }

                // Count CTR and UBIT updates (aggregate, not per-offset)
                val<1> ctr_update_cond = was_used | do_alloc;
                val<1> ubit_update_cond = do_alloc | inc_useful;
                perf_tage_ctr_updates[w][j]   += static_cast<u64>(ctr_update_cond);
                perf_tage_useful_updates[w][j] += static_cast<u64>(ubit_update_cond);
                perf_useful_alloc += static_cast<u64>(do_alloc);   
                perf_useful_inc  += static_cast<u64>(inc_useful);
                perf_useful_reset+= static_cast<u64>(should_reset);
#endif
#endif

                // Merge writes: CTR write when used or allocated
                execute_if(was_used | do_alloc, [&](){
                    gctr[w][j].write(gindex[j], new_ctr, extra_cycle);
                });

                // TAG write only on allocation
                execute_if(do_alloc, [&](){
#ifdef HASH_TAG
                    val<TAGW> new_tag = htag[j] ^ val<HTAGBITS>{last_offset};
#else
                    val<TAGW> new_tag = concat(val<LOGLINEINST>{last_offset},htag[j]);
#endif
                    gtag[w][j].write(gindex[j], new_tag,extra_cycle);
                });

                // UBIT write when allocated or incremented
#ifdef RESET_UBITS
                // When reset is triggered, write 0 to all u-bits; otherwise write normal value
                val<UBIT> ubit_value = select(should_reset, val<UBIT>{0}, final_u);
                // execute_if(should_reset,[&](){ubit[w][j].reset();});
                execute_if(do_alloc | inc_useful | should_reset, [&](){
                    ubit[w][j].write(gindex[j], ubit_value, extra_cycle);
                });
#else
                execute_if(do_alloc | inc_useful, [&](){
                    ubit[w][j].write(gindex[j], final_u, extra_cycle);
                });
#endif
            }
        }
#ifdef CHEATING_MODE
        assert(allocate.concat().ones()<=1);
#endif
        //TODO:fix it
#ifdef USE_ALT_PRED
        // USE_ALT_ON_NA update - aggregate conditions using arrays

        // Helper lambda to compute common USE_ALT_PRED conditions
        auto compute_use_alt_conditions = [&](u64 w){
            struct Conditions {
                val<1> has_provider;
                val<1> provider_is_weak;
                val<1> provider_pred;
                val<1> alt_pred;
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

            arr<val<1>,LINEINST> alt_pred_arr = [&](u64 offset){
                return is_branch[offset] & alt[w][offset];
            };
            val<1> alt_pred = alt_pred_arr.fold_or();

            // arr<val<1>,LINEINST> base_arr = [&](u64 offset){
            //     return is_branch[offset] & readb[offset];
            // };
            // val<1> base_pred = base_arr.fo1().fold_or();

            // arr<val<NUMG>,LINEINST> alt_mask = [&](u64 offset){
            //     return select(is_branch[offset], val<NUMG>{match_alt[w][offset]}, val<NUMG>{0});
            // };
            // val<1> has_alt = alt_mask.fo1().fold_or() != hard<0>{};
            // val<1> alt_or_base = select(has_alt, alt_pred, base_pred);
            val<1> actual_dir = actualdirs.concat() != hard<0>{};

            val<1> has_provider = primary_mask[w] != hard<0>{};
            val<1> provider_is_weak = is_weak(provider_pred, prov_ctr);
            // removed unused alt_or_base_match

            return Conditions{has_provider, provider_is_weak, provider_pred, alt_pred, actual_dir};
        };

        arr<val<1>,NUMWAYS> inc_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & (cond.alt_pred != cond.provider_pred) & (cond.alt_pred == cond.actual_dir);
        };

        arr<val<1>,NUMWAYS> dec_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & (cond.alt_pred != cond.provider_pred) & (cond.alt_pred != cond.actual_dir);
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
                val<1> alt_hit = (val<NUMG>{match_alt[w][offset]} != hard<0>{});

                val<1> is_prov = (hc_t == val<2>{1}) & prov_hit;
                val<1> is_alt = (hc_t == val<2>{2}) & alt_hit;

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
        dump_mispred_branch_db_csv();
#ifdef PERF_COUNTERS
        print_perf_counters();
#endif
#endif
    }
};