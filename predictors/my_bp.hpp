// this is a basic TAGE, not necessarily well optimized

#include <cmath>
#include <iostream>
#include <iomanip>
#define USE_ALT_PRED
#define RESET_UBITS

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include "tutorial/energy_monitor.hpp"
using namespace hcm;
#ifdef DEBUG_ENERGY
    struct energy_monitor monitor;
#endif
// #define PERF_COUNTERS
template<u64 LOGLB=6, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGW=11, u64 GHIST=100, u64 LOGP1=14, u64 GHIST1=6
, u64 NUMBANKS=4,u64 NUMWAYS=1,u64 CTRBIT=3,u64 UBIT=2>
struct my_bp : predictor {
    static_assert(LOGLB>2);
    static_assert(NUMG>0);
    static constexpr u64 MINHIST = 2;
    static constexpr u64 USE_ALT_PRED_BITS = 4;
    static constexpr u64 UCTRBITS = 8;
    static constexpr u64 PATHBITS = 6;

    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGB > LOGLINEINST);
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 bindex_bits = LOGB-LOGLINEINST;

    static_assert(TAGW > LOGLINEINST);
    static constexpr u64 HTAGBITS = TAGW;

    static constexpr u64 NUMGSETS = (1<<LOGB)/NUMBANKS/NUMWAYS;
    static constexpr u64 NUMBSETS = (1<<(bindex_bits))/NUMBANKS;
    static constexpr u64 NUMP1SETS = (1<<(index1_bits))/NUMBANKS;
    static constexpr u64 GINDEXBITS = static_cast<u64>(std::log2(static_cast<double>(NUMGSETS)));
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
    arr<reg<HTAGBITS>,NUMG> readt[NUMWAYS];
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

#ifdef USE_ALT_PRED
    reg<USE_ALT_PRED_BITS> use_alt_on_na;
    arr<reg<1>,LINEINST> newly_alloc;
    arr<reg<1>,LINEINST> use_provider[NUMWAYS];
#endif

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
    u64 perf_tage_correct[NUMWAYS][NUMG] = {};
    u64 perf_tage_alloc[NUMWAYS][NUMG] = {};

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

        // TAGE statistics per way per table
        for (u64 w=0; w<NUMWAYS; w++) {
            std::cerr << "\n┌─ TAGE Way " << w << " ────────────────────────────────────────────────────┐\n";
            std::cerr << "│ Table │ HistLen │ Predictions │ Correct │ Allocations │ Accuracy │\n";
            std::cerr << "├───────┼─────────┼─────────────┼─────────┼─────────────┼──────────┤\n";

            for (u64 j=0; j<NUMG; j++) {
                std::cerr << "│ " << std::setw(5) << std::left << j << " │ ";
                std::cerr << std::setw(7) << std::right << gfolds.HLEN[j] << " │ ";
                std::cerr << std::setw(11) << std::right << perf_tage_predictions[w][j] << " │ ";
                std::cerr << std::setw(7) << std::right << perf_tage_correct[w][j] << " │ ";
                std::cerr << std::setw(11) << std::right << perf_tage_alloc[w][j] << " │ ";

                if (perf_tage_predictions[w][j] > 0) {
                    double accuracy = (100.0 * perf_tage_correct[w][j]) / perf_tage_predictions[w][j];
                    std::cerr << std::fixed << std::setprecision(1) << std::setw(7) << std::right << accuracy << "% │\n";
                } else {
                    std::cerr << "    N/A │\n";
                }
            }
            std::cerr << "└───────┴─────────┴─────────────┴─────────┴─────────────┴──────────┘\n";
        }

        std::cerr << "\n";
    }
#endif

    u64 num_branch = 0;
    u64 block_size = 0;

    arr<reg<LOGLINEINST>,LINEINST> branch_offset;
    arr<reg<1>,LINEINST> branch_dir;
    reg<LINEINST> inst_oh;

    rwram<1,NUMP1SETS,NUMBANKS> table1_pred[LINEINST] {"P1 pred"};

    //205ps
    rwram<TAGW,NUMGSETS,NUMBANKS> gtag[NUMWAYS][NUMG] {"tags"};
    //183ps
    rwram<CTRBIT,NUMGSETS,NUMBANKS> gctr[NUMWAYS][NUMG] {"ctr"};
    //
    rwram<UBIT,NUMGSETS,NUMBANKS> ubit[NUMWAYS][NUMG] {"uctr"};

    //116ps
    rwram<1,NUMBSETS,NUMBANKS> bim_hi[LINEINST] {"bpred"};
    rwram<1,NUMBSETS,NUMBANKS> bim_low[LINEINST] {"bhyst"};

    zone UPDATE_ONLY;
    rwram<1,NUMP1SETS,NUMBANKS> table1_hyst[LINEINST] {"P1 hyst"};


    my_bp()
    {
        constexpr u64 lineinst = 1 << LOGLINEINST;
        constexpr u64 total_bits = (1ULL << index1_bits) * lineinst
            + (1ULL << LOGG) * (TAGW + CTRBIT + UBIT) * NUMG * NUMWAYS
            + (1ULL << bindex_bits) * 2 * lineinst;
        std::cerr << "TAGE history lengths: ";
        for (u64 i=0; i<NUMG; i++) std::cerr << gfolds.HLEN[i] << " ";
        std::cerr << std::endl;
        std::cerr << "Total storage: " << total_bits << " bits (" << (total_bits / 8192.0) << " KB)" << std::endl;
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
        val<std::max(index1_bits,GHIST1)> lineaddr = inst_pc >> LOGLB;
        lineaddr.fanout(hard<2>{});
        if constexpr (GHIST1 <= index1_bits) {
            index1 = lineaddr ^ (val<index1_bits>{path_history}<<(index1_bits-GHIST1));
        } else {
            index1 = path_history.make_array(val<index1_bits>{}).append(lineaddr).fold_xor();
        }
        index1.fanout(hard<LINEINST>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            readp1[offset] = table1_pred[offset].read(index1);
        }
        readp1.fanout(hard<2>{});
        p1 = readp1.concat();
        p1.fanout(hard<LINEINST>{});
        return (inst_oh & p1) != hard<0>{};
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        return ((inst_oh<<block_size) & p1) != hard<0>{};
    };

    val<1> is_weak(val<1> pred, val<CTRBIT-1> ctr){
        return select(pred, ctr!=ctr.maxval, ctr!=ctr.minval);
    }

    val<1> predict2(val<64> inst_pc)
    {
        val<HTAGBITS> raw_tag = inst_pc >> (2 + GINDEXBITS);
        val<GINDEXBITS> raw_gindex = inst_pc >> 2;
        val<bindex_bits> bindex_val = inst_pc >> LOGLB;

        raw_tag.fanout(hard<1+NUMG*2>{});
        raw_gindex.fanout(hard<NUMG>{});
        gfolds.fanout(hard<2>{});
        bindex_val.fanout(hard<LINEINST>{});

        for (u64 i=0; i<NUMG; i++) {
            gindex[i] = raw_gindex ^ gfolds.template get<0>(i);
        }
        gindex.fanout(hard<3>{});

        for (u64 i=0; i<NUMG; i++) {
            htag[i] = raw_tag ^ gfolds.template get<1>(i);
        }
        htag.fanout(hard<2>{});

        for (u64 offset=0; offset<LINEINST; offset++) {
            readb[offset] = bim_hi[offset].read(bindex_val);
            readb_low[offset] = bim_low[offset].read(bindex_val);
        }
        readb.fanout(hard<2>{});
        readb_low.fanout(hard<2>{});

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
            readt[w].fanout(hard<LINEINST+1>{});
            readu[w].fanout(hard<NUMG>{});
        }

        //20ps
        arr<val<NUMG>,NUMWAYS> umask = [&](u64 w){
            arr<val<1>,NUMG> u = [&](u64 j){ return readu[w][j]!=hard<0>{}; };
            return u.concat();
        };

        //TODO: fix this
        for(u64 way=0;way<NUMWAYS;way++){
            notumask[way] = ~umask[way];
        }
        // notumask.fanout(hard<2>{});

        //12ps 2000fJ
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                arr<val<1>,NUMG> wayhit = [&](u64 j){
                    return readt[w][j] == (htag[j] ^ val<HTAGBITS>{offset});
                };
                match[w][offset] = wayhit.concat();
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) match[w].fanout(hard<2>{});

        //20ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                match_provider[w][offset] = match[w][offset].one_hot();
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) match_provider[w].fanout(hard<3>{});

        arr<val<NUMG+1>,NUMWAYS> gpreds = [&](u64 w){ return readctr_pred[w].concat(); };
        // arr<val<NUMG+1>,LINEINST> preds = [&](u64 offset){return concat(readb[offset],gpreds);};
        gpreds.fanout(hard<2>{});
        //30ps
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 offset=0; offset<LINEINST; offset++) {
                provider[w][offset] = (match_provider[w][offset] & gpreds[w]) != hard<0>{};
            }
        }
        for (u64 w=0; w<NUMWAYS; w++) provider[w].fanout(hard<2>{});

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

#ifdef USE_ALT_PRED
        use_alt_on_na.fanout(hard<2>{});
        arr<val<1>,USE_ALT_PRED_BITS> use_alt_on_na_array = use_alt_on_na.make_array(val<1>{});
        val<1> use_alt_on_na_pos = use_alt_on_na_array[USE_ALT_PRED_BITS-1];
        for (u64 w=0; w<NUMWAYS; w++) {
            for (u64 inst=0; inst<LINEINST; inst++) {
                arr<val<1>,NUMG> prov_oh = val<NUMG>{match_provider[w][inst]}.make_array(val<1>{});
                arr<val<1>,NUMG> alt_oh = val<NUMG>{match_alt[w][inst]}.make_array(val<1>{});
                arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                    return select(prov_oh[j],readctr_cnt[w][j],val<CTRBIT-1>{0});
                };
                val<CTRBIT-1> prov_ctr  = prov_ctrs.fold_or();
                // provider & alt 
                use_provider[w][inst] = prov_oh.fold_or() & 
                select(alt_oh.fold_or(),(~is_weak(provider[w][inst], prov_ctr) | ~use_alt_on_na_pos),val<1>{1});
            }
        }
        // use_provider.fanout(hard<2>{});

        p2 = arr<val<1>,LINEINST>{[&](u64 offset){
            arr<val<1>,NUMWAYS> final_pred = [&](u64 w){
                val<1> has_alt = match_alt[w][offset] != hard<0>{};
                return select(use_provider[w][offset], provider[w][offset],select(has_alt,alt[w][offset],readb[offset]));
            };
            return final_pred.fo1().fold_or();
        }}.concat();

#else
        p2 = arr<val<1>,LINEINST>{[&](u64 offset){
            arr<val<1>,NUMWAYS> final_pred = [&](u64 w){ return provider[w][offset]; };
            return select(final_pred.fo1().fold_or(), val<1>{1}, readb[offset]);
        }}.concat();
#endif

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
        val<LOGLINEINST> offset = branch_pc.fo1() >> 2;
        branch_offset[num_branch] = offset;
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

        mispredict.fanout(hard<3>{});
        val<1> correct_pred = ~mispredict;
        // correct_pred.fanout(hard<2>{});
        gindex.fanout(hard<3>{});
        // htag.fanout(hard<1>{});
        for (u64 w=0; w<NUMWAYS; w++) {
            match_provider[w].fanout(hard<5>{});
            match_alt[w].fanout(hard<3>{});
            provider[w].fanout(hard<3>{});
            alt[w].fanout(hard<2>{});
            use_provider[w].fanout(hard<5>{});
            readctr_cnt[w].fanout(hard<3>{});
            // notumask[w].fanout(hard<1>{});
        }
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<LINEINST+NUMWAYS>{});
        gfolds.fanout(hard<2>{});
        index1.fanout(hard<LINEINST*3>{});
        // p1.fanout(hard<2>{});
        p2.fanout(hard<2*LINEINST>{});
        readp1.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST*3>{});
        readb.fanout(hard<2>{});
        readb_low.fanout(hard<2>{});


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
        is_branch.fanout(hard<6>{});

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
        actualdirs.fanout(hard<6>{});

        // Aggregate match_provider across branches (like tage.hpp actual_match1)
        arr<val<NUMG>,NUMWAYS> primary_mask = [&](u64 w){
            arr<val<NUMG>,LINEINST> m = [&](u64 offset){
                return val<NUMG>{is_branch[offset].replicate(hard<NUMG>{}).concat()} & 
                val<NUMG>{match_provider[w][offset]};
            };
            return m.fo1().fold_or();
        };
        primary_mask.fanout(hard<5>{});

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

        // Compute actual branch directions for update purposes
        // arr<val<1>,LINEINST> branch_taken = [&](u64 offset) {
        //     // Extract the actual direction from branch_dir array
        //     return val<1>{branch_dir[offset]};
        // };

        // Calculate bimodal predictions
        val<LINEINST> bim_predictions = readb.concat();

        // Determine which updates are needed for extra cycle
        val<1> some_tage_update = primary_mask.fold_or() != hard<0>{};
        val<1> some_bim_update = ((bim_predictions ^ actualdirs.concat()) & branch_mask) != hard<0>{};
        val<1> some_p1_update = (disagree_mask != hard<0>{});
        val<1> extra_cycle = some_tage_update | mispredict | some_bim_update | some_p1_update;
        extra_cycle.fanout(hard<7>{});
        need_extra_cycle(extra_cycle);

        // ==================== Bimodal Update ====================
        // Bimodal update conditions:
        // 1. pred wrong
        // 2. pred right and weak (hysteresis == 0)

        // Perform updates for each offset
        for (u64 offset = 0; offset < LINEINST; offset++) {
            val<1> branch_valid = is_branch[offset];
            val<1> bim_pred = readb[offset];
            val<1> curr_hyst = readb_low[offset];  // Use stored hysteresis
            val<1> actual_dir = actualdirs[offset];
            val<1> pred_correct = (bim_pred == actual_dir);

            // Compute conditions
            val<1> pred_weak = (curr_hyst == hard<0>{});
            val<1> should_update = branch_valid & (~pred_correct | (pred_correct & pred_weak));

            // Write (after extra_cycle), merged into single execute_if
            execute_if(should_update, [&](){
                val<1> new_pred = update_ctr(bim_pred, actual_dir);
                val<1> new_hyst = update_ctr(curr_hyst, pred_correct);

                bim_hi[offset].write(bindex, new_pred, extra_cycle);
                bim_low[offset].write(bindex, new_hyst, extra_cycle);
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
        

        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match.fo1().fold_or()-1);
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

        val<NUMG> candallocmask = postmask & notumask.fold_or() & ~(skip_alloc.replicate(hard<NUMG>{}).concat());
        // candallocmask.fanout(hard<1>{});

        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask_low = collamask.one_hot();
        collamask_low.fanout(hard<2>{});
        val<NUMG> collamask_high = (collamask ^ collamask_low).one_hot();
        val<NUMG> collamask_final = select(val<2>{static_cast<u64>(std::rand())%2}==hard<0>{}, collamask_high.fo1(), collamask_low);
        arr<val<1>,NUMG> allocate = collamask_final.fo1().reverse().make_array(val<1>{});
        // allocate.fanout(hard<1>{});


        u64 way_sel = static_cast<u64>(std::rand()) % NUMWAYS;

        // ==================== TAGE Update ====================
        // CTR update, tag update, and ubit update
        // All reads before extra_cycle, all writes after extra_cycle

        //one pred max 2 way
        for (u64 w=0; w<NUMWAYS; w++) {
            // Compute provider_correct and alt_wrong for this way
            arr<val<1>,LINEINST> prov_correct_arr = [&](u64 offset){
                return is_branch[offset] & (provider[w][offset]==actualdirs[offset]);
            };

            arr<val<1>,LINEINST> alt_wrong_arr = [&](u64 offset){
                return is_branch[offset] & (alt[w][offset]!=actualdirs[offset]);
            };

            // val<1> actual = actualdirs != hard<0>{};
            val<1> provider_correct = prov_correct_arr.fold_or();
            val<1> alt_wrong = alt_wrong_arr.fold_or();

            // Compute masks for CTR update
            arr<val<NUMG>,LINEINST> provider_used_mask = [&](u64 offset){
                return select(is_branch[offset] & use_provider[w][offset], val<NUMG>{match_provider[w][offset]}, val<NUMG>{0});
            };
            arr<val<NUMG>,LINEINST> alt_used_mask = [&](u64 offset){
                return select(is_branch[offset] & ~use_provider[w][offset], val<NUMG>{match_alt[w][offset]}, val<NUMG>{0});
            };

            //at least one pred ,or bim pred fix it
            arr<val<NUMG>,LINEINST> dir_mask = [&](u64 offset){
                val<NUMG> dir_rep = val<1>{actualdirs[offset]}.replicate(hard<NUMG>{}).concat();
                val<NUMG> prov_dir = select(is_branch[offset] & use_provider[w][offset], val<NUMG>{match_provider[w][offset]} & dir_rep, val<NUMG>{0});
                val<NUMG> alt_dir = select(is_branch[offset] & ~use_provider[w][offset], val<NUMG>{match_alt[w][offset]} & dir_rep, val<NUMG>{0});
                return prov_dir | alt_dir ;
            };


            val<NUMG> provider_used = provider_used_mask.fo1().fold_or();
            val<NUMG> alt_used = alt_used_mask.fo1().fold_or();
            val<NUMG> actual_dir = dir_mask.fo1().fold_or();

            arr<val<1>,NUMG> prov_used_bits = provider_used.make_array(val<1>{});
            arr<val<1>,NUMG> alt_used_bits = alt_used.make_array(val<1>{});
            arr<val<1>,NUMG> actual_bits = actual_dir.make_array(val<1>{});

            for (u64 j=0; j<NUMG; j++) {
                val<1> was_used = prov_used_bits[j] | alt_used_bits[j];
                val<1> actual_val = actual_bits[j];
                val<1> do_alloc = (w == way_sel) ? allocate[j] : val<1>{0};

                // Use stored CTR and UBIT values from predict phase
                val<CTRBIT> current_ctr = concat(readctr_pred[w][j], readctr_cnt[w][j]);
                val<UBIT> cur_u = readu[w][j];

                // CTR update: increment if correct, decrement if wrong, init on alloc
                val<CTRBIT> init_ctr = val<CTRBIT>{1 << (CTRBIT - 1)};
                val<CTRBIT> updated_ctr = update_ctr(current_ctr, actual_val);
                val<CTRBIT> new_ctr = select(do_alloc, init_ctr, updated_ctr);

                // UBIT: alloc -> 0, provider right & alt wrong -> +1
                arr<val<1>,NUMG> primary = primary_mask[w].make_array(val<1>{});
                val<1> has_provider = primary[j];
                val<1> inc_useful = has_provider & provider_correct & alt_wrong;

                val<UBIT> new_u = select(cur_u == hard<cur_u.maxval>{}, cur_u, val<UBIT>{cur_u + 1});
                val<UBIT> final_u = select(do_alloc, val<UBIT>{0}, new_u);

                // Merge writes: CTR write when used or allocated
                execute_if(was_used | do_alloc, [&](){
                    gctr[w][j].write(gindex[j], new_ctr, extra_cycle);
                });

                // TAG write only on allocation
                execute_if(do_alloc, [&](){
                    val<TAGW> new_tag = htag[j] ^ val<HTAGBITS>{last_offset};
                    gtag[w][j].write(gindex[j], new_tag, extra_cycle);
                });

                // UBIT write when allocated or incremented
                execute_if(do_alloc | inc_useful, [&](){
                    ubit[w][j].write(gindex[j], final_u, extra_cycle);
                });
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
                val<1> alt_or_base_match;
            };

            arr<val<1>,NUMG> primary_bits = primary_mask[w].make_array(val<1>{});
            arr<val<CTRBIT-1>,NUMG> prov_ctrs = [&](u64 j){
                return select(primary_bits[j], readctr_cnt[w][j], val<CTRBIT-1>{0});
            };
            val<CTRBIT-1> prov_ctr = prov_ctrs.fo1().fold_or();

            arr<val<1>,LINEINST> prov_pred = [&](u64 offset){
                return is_branch[offset] & provider[w][offset];
            };
            val<1> provider_pred = prov_pred.fo1().fold_or();

            arr<val<1>,LINEINST> alt_pred_arr = [&](u64 offset){
                return is_branch[offset] & alt[w][offset];
            };
            val<1> alt_pred = alt_pred_arr.fo1().fold_or();

            arr<val<1>,LINEINST> base_arr = [&](u64 offset){
                return is_branch[offset] & readb[offset];
            };
            val<1> base_pred = base_arr.fo1().fold_or();

            arr<val<NUMG>,LINEINST> alt_mask = [&](u64 offset){
                return select(is_branch[offset], val<NUMG>{match_alt[w][offset]}, val<NUMG>{0});
            };
            val<1> has_alt = alt_mask.fo1().fold_or() != hard<0>{};
            val<1> alt_or_base = select(has_alt, alt_pred, base_pred);
            val<1> actual = actualdirs.concat() != hard<0>{};

            val<1> has_provider = primary_mask[w] != hard<0>{};
            val<1> provider_is_weak = is_weak(provider_pred, prov_ctr);
            val<1> alt_or_base_match = alt_or_base == actual;

            return Conditions{has_provider, provider_is_weak, alt_or_base_match};
        };

        arr<val<1>,NUMWAYS> inc_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & cond.alt_or_base_match;
        };

        arr<val<1>,NUMWAYS> dec_use_alt_arr = [&](u64 w){
            auto cond = compute_use_alt_conditions(w);
            return cond.has_provider & cond.provider_is_weak & ~cond.alt_or_base_match;
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
        // arr<val<1>,LINEINST> p2_split = p2.make_array(val<1>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(p1_weak[offset].fo1(), [&](){
                table1_pred[offset].write(index1, actualdirs[offset],extra_cycle);
            });
        }

        // Update P1 hysteresis for all branches
        for (u64 offset=0; offset<LINEINST; offset++) {
            val<1> pred_correct = readp1[offset] == actualdirs[offset];
            execute_if(is_branch[offset], [&](){
                table1_hyst[offset].write(index1, pred_correct,extra_cycle);
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
                // Count predictions when this table was provider
                arr<val<1>,NUMG> prov_bits = primary_mask[w].make_array(val<1>{});
                val<1> was_provider = prov_bits[j];
                perf_tage_predictions[w][j] += static_cast<u64>(was_provider);

                // Check if prediction was correct
                val<1> pred_val = readctr_pred[w][j];
                arr<val<1>,LINEINST> correct_arr = [&](u64 offset){
                    arr<val<1>,NUMG> match_bits = match_provider[w][offset].make_array(val<1>{});
                    return match_bits[j] & (pred_val == actualdirs[offset]);
                };
                val<1> was_correct = correct_arr.fo1().fold_or();
                perf_tage_correct[w][j] += static_cast<u64>(was_provider & was_correct);

                // Count allocations
                val<1> do_alloc = (w == way_sel) ? allocate[j] : val<1>{0};
                perf_tage_alloc[w][j] += static_cast<u64>(do_alloc);
            }
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