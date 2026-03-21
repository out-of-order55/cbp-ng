#define USE_META
#define RESET_UBITS

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include <iomanip>

using namespace hcm;

template<u64 LOGLB=6, u64 LOGN=3, u64 NUMG=8, u64 LOGG=11, u64 LOGB=12, u64 TAGHBITS=8, u64 GHIST=100, u64 LOGP1=14, u64 GHIST1=6>
struct tage64b8 : predictor {
    static_assert(LOGLB == 6, "tage64b8 is fixed to 64B blocks");
    static_assert(LOGN == 3, "tage64b8 is fixed to 8 branches per block");
    static_assert(NUMG > 0);
    static_assert(LOGP1 > LOGN);

    static constexpr u64 N = u64(1) << LOGN;
    static constexpr u64 LOGLINEINST = LOGLB - 2;
    static constexpr u64 LINEINST = u64(1) << LOGLINEINST;
    static constexpr u64 MINHIST = 2;
    static constexpr u64 PATHBITS = 6;
    static constexpr u64 TAGW = TAGHBITS + LOGN; // bank_id + hash_tag
    static constexpr u64 LOGNUMG = (NUMG <= 1) ? 1 : std::bit_width(NUMG - 1);
    static constexpr u64 index1_bits = LOGP1 - LOGN;
    static constexpr u64 bindex_bits = LOGB - LOGN;
    static constexpr u64 METABITS = 4;
    static constexpr u64 UCTRBITS = 8;
#ifdef USE_META
    static constexpr u64 METAPIPE = 2;
#endif

    geometric_folds<NUMG,MINHIST,GHIST,LOGG,TAGHBITS> gfolds;
    reg<1> true_block = 1;

    // P1 state (8-lane package)
    reg<GHIST1> global_history1;
    reg<index1_bits> index1;
    reg<N> X; // gshareN-style bank permutation base (one-hot)
    arr<reg<LOGN>,N> bank_id_rank; // rank -> physical bank id
    arr<reg<1>,N> readp1; // unordered (by bank id)
    arr<reg<1>,N> p1_rank;

    // block control
    reg<LINEINST> block_entry;
    reg<N> rank;
    u64 num_branch = 0;
    u64 block_size = 0;
    arr<reg<1>,N> branch_dir;

    // P2 index/tag/cache
    reg<bindex_bits> bindex;
    arr<reg<LOGG>,NUMG> gindex;
    arr<reg<TAGHBITS>,NUMG> htag;

    arr<reg<TAGW>,NUMG> readt;
    arr<reg<1>,NUMG> readc;
    arr<reg<2>,NUMG> readh;
    arr<reg<1>,NUMG> readu;
    arr<reg<1>,N> readb; // unordered (by bank id)

    arr<reg<1>,N> p2_rank;
    arr<reg<1>,N> pred1_rank;
    arr<reg<1>,N> pred2_rank;
    arr<reg<1>,N> alt_valid_rank;
#ifdef USE_META
    arr<reg<1>,N> newly_alloc_rank;
    arr<reg<METABITS,i64>,METAPIPE> meta;
#endif
    arr<reg<1>,N> provider_hit_rank;
    arr<reg<LOGNUMG>,N> provider_idx_rank;
    arr<reg<1>,N> alt_provider_hit_rank;
    arr<reg<LOGNUMG>,N> alt_provider_idx_rank;
    arr<reg<NUMG>,N> provider_mask_rank;

#ifdef RESET_UBITS
    reg<UCTRBITS> uctr;
#endif

#ifdef PERF_COUNTERS
    u64 perf_predictions = 0;
    u64 perf_correct = 0;
    u64 perf_bimodal_used = 0;
    u64 perf_bimodal_correct = 0;
    u64 perf_alt_used = 0;
    u64 perf_alt_correct = 0;
    u64 perf_alt_bimodal_used = 0;
    u64 perf_alt_bimodal_correct = 0;
    u64 perf_alt_provider_used[NUMG] = {};
    u64 perf_alt_provider_correct[NUMG] = {};
    u64 perf_provider_used[NUMG] = {};
    u64 perf_provider_correct[NUMG] = {};
    u64 perf_table_reads[NUMG] = {};
    u64 perf_table_hits[NUMG] = {};
    u64 perf_table_alloc[NUMG] = {};
    u64 perf_alloc_failures = 0;
    u64 perf_alloc_fail_highest = 0;
    u64 perf_alloc_fail_noubit = 0;
    u64 perf_extra_cycle_total = 0;
    u64 perf_extra_cycle_badpred = 0;
    u64 perf_extra_cycle_mispredict = 0;
    u64 perf_extra_cycle_p1_update = 0;
    u64 perf_predict2_calls = 0;
    u64 perf_reuse_predict2_calls = 0;
    u64 perf_update_condbr_calls = 0;
    u64 perf_predict2_true_block = 0;
    u64 perf_predict2_reused_block = 0;
    u64 perf_reuse_keep_predict2 = 0;
    u64 perf_reuse_stop_predict2 = 0;
    u64 perf_reuse_keep_reuse_predict2 = 0;
    u64 perf_reuse_stop_reuse_predict2 = 0;
    u64 perf_reuse_stop_update_line_end = 0;
    u64 perf_reuse_stop_update_budget = 0;
    u64 perf_blocks_with_branch = 0;
    u64 perf_blocks_budget8 = 0;
    u64 perf_block_branch_hist[N+1] = {};
    u64 perf_meta_updates = 0;
    u64 perf_meta_good_alt = 0;
    u64 perf_meta_bad_alt = 0;
#endif

    // P1 tables: 8 predictions per block
    ram<val<1>,(1<<index1_bits)> table1_pred[N] {"P1 pred"};
    zone UPDATE_ONLY;
    ram<val<1>,(1<<index1_bits)> table1_hyst[N] {"P1 hyst"};

    // P2 tables: single-lane entries (no bank dimension)
    ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"gtag"};
    ram<val<1>,(1<<LOGG)> gpred[NUMG] {"gpred"};
    rwram<2,(1<<LOGG),4> ghyst[NUMG] {"ghyst"};
    rwram<1,(1<<LOGG),4> ubit[NUMG] {"ubit"};

    // bimodal limited to 8 lanes
    ram<val<1>,(1<<bindex_bits)> bim[N] {"bim8"};
    zone UPDATE_ONLY_BIM;
    ram<val<1>,(1<<bindex_bits)> bhyst[N] {"bhyst8"};

#ifdef PERF_COUNTERS
    void print_perf_counters()
    {
        std::cerr << "\n=== TAGE64B8 PERF ===\n";
        std::cerr << "predictions          : " << perf_predictions << "\n";
        std::cerr << "correct              : " << perf_correct << "\n";
        if (perf_predictions != 0) {
            std::cerr << "accuracy             : "
                      << std::setprecision(4)
                      << (100.0 * double(perf_correct) / double(perf_predictions)) << "%\n";
        }
        std::cerr << "bimodal used/correct : " << perf_bimodal_used << " / " << perf_bimodal_correct << "\n";
        std::cerr << "alt used/correct     : " << perf_alt_used << " / " << perf_alt_correct << "\n";
        std::cerr << "alt bim used/correct : " << perf_alt_bimodal_used << " / " << perf_alt_bimodal_correct << "\n";
        std::cerr << "alloc failures       : " << perf_alloc_failures
                  << " (highest=" << perf_alloc_fail_highest
                  << ", no_ubit=" << perf_alloc_fail_noubit << ")\n";
        std::cerr << "extra cycles total   : " << perf_extra_cycle_total
                  << " (badpred=" << perf_extra_cycle_badpred
                  << ", mispred=" << perf_extra_cycle_mispredict
                  << ", p1_update=" << perf_extra_cycle_p1_update << ")\n";
        std::cerr << "predict2/reuse2/upd  : " << perf_predict2_calls
                  << " / " << perf_reuse_predict2_calls
                  << " / " << perf_update_condbr_calls << "\n";
        std::cerr << "predict2 true/reused : " << perf_predict2_true_block
                  << " / " << perf_predict2_reused_block << "\n";
        std::cerr << "reuse keep/stop p2   : " << perf_reuse_keep_predict2
                  << " / " << perf_reuse_stop_predict2 << "\n";
        std::cerr << "reuse keep/stop rp2  : " << perf_reuse_keep_reuse_predict2
                  << " / " << perf_reuse_stop_reuse_predict2 << "\n";
        std::cerr << "reuse stop upd(le/b) : " << perf_reuse_stop_update_line_end
                  << " / " << perf_reuse_stop_update_budget << "\n";
        std::cerr << "blocks(branch/b8)    : " << perf_blocks_with_branch
                  << " / " << perf_blocks_budget8 << "\n";
        std::cerr << "meta updates(g/b)    : " << perf_meta_updates
                  << " (" << perf_meta_good_alt << "/" << perf_meta_bad_alt << ")\n";
        std::cerr << "block-branch hist    :";
        for (u64 k=1; k<=N; k++) {
            std::cerr << " " << k << ":" << perf_block_branch_hist[k];
        }
        std::cerr << "\n";
        for (u64 j=0; j<NUMG; j++) {
            std::cerr << "T" << j
                      << " used/correct=" << perf_provider_used[j] << "/" << perf_provider_correct[j]
                      << " alt_used/correct=" << perf_alt_provider_used[j] << "/" << perf_alt_provider_correct[j]
                      << " reads/hits/alloc=" << perf_table_reads[j] << "/" << perf_table_hits[j] << "/" << perf_table_alloc[j]
                      << "\n";
        }
        std::cerr << "=== END PERF ===\n";
    }

    ~tage64b8()
    {
        print_perf_counters();
    }
#endif

    val<1> line_end()
    {
        return block_entry >> (LINEINST-block_size);
    }

    val<1> last_pred()
    {
        return rank.rotate_left(num_branch);
    }

    val<1> predict1(val<64> inst_pc)
    {
        block_entry = select(true_block,
                             val<LOGLINEINST>{inst_pc>>2}.decode().concat(),
                             block_entry<<block_size);
        rank = select(true_block,
                      val<N>{1},
                      rank.rotate_left(num_branch));
        X = select(true_block,
                   val<LOGN>{inst_pc>>2}.decode().concat(),
                   X.rotate_left(num_branch));

        execute_if(true_block, [&](){
            val<std::max(index1_bits,GHIST1)> lineaddr = inst_pc >> LOGLB;
            if constexpr (GHIST1 <= index1_bits) {
                index1 = lineaddr ^ (val<index1_bits>{global_history1} << (index1_bits-GHIST1));
            } else {
                index1 = global_history1.make_array(val<index1_bits>{}).append(lineaddr).fold_xor();
            }

            for (u64 r=0; r<N; r++) {
                readp1[r] = table1_pred[r].read(index1);
            }
        });

        val<N> p1_unordered = readp1.concat();
        for (u64 r=0; r<N; r++) {
            val<N> bank_mask = X.rotate_left(r);
            p1_rank[r] = (p1_unordered & bank_mask) != hard<0>{};
            bank_id_rank[r] = arr<val<LOGN>,N>{[&](u64 b) -> val<LOGN> {
                val<1> hit_b = (bank_mask >> b) & val<1>{1};
                return select(hit_b, val<LOGN>{b}, val<LOGN>{0});
            }}.fold_or();
        }

        // keep previous block branch count for next-cycle true_block carry logic
        block_size = 1;
        return p1_rank[num_branch];
    }

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        val<1> pred = p1_rank[num_branch];
        return pred;
    }

    val<1> predict2(val<64> inst_pc)
    {
        execute_if(true_block, [&](){
            val<std::max(bindex_bits,LOGG)> lineaddr = inst_pc >> LOGLB;
            bindex = lineaddr;

            for (u64 j=0; j<NUMG; j++) {
                gindex[j] = lineaddr ^ gfolds.template get<0>(j);
                htag[j] = val<TAGHBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(j);
                readt[j] = gtag[j].read(gindex[j]);
                readc[j] = gpred[j].read(gindex[j]);
                readh[j] = ghyst[j].read(gindex[j]);
                readu[j] = ubit[j].read(gindex[j]);
            }

            for (u64 r=0; r<N; r++) {
                readb[r] = bim[r].read(bindex);
            }
        });

        val<NUMG> gpreds = readc.concat();
        val<NUMG> weakctr = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return readh[j] == hard<0>{};
        }}.concat();
        val<NUMG> coldctr = ~readu.concat() & weakctr;
#ifdef USE_META
        val<1> metasign = meta[METAPIPE-1] >= hard<0>{};
#endif

        for (u64 r=0; r<N; r++) {
            val<NUMG> hash_match = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
                return val<TAGHBITS>{readt[j]} == htag[j];
            }}.concat();
            val<NUMG> rank_match = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
                return val<LOGN>{readt[j] >> TAGHBITS} == bank_id_rank[r];
            }}.concat();
            val<NUMG> table_match = hash_match & rank_match;

            val<NUMG+1> match = concat(val<1>{1}, table_match);
            val<NUMG+1> match1 = match.one_hot();
            val<NUMG+1> match2 = (match ^ match1).one_hot();
            val<1> bim_pred_r = (readb.concat() & X.rotate_left(r)) != hard<0>{};
            val<NUMG+1> preds = concat(bim_pred_r, gpreds);

            val<NUMG> pmask = val<NUMG+1>{match1} & hard<(u64(1) << NUMG)-1>{};
            val<NUMG> amask = val<NUMG+1>{match2} & hard<(u64(1) << NUMG)-1>{};
            pred1_rank[r] = (match1 & preds) != hard<0>{};
            pred2_rank[r] = (match2 & preds) != hard<0>{};
            alt_valid_rank[r] = match2 != hard<0>{};
            provider_mask_rank[r] = pmask;
            provider_hit_rank[r] = pmask != hard<0>{};
            provider_idx_rank[r] = arr<val<LOGNUMG>,NUMG>{[&](u64 j) -> val<LOGNUMG> {
                val<1> hit_j = (pmask >> j) & val<1>{1};
                return select(hit_j, val<LOGNUMG>{j}, val<LOGNUMG>{0});
            }}.fold_or();
            alt_provider_hit_rank[r] = amask != hard<0>{};
            alt_provider_idx_rank[r] = arr<val<LOGNUMG>,NUMG>{[&](u64 j) -> val<LOGNUMG> {
                val<1> hit_j = (amask >> j) & val<1>{1};
                return select(hit_j, val<LOGNUMG>{j}, val<LOGNUMG>{0});
            }}.fold_or();
#ifdef USE_META
            newly_alloc_rank[r] = (pmask & coldctr) != hard<0>{};
            val<1> use_alt = metasign & newly_alloc_rank[r] & alt_valid_rank[r];
            p2_rank[r] = select(use_alt, pred2_rank[r], pred1_rank[r]);
#else
            p2_rank[r] = pred1_rank[r];
#endif
        }

        val<1> taken = p2_rank[num_branch];
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_predict2_calls++;
        if (static_cast<bool>(true_block)) perf_predict2_true_block++;
        else perf_predict2_reused_block++;
#endif
#endif
        val<1> keep_reuse = ~val<1>{block_entry>>(LINEINST-1)};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        if (static_cast<bool>(keep_reuse)) perf_reuse_keep_predict2++;
        else perf_reuse_stop_predict2++;
#endif
#endif
        reuse_prediction(keep_reuse);
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = p2_rank[num_branch];
        val<1> keep_reuse = ~val<1>{block_entry>>(LINEINST-1-block_size)};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_reuse_predict2_calls++;
        if (static_cast<bool>(keep_reuse)) perf_reuse_keep_reuse_predict2++;
        else perf_reuse_stop_reuse_predict2++;
#endif
#endif
        reuse_prediction(keep_reuse);
        block_size++;
        return taken;
    }

    void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken, [[maybe_unused]] val<64> next_pc)
    {
        assert(num_branch < N);
        branch_dir[num_branch] = taken;
        num_branch++;
        val<1> line_end_now = line_end();
        val<1> budget_end = last_pred();
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_update_condbr_calls++;
        if (static_cast<bool>(line_end_now)) perf_reuse_stop_update_line_end++;
        if (static_cast<bool>(budget_end)) perf_reuse_stop_update_budget++;
#endif
#endif
        reuse_prediction(~(line_end_now | budget_end));
    }

    void update_cycle(instruction_info &block_end_info)
    {
        val<1> &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc = block_end_info.next_pc;

        if (num_branch == 0) {
            val<1> line_end = block_entry >> (LINEINST - block_size);
            val<1> actual_block = ~(true_block & line_end);
            execute_if(actual_block, [&](){
                val<PATHBITS> path = next_pc >> 2;
                global_history1 = (global_history1 << 1) ^ val<GHIST1>{path};
                gfolds.update(path);
                true_block = 1;
            });
            return;
        }

        val<1> correct_pred = ~mispredict;
        u64 lr = num_branch - 1;
        val<LOGN> last_bank = bank_id_rank[lr];
        val<1> last_actual = branch_dir[lr];

        u64 valid_u64 = (u64(1) << num_branch) - 1;
        val<N> branch_mask = val<N>{valid_u64};
        arr<val<1>,N> is_branch = branch_mask.make_array(val<1>{});
        arr<val<1>,N> branch_taken = arr<val<1>,N>{[&](u64 r) -> val<1> {
            return branch_dir[r];
        }};
#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        perf_blocks_with_branch++;
        if (num_branch <= N) perf_block_branch_hist[num_branch]++;
        if (num_branch == N) perf_blocks_budget8++;
#endif
#endif

        // build bank-ordered views (bank_id -> valid/dir/pred2)
        arr<val<1>,N> bank_valid = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return arr<val<1>,N>{[&](u64 r) -> val<1> {
                return is_branch[r] & (bank_id_rank[r] == val<LOGN>{b});
            }}.concat() != hard<0>{};
        }};
        arr<val<1>,N> branch_taken_bank = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return arr<val<1>,N>{[&](u64 r) -> val<1> {
                return is_branch[r] & (bank_id_rank[r] == val<LOGN>{b}) & branch_taken[r];
            }}.concat() != hard<0>{};
        }};
        arr<val<1>,N> pred2_bank = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return arr<val<1>,N>{[&](u64 r) -> val<1> {
                return is_branch[r] & (bank_id_rank[r] == val<LOGN>{b}) & pred2_rank[r];
            }}.concat() != hard<0>{};
        }};

        // primary provider per table (tage: fold-or over actual match1)
        arr<val<1>,NUMG> primary = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return arr<val<1>,N>{[&](u64 r) -> val<1> {
                val<1> hit_j = (provider_mask_rank[r] >> j) & val<1>{1};
                return is_branch[r] & hit_j;
            }}.concat() != hard<0>{};
        }};
        val<NUMG> primary_mask = primary.concat();

        // allocation candidates (tage style), but offset field is replaced by bank_id
        arr<val<1>,NUMG> last_tagcmp = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return readt[j] == concat(last_bank, htag[j]);
        }};
        val<NUMG+1> last_match1 = last_tagcmp.fo1().append(1).concat().one_hot();
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
        val<NUMG> postmask = mispmask & val<NUMG>{last_match1 - 1};
        val<NUMG> notumask = ~readu.concat();
        val<NUMG> candallocmask = postmask & notumask;
        val<NUMG> collamask = candallocmask.reverse();
        val<NUMG> collamask1 = collamask.one_hot();
        val<NUMG> collamask2 = (collamask ^ collamask1).one_hot();
        val<NUMG> collamask12 = select(val<2>{std::rand()} == hard<0>{}, collamask2, collamask1);
        arr<val<1>,NUMG> allocate = collamask12.reverse().make_array(val<1>{});

        val<1> noalloc = candallocmask == hard<0>{};
        val<NUMG> uclearmask = postmask & noalloc.replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.make_array(val<1>{});

        // entry owner bank_id is encoded in tag upper bits
        arr<val<LOGN>,NUMG> owner_bank = arr<val<LOGN>,NUMG>{[&](u64 j) -> val<LOGN> {
            return readt[j] >> TAGHBITS;
        }};
        arr<val<1>,NUMG> owner_valid = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return bank_valid.select(owner_bank[j]);
        }};
        arr<val<1>,NUMG> bdir = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            val<1> tagged_dir = select(owner_valid[j], branch_taken_bank.select(owner_bank[j]), val<1>{0});
            return select(allocate[j], last_actual, tagged_dir);
        }};
        arr<val<1>,NUMG> badpred1 = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return readc[j] != bdir[j];
        }};
        arr<val<1>,NUMG> altdiffer = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return readc[j] != pred2_bank.select(owner_bank[j]);
        }};
        arr<val<1>,NUMG> goodpred = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return (owner_bank[j] != last_bank) | correct_pred;
        }};
        arr<val<1>,NUMG> g_weak = arr<val<1>,NUMG>{[&](u64 j) -> val<1> {
            return primary[j] & badpred1[j] & (readh[j] == hard<0>{});
        }};

        val<N> disagree_mask = (p1_rank.concat() ^ p2_rank.concat()) & branch_mask;
        arr<val<1>,N> primary_wrong = arr<val<1>,N>{[&](u64 r) -> val<1> {
            return pred1_rank[r] != branch_taken[r];
        }};
        arr<val<1>,N> bim_primary = arr<val<1>,N>{[&](u64 r) -> val<1> {
            return is_branch[r] & ~provider_hit_rank[r];
        }};
        arr<val<N>,N> bank_rmask = arr<val<N>,N>{[&](u64 b) -> val<N> {
            return arr<val<1>,N>{[&](u64 r) -> val<1> {
                return is_branch[r] & (bank_id_rank[r] == val<LOGN>{b});
            }}.concat();
        }};
        arr<val<1>,N> access_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return bank_rmask[b] != hard<0>{};
        }};
        arr<val<1>,N> disagree_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return (bank_rmask[b] & disagree_mask) != hard<0>{};
        }};
        arr<val<1>,N> p2_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return (bank_rmask[b] & p2_rank.concat()) != hard<0>{};
        }};
        arr<val<1>,N> taken_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return (bank_rmask[b] & branch_taken.concat()) != hard<0>{};
        }};
        arr<val<1>,N> bim_primary_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return (bank_rmask[b] & bim_primary.concat()) != hard<0>{};
        }};
        arr<val<1>,N> primary_wrong_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return (bank_rmask[b] & primary_wrong.concat()) != hard<0>{};
        }};
        arr<val<1>,N> p1_weak_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return execute_if(disagree_b[b], [&]() -> val<1> {
                return ~table1_hyst[b].read(index1);
            });
        }};
        arr<val<1>,N> b_weak_b = arr<val<1>,N>{[&](u64 b) -> val<1> {
            return execute_if(bim_primary_b[b] & primary_wrong_b[b], [&]() -> val<1> {
                return ~bhyst[b].read(bindex);
            });
        }};

        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
        val<1> extra_cycle = some_badpred1 | mispredict | (disagree_mask != hard<0>{});
        need_extra_cycle(extra_cycle);

#ifdef USE_META
        arr<val<1>,N> altdiff_rank = arr<val<1>,N>{[&](u64 r) -> val<1> {
            return alt_valid_rank[r] & (pred2_rank[r] != pred1_rank[r]);
        }};
        arr<val<2,i64>,N> meta_incr = arr<val<2,i64>,N>{[&](u64 r) -> val<2,i64> {
            val<1> update_meta = is_branch[r] & altdiff_rank[r] & newly_alloc_rank[r];
            val<1> bad_pred2 = pred2_rank[r] != branch_taken[r];
            return select(update_meta, concat(bad_pred2, val<1>{1}), val<2>{0});
        }};
        for (u64 i=METAPIPE-1; i!=0; i--) {
            meta[i] = meta[i-1];
        }
        auto newmeta = meta[0] + meta_incr.fo1().fold_add();
        using meta_t = valt<decltype(meta[0])>;
        meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                         select(newmeta < meta_t::minval, meta_t{meta_t::minval}, meta_t{newmeta}));
#endif

#ifdef CHEATING_MODE
#ifdef PERF_COUNTERS
        for (u64 j=0; j<NUMG; j++) {
            perf_table_reads[j] += num_branch;
            if (static_cast<bool>(allocate[j])) {
                perf_table_alloc[j]++;
            }
        }
        if (static_cast<bool>(extra_cycle)) perf_extra_cycle_total++;
        if (static_cast<bool>(some_badpred1)) perf_extra_cycle_badpred++;
        if (static_cast<bool>(mispredict)) perf_extra_cycle_mispredict++;
        if (static_cast<bool>(disagree_mask != hard<0>{})) perf_extra_cycle_p1_update++;
        if (static_cast<bool>(mispredict & noalloc)) {
            perf_alloc_failures++;
            if (static_cast<bool>(postmask == hard<0>{})) perf_alloc_fail_highest++;
            else perf_alloc_fail_noubit++;
        }
#ifdef USE_META
        for (u64 r=0; r<num_branch; r++) {
            bool update_meta = static_cast<bool>(is_branch[r] & alt_valid_rank[r] & (pred2_rank[r] != pred1_rank[r]) & newly_alloc_rank[r]);
            if (update_meta) {
                perf_meta_updates++;
                bool bad_pred2 = static_cast<bool>(pred2_rank[r] != branch_taken[r]);
                if (bad_pred2) perf_meta_bad_alt++;
                else perf_meta_good_alt++;
            }
        }
#endif
        for (u64 r=0; r<num_branch; r++) {
            bool actual = static_cast<bool>(branch_taken[r]);
            bool final_pred = static_cast<bool>(p2_rank[r]);
            perf_predictions++;
            if (final_pred == actual) perf_correct++;
            bool provider_hit = static_cast<bool>(provider_hit_rank[r]);
            bool alt_used = static_cast<bool>(alt_valid_rank[r] & (p2_rank[r] != pred1_rank[r]));
            if (alt_used) {
                perf_alt_used++;
                if (final_pred == actual) perf_alt_correct++;
                bool alt_provider_hit = static_cast<bool>(alt_provider_hit_rank[r]);
                if (alt_provider_hit) {
                    u64 aj = static_cast<u64>(alt_provider_idx_rank[r]);
                    if (aj < NUMG) {
                        perf_alt_provider_used[aj]++;
                        if (final_pred == actual) perf_alt_provider_correct[aj]++;
                    }
                } else {
                    perf_alt_bimodal_used++;
                    if (final_pred == actual) perf_alt_bimodal_correct++;
                }
            } else if (provider_hit) {
                u64 pj = static_cast<u64>(provider_idx_rank[r]);
                if (pj < NUMG) {
                    perf_provider_used[pj]++;
                    perf_table_hits[pj]++;
                    if (final_pred == actual) perf_provider_correct[pj]++;
                }
            } else {
                perf_bimodal_used++;
                if (final_pred == actual) perf_bimodal_correct++;
            }
        }
#endif
#endif

        // P1/Bim updates are issued on physical banks (unordered tables)
        for (u64 b=0; b<N; b++) {
            val<1> wr_p1pred = disagree_b[b] & p1_weak_b[b];
            val<1> p1pred_data = p2_b[b];
            execute_if(wr_p1pred, [&](){
                table1_pred[b].write(index1, p1pred_data);
            });

            val<1> wr_p1h = access_b[b];
            val<1> p1h_data = ~disagree_b[b];
            execute_if(wr_p1h, [&](){
                table1_hyst[b].write(index1, p1h_data);
            });

            val<1> wr_bim = bim_primary_b[b] & primary_wrong_b[b] & b_weak_b[b];
            val<1> bim_data = taken_b[b];
            execute_if(wr_bim, [&](){
                bim[b].write(bindex, bim_data);
            });

            val<1> wr_bh = bim_primary_b[b];
            val<1> bh_data = ~primary_wrong_b[b];
            execute_if(wr_bh, [&](){
                bhyst[b].write(bindex, bh_data);
            });
        }

        // P2 table updates: one write entry per RAM instance
        for (u64 j=0; j<NUMG; j++) {
            execute_if(allocate[j], [&](){
                gtag[j].write(gindex[j], concat(last_bank, htag[j]));
            });

            val<1> wr_pred = g_weak[j] | allocate[j];
            execute_if(wr_pred, [&](){
                val<1> new_pred = select(allocate[j], last_actual, bdir[j]);
                gpred[j].write(gindex[j], new_pred);
            });

            val<1> wr_hyst = primary[j] | allocate[j];
            execute_if(wr_hyst, [&](){
                val<2> new_h = select(allocate[j], val<2>{0}, update_ctr(readh[j], ~badpred1[j]));
                ghyst[j].write(gindex[j], new_h, extra_cycle);
            });

            val<1> update_u = primary[j] & altdiffer[j];
            val<1> wr_u = update_u | allocate[j] | uclear[j];
            execute_if(wr_u, [&](){
                val<1> new_u = goodpred[j] & ~allocate[j] & ~uclear[j];
                ubit[j].write(gindex[j], new_u, extra_cycle);
            });
        }

#ifdef RESET_UBITS
        val<NUMG> allocmask1 = collamask1.reverse();
        val<1> faralloc = (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
        val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
        uctr = select(correct_pred,
                      uctr,
                      select(uctrsat, val<decltype(uctr)::size>{0}, update_ctr(uctr, faralloc)));
        execute_if(uctrsat, [&](){
            for (auto &uram : ubit) uram.reset();
        });
#endif

        true_block = correct_pred | last_actual | last_pred() | line_end();
        execute_if(true_block, [&](){
            val<PATHBITS> path = next_pc >> 2;
            global_history1 = (global_history1 << 1) ^ val<GHIST1>{path};
            gfolds.update(path);
        });

        num_branch = 0;
    }
};
