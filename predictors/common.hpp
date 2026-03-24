#ifndef COMMON_H
#define COMMON_H

#include "../harcom.hpp"

using namespace hcm;


// up-down saturating counter update function
template<u64 N, typename T>
[[nodiscard]] val<N,T> update_ctr(val<N,T> ctr, val<1> incr)
{
    ctr.fanout(hard<6>{});
    val<N,T> incsat = select(ctr==hard<ctr.maxval>{},ctr,val<N,T>{ctr+1});
    val<N,T> decsat = select(ctr==hard<ctr.minval>{},ctr,val<N,T>{ctr-1});
    return select(incr.fo1(),incsat.fo1(),decsat.fo1());
}

#ifdef PERF_COUNTERS
template<u64 LIFE_BINS>
inline void wb_record_head_lifetime(bool did_dequeue,
                                    u64 head_wait_cycles,
                                    u64 &samples,
                                    u64 &sum,
                                    u64 &maxv,
                                    std::array<u64,LIFE_BINS+1> &hist)
{
    if (!did_dequeue) return;
    // report only queue-head residence time until RAM write, in cycles.
    u64 lifetime = head_wait_cycles + 1;
    u64 bin = std::min(lifetime, LIFE_BINS);
    samples++;
    sum += lifetime;
    maxv = std::max(maxv, lifetime);
    hist[bin]++;
}
#endif


// banked RAM for doing (almost) 1 read-modify-write per cycle
template<u64 N, u64 M, u64 B>
struct rwram {
    // M entries of N bits (total), B banks
    static_assert(std::has_single_bit(M)); // number of entries is power of two
    static constexpr u64 A = std::bit_width(M-1); // address bits
    static_assert(B>=2 && B<=64);
    static_assert(std::has_single_bit(B)); // number of banks is power of two
    static constexpr u64 E = M/B; // entries per bank
    static_assert(E>1);
    static constexpr u64 L = std::bit_width(E-1); // local address bits
    static constexpr u64 I = std::bit_width(B-1); // bank ID bits
    static_assert(A==L+I);

    ram<val<N>,E> bank[B];
    reg<B> read_bank;

    // buffered write
    reg<B> write_bank;
    reg<L> write_localaddr;
    reg<N> write_data;

    rwram(const char *label="") : bank{label} {}

    val<N> read(val<A> addr)
    {
        auto [localaddr,bankid] = split<L,I>(addr.fo1());
        localaddr.fanout(hard<B>{});
        arr<val<1>,B> banksel = bankid.fo1().decode();
        banksel.fanout(hard<2>{});
        arr<val<N>,B> data = [&] (u64 i) -> val<N> {
            return execute_if(banksel[i],[&](){return bank[i].read(localaddr);});
        };
        read_bank = banksel.concat();
        return data.fo1().fold_or();
    }

    void write(val<A> addr, val<N> data, val<1> noconflict)
    {
        // if noconflict is set, there was no read in this cycle: we do the write immediately;
        // we do the buffered write if no conflict with the read or the current write
        auto [localaddr,bankid] = split<L,I>(addr.fo1());
        data.fanout(hard<B+1>{});
        noconflict.fanout(hard<B+2>{});
        val<B> banksel = bankid.fo1().decode().concat();
        banksel.fanout(hard<2>{});
        val<B> noconflict_mask = noconflict.replicate(hard<B>{}).concat();
        noconflict_mask.fanout(hard<2>{});
        val<B> current_write = banksel & noconflict_mask;
        current_write.fanout(hard<3>{});
        arr<val<1>,B> current_write_split = current_write.make_array(val<1>{});
        current_write_split.fanout(hard<3>{});
        arr<val<1>,B> write_bank_split = write_bank.make_array(val<1>{});
        arr<val<1>,B> read_bank_split = read_bank.make_array(val<1>{});
        for (u64 i=0; i<B; i++) {
            execute_if(current_write_split[i] | (write_bank_split[i].fo1() & ~read_bank_split[i].fo1()), [&](){
                val<L> a = select(current_write_split[i],localaddr,write_localaddr);
                val<N> d = select(current_write_split[i],data,write_data);
                bank[i].write(a.fo1(),d.fo1());
            });
        }
        // buffer the current write if not done
        // keep the previous write if not done and the current write is done
        // otherwise invalidate buffered write
        val<1> buffered_done = (write_bank & (current_write | read_bank)) == hard<0>{};
        execute_if(buffered_done.fo1() | ~noconflict, [&](){
            write_bank = banksel & ~noconflict_mask;
            execute_if(~noconflict,[&](){
                write_localaddr = localaddr;
                write_data = data;
            });
        });
    }

    void reset()
    {
        for (u64 i=0; i<B; i++) {
            bank[i].reset();
        }
    }
};


// write-buffered RAM wrapper
// read path always goes directly to backend RAM (no WB lookup).
// write path supports WAW merge (same-address write updates buffered entry).
// one write request per cycle per instance.
template<u64 N, u64 M, u64 D = 4, arith T = u64>
struct wb_ram {
    static_assert(std::has_single_bit(M));
    static constexpr u64 A = std::bit_width(M-1);
    static_assert(D >= 2);
    static_assert(std::has_single_bit(D));
    static constexpr u64 Q = std::bit_width(D-1);
    static constexpr u64 QC = std::bit_width(D);

    ram<val<N,T>,M> mem;

    arr<reg<1>,D> q_valid;
    arr<reg<A>,D> q_addr;
    arr<reg<N,T>,D> q_data;
    reg<Q> q_head;
    reg<Q> q_tail;
    reg<QC> q_count;
    reg<4> q_head_wait;

    reg<1> write_stall;

#ifdef PERF_COUNTERS
    static constexpr u64 WB_LIFETIME_BINS = 64;
    u64 perf_wb_rd_req = 0;
    u64 perf_wb_rd_pending_hit = 0;
    u64 perf_wb_wr_req = 0;
    u64 perf_wb_wr_hit_merge = 0;
    u64 perf_wb_wr_enq = 0;
    u64 perf_wb_wr_drop_full = 0;
    u64 perf_wb_drain_req = 0;
    u64 perf_wb_drain_do = 0;
    u64 perf_wb_drain_block_wr = 0;
    u64 perf_wb_depth_peak = 0;
    u64 perf_wb_block_cycles = 0;
    std::array<u64,D+1> perf_wb_block_depth = {};
    u64 perf_wb_lifetime_samples = 0;
    u64 perf_wb_lifetime_sum = 0;
    u64 perf_wb_lifetime_max = 0;
    std::array<u64,WB_LIFETIME_BINS+1> perf_wb_lifetime_hist = {};
#endif

    wb_ram(const char *label="") : mem{label} {}

    val<N,T> read(val<A> addr)
    {
#ifdef PERF_COUNTERS
        perf_wb_rd_req++;
        bool pending_hit = false;
        for (u64 i = 0; i < D; i++) {
            bool v = static_cast<bool>(val<1>{q_valid[i]});
            bool h = static_cast<bool>(val<A>{q_addr[i]} == addr);
            pending_hit |= (v & h);
        }
        perf_wb_rd_pending_hit += static_cast<u64>(pending_hit);
#endif
        return mem.read(addr);
    }

    // noconflict=1 means no backend read conflict this cycle.
    // drain is allowed whenever queue is non-empty and noconflict=1.
    void write(val<A> addr, val<N,T> data, val<1> we, val<1> noconflict)
    {
        addr.fanout(hard<2*D+2>{});
        data.fanout(hard<D+3>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<D+8>{});
        q_head.fanout(hard<5+D>{});
        q_tail.fanout(hard<2+D>{});
        q_count.fanout(hard<3>{});
        q_valid.fanout(hard<2>{});
        q_addr.fanout(hard<3>{});
        q_data.fanout(hard<2>{});
        val<1> can_drain = noconflict.fo1();
        can_drain.fanout(hard<8>{});
        val<1> has_pending = q_count != val<QC>{0};
        has_pending.fanout(hard<8>{});
        val<1> do_dequeue = has_pending & can_drain;
        do_dequeue.fanout(hard<D+12>{});
        val<1> direct_write = wr_en & ~has_pending & can_drain;
        direct_write.fanout(hard<4>{});
        val<4> q_head_wait_prev = q_head_wait;
        val<1> wait_active = has_pending & ~do_dequeue;
        q_head_wait.fanout(hard<4>{});
        val<4> wait_inc = select(q_head_wait == hard<15>{}, val<4>{q_head_wait}, val<4>{q_head_wait + hard<1>{}});
        q_head_wait = select(wait_active, wait_inc, val<4>{0});

        // only build compare tree when there is an incoming write
        val<D> hit_mask = execute_if(wr_en, [&](){
            arr<val<1>,D> hit = arr<val<1>,D>{[&](u64 i){
                return val<1>{q_valid[i]} & (val<A>{q_addr[i]} == addr);
            }};
            return hit.concat();
        });
        val<1> has_hit = hit_mask != val<D>{0};
        has_hit.fanout(hard<6>{});
        val<D> hit_oh = hit_mask.one_hot();
        arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});
        val<D> head_mask = val<Q>{q_head}.decode().concat();
        val<1> hit_head = (hit_mask & head_mask) != val<D>{0};

        // independent RAM-drain stage
        val<A> deq_addr = q_addr.select(q_head);
        val<N,T> deq_data_old = q_data.select(q_head);
        val<N,T> deq_data = select(hit_head, data, deq_data_old);
        val<A> ram_write_addr = select(do_dequeue, deq_addr, addr);
        val<N,T> ram_write_data = select(do_dequeue, deq_data, data);
        execute_if(do_dequeue | direct_write, [&](){
            mem.write(ram_write_addr, ram_write_data);
        });

        // queue update stage
        val<1> stall_now = execute_if(wr_en | do_dequeue, [&](){
            val<QC> count_after_deq = select(do_dequeue, val<QC>{q_count - val<QC>{1}}, val<QC>{q_count});
            count_after_deq.fanout(hard<3>{});
            val<1> full_after_deq = count_after_deq == val<QC>{D};
            count_after_deq.fanout(hard<2>{});
            val<1> do_enqueue = wr_en & ~has_hit & ~full_after_deq & ~direct_write;
            do_enqueue.fanout(hard<D+8>{});
            val<1> stall_inner = wr_en & ~has_hit & full_after_deq;

#ifdef PERF_COUNTERS
            perf_wb_wr_req += static_cast<u64>(wr_en);
            perf_wb_wr_hit_merge += static_cast<u64>(has_hit);
            perf_wb_wr_enq += static_cast<u64>(do_enqueue);
            perf_wb_wr_drop_full += static_cast<u64>(stall_inner);
            perf_wb_drain_req += static_cast<u64>(has_pending & can_drain);
            perf_wb_drain_do += static_cast<u64>(do_dequeue);
            perf_wb_drain_block_wr += static_cast<u64>(has_pending & ~can_drain);
            if (static_cast<u64>(stall_inner) != 0) {
                u64 depth_now = static_cast<u64>(q_count);
                depth_now = std::min(depth_now, D);
                perf_wb_block_cycles++;
                perf_wb_block_depth[depth_now]++;
            }
            wb_record_head_lifetime<WB_LIFETIME_BINS>(static_cast<u64>(do_dequeue) != 0,
                                                      static_cast<u64>(q_head_wait_prev),
                                                      perf_wb_lifetime_samples,
                                                      perf_wb_lifetime_sum,
                                                      perf_wb_lifetime_max,
                                                      perf_wb_lifetime_hist);
#endif

            val<Q> q_head_next = select(do_dequeue, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
            val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
            val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_deq + val<QC>{1}}, count_after_deq);

#ifdef PERF_COUNTERS
            u64 depth_now = static_cast<u64>(q_count_next);
            perf_wb_depth_peak = std::max(perf_wb_depth_peak, depth_now);
#endif

            for (u64 i = 0; i < D; i++) {
                val<1> is_head = q_head == val<Q>{i};
                val<1> is_tail = q_tail == val<Q>{i};
                val<1> clear_slot = do_dequeue & is_head;
                val<1> fill_slot = do_enqueue & is_tail;
                fill_slot.fanout(hard<3>{});

                q_valid[i] = select(fill_slot, val<1>{1},
                             select(clear_slot, val<1>{0}, val<1>{q_valid[i]}));
                q_addr[i] = select(fill_slot, addr, val<A>{q_addr[i]});
                q_data[i] = select(hit_oh_split[i] | fill_slot, data, val<N,T>{q_data[i]});
            }

            q_head = q_head_next;
            q_tail = q_tail_next;
            q_count = q_count_next;
            return stall_inner;
        });
        write_stall = stall_now;
    }

    val<1> stalled() const
    {
        return write_stall;
    }

    val<1> full() const
    {
        return q_count == val<QC>{D};
    }

    // Head request has waited >8 cycles without draining.
    val<1> drain_urgent()
    {
        return (q_count != val<QC>{0}) & (q_head_wait > hard<8>{});
    }

#ifdef PERF_COUNTERS
    u64 wb_rd_req() const { return perf_wb_rd_req; }
    u64 wb_rd_pending_hit() const { return perf_wb_rd_pending_hit; }
    u64 wb_wr_req() const { return perf_wb_wr_req; }
    u64 wb_wr_hit_merge() const { return perf_wb_wr_hit_merge; }
    u64 wb_wr_enq() const { return perf_wb_wr_enq; }
    u64 wb_wr_drop_full() const { return perf_wb_wr_drop_full; }
    u64 wb_drain_req() const { return perf_wb_drain_req; }
    u64 wb_drain_do() const { return perf_wb_drain_do; }
    u64 wb_drain_block_wr() const { return perf_wb_drain_block_wr; }
    u64 wb_depth_peak() const { return perf_wb_depth_peak; }
    u64 wb_block_cycles() const { return perf_wb_block_cycles; }
    u64 wb_block_depth(u64 d) const { return (d <= D) ? perf_wb_block_depth[d] : 0; }
    u64 wb_lifetime_samples() const { return perf_wb_lifetime_samples; }
    u64 wb_lifetime_sum() const { return perf_wb_lifetime_sum; }
    u64 wb_lifetime_max() const { return perf_wb_lifetime_max; }
    u64 wb_lifetime_bin(u64 b) const { return (b <= WB_LIFETIME_BINS) ? perf_wb_lifetime_hist[b] : 0; }
#endif

    void reset()
    {
        mem.reset();
    }
};

// write-buffered RAM wrapper with masked vector write.
// Storage is one aggregated row (L lanes, each W bits) per address.
// read path always goes directly to backend RAM (no WB lookup).
// write path supports WAW merge on same address with lane mask merge.
template<u64 W, u64 L, u64 M, u64 D = 4, arith T = i64>
struct wb_mask_ram {
    static_assert(std::has_single_bit(M));
    static constexpr u64 A = std::bit_width(M-1);
    static_assert(D >= 2);
    static_assert(std::has_single_bit(D));
    static_assert(L >= 1 && L <= 64);
    static constexpr u64 Q = std::bit_width(D-1);
    static constexpr u64 QC = std::bit_width(D);

    ram<arr<val<W,T>,L>,M> mem;

    arr<reg<1>,D> q_valid;
    arr<reg<A>,D> q_addr;
    arr<reg<L>,D> q_wmask;
    arr<reg<W,T>,D> q_data[L];
    reg<Q> q_head;
    reg<Q> q_tail;
    reg<QC> q_count;
    reg<4> q_head_wait;

    reg<1> write_stall;

#ifdef PERF_COUNTERS
    static constexpr u64 WB_LIFETIME_BINS = 64;
    u64 perf_wb_rd_req = 0;
    u64 perf_wb_rd_pending_hit = 0;
    u64 perf_wb_wr_req = 0;
    u64 perf_wb_wr_hit_merge = 0;
    u64 perf_wb_wr_enq = 0;
    u64 perf_wb_wr_drop_full = 0;
    u64 perf_wb_drain_req = 0;
    u64 perf_wb_drain_do = 0;
    u64 perf_wb_drain_block_wr = 0;
    u64 perf_wb_depth_peak = 0;
    u64 perf_wb_block_cycles = 0;
    std::array<u64,D+1> perf_wb_block_depth = {};
    u64 perf_wb_lifetime_samples = 0;
    u64 perf_wb_lifetime_sum = 0;
    u64 perf_wb_lifetime_max = 0;
    std::array<u64,WB_LIFETIME_BINS+1> perf_wb_lifetime_hist = {};
#endif

    wb_mask_ram(const char *label="") : mem{label} {}

    arr<val<W,T>,L> read(val<A> addr)
    {
#ifdef PERF_COUNTERS
        perf_wb_rd_req++;
        bool pending_hit = false;
        for (u64 i = 0; i < D; i++) {
            bool v = static_cast<bool>(val<1>{q_valid[i]});
            bool h = static_cast<bool>(val<A>{q_addr[i]} == addr);
            pending_hit |= (v & h);
        }
        perf_wb_rd_pending_hit += static_cast<u64>(pending_hit);
#endif
        return mem.read(addr);
    }

    // noconflict=1 means no backend read conflict this cycle.
    // wmask bit=1 means corresponding lane should be updated.
    // taken_mask bit controls up/down direction for update_ctr on each masked lane.
    void write(val<A> addr, arr<val<W,T>,L> data, val<L> wmask, val<L> taken_mask, val<1> we, val<1> noconflict)
    {
        addr.fanout(hard<2*D+2>{});
        data.fanout(hard<L+3>{});
        wmask.fanout(hard<D+4>{});
        val<1> wr_en = we.fo1() & (wmask != val<L>{0});
        wr_en.fanout(hard<D+8>{});
        q_head.fanout(hard<5+D>{});
        q_tail.fanout(hard<2+D>{});
        q_count.fanout(hard<3>{});
        q_valid.fanout(hard<2>{});
        q_addr.fanout(hard<3>{});
        q_wmask.fanout(hard<2>{});
        for (u64 lane = 0; lane < L; lane++) {
            q_data[lane].fanout(hard<2>{});
        }
        val<1> can_drain = noconflict.fo1();
        can_drain.fanout(hard<8>{});
        val<1> has_pending = q_count != val<QC>{0};
        has_pending.fanout(hard<8>{});
        val<1> do_dequeue = has_pending & can_drain;
        do_dequeue.fanout(hard<L+D+12>{});
        val<1> direct_write = wr_en & ~has_pending & can_drain;
        direct_write.fanout(hard<4>{});
        val<4> q_head_wait_prev = q_head_wait;
        val<1> wait_active = has_pending & ~do_dequeue;
        q_head_wait.fanout(hard<4>{});
        val<4> wait_inc = select(q_head_wait == hard<15>{}, val<4>{q_head_wait}, val<4>{q_head_wait + hard<1>{}});
        q_head_wait = select(wait_active, wait_inc, val<4>{0});

        arr<val<1>,D> hit = arr<val<1>,D>{[&](u64 i){
            return wr_en & val<1>{q_valid[i]} & (val<A>{q_addr[i]} == addr);
        }};
        val<D> hit_mask = hit.concat();
        val<1> has_hit = hit_mask != val<D>{0};
        has_hit.fanout(hard<6>{});
        val<D> hit_oh = hit_mask.one_hot();
        arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});
        val<D> head_mask = val<Q>{q_head}.decode().concat();
        val<1> hit_head = (hit_mask & head_mask) != val<D>{0};
        arr<val<1>,L> wmask_bits = wmask.make_array(val<1>{});
        arr<val<1>,L> taken_bits = taken_mask.make_array(val<1>{});
        wmask_bits.fanout(hard<D+3>{});
        taken_bits.fanout(hard<D+3>{});

        val<A> deq_addr = q_addr.select(q_head);
        arr<val<W,T>,L> deq_data_old = arr<val<W,T>,L>{[&](u64 lane){
            return q_data[lane].select(q_head);
        }};
        arr<val<W,T>,L> deq_data = arr<val<W,T>,L>{[&](u64 lane){
            val<W,T> deq_lane_new = update_ctr(deq_data_old[lane], taken_bits[lane]);
            return select(hit_head & wmask_bits[lane], deq_lane_new, deq_data_old[lane]);
        }};
        arr<val<W,T>,L> req_data_updated = arr<val<W,T>,L>{[&](u64 lane){
            val<W,T> req_lane_new = update_ctr(data[lane], taken_bits[lane]);
            return select(wmask_bits[lane], req_lane_new, data[lane]);
        }};
        val<A> ram_write_addr = select(do_dequeue, deq_addr, addr);
        arr<val<W,T>,L> ram_write_data = arr<val<W,T>,L>{[&](u64 lane){
            return select(do_dequeue, deq_data[lane], req_data_updated[lane]);
        }};
        execute_if(do_dequeue | direct_write, [&](){
            mem.write(ram_write_addr, ram_write_data);
        });

        val<1> stall_now = execute_if(wr_en | do_dequeue, [&](){
            val<QC> count_after_deq = select(do_dequeue, val<QC>{q_count - val<QC>{1}}, val<QC>{q_count});
            val<1> full_after_deq = count_after_deq == val<QC>{D};
            val<1> do_enqueue = wr_en & ~has_hit & ~full_after_deq & ~direct_write;
            do_enqueue.fanout(hard<D+8>{});
            val<1> stall_inner = wr_en & ~has_hit & full_after_deq;

#ifdef PERF_COUNTERS
            perf_wb_wr_req += static_cast<u64>(wr_en);
            perf_wb_wr_hit_merge += static_cast<u64>(has_hit);
            perf_wb_wr_enq += static_cast<u64>(do_enqueue);
            perf_wb_wr_drop_full += static_cast<u64>(stall_inner);
            perf_wb_drain_req += static_cast<u64>(has_pending & can_drain);
            perf_wb_drain_do += static_cast<u64>(do_dequeue);
            perf_wb_drain_block_wr += static_cast<u64>(has_pending & ~can_drain);
            if (static_cast<u64>(stall_inner) != 0) {
                u64 depth_now = static_cast<u64>(q_count);
                depth_now = std::min(depth_now, D);
                perf_wb_block_cycles++;
                perf_wb_block_depth[depth_now]++;
            }
            wb_record_head_lifetime<WB_LIFETIME_BINS>(static_cast<u64>(do_dequeue) != 0,
                                                      static_cast<u64>(q_head_wait_prev),
                                                      perf_wb_lifetime_samples,
                                                      perf_wb_lifetime_sum,
                                                      perf_wb_lifetime_max,
                                                      perf_wb_lifetime_hist);
#endif

            val<Q> q_head_next = select(do_dequeue, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
            val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
            val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_deq + val<QC>{1}}, count_after_deq);

#ifdef PERF_COUNTERS
            u64 depth_now = static_cast<u64>(q_count_next);
            perf_wb_depth_peak = std::max(perf_wb_depth_peak, depth_now);
#endif

            for (u64 i = 0; i < D; i++) {
                val<1> is_head = q_head == val<Q>{i};
                val<1> is_tail = q_tail == val<Q>{i};
                val<1> clear_slot = do_dequeue & is_head;
                val<1> fill_slot = do_enqueue & is_tail;
                fill_slot.fanout(hard<L+4>{});
                val<1> merge_slot = hit_oh_split[i];
                merge_slot.fanout(hard<L+2>{});

                q_valid[i] = select(fill_slot, val<1>{1},
                             select(clear_slot, val<1>{0}, val<1>{q_valid[i]}));
                q_addr[i] = select(fill_slot, addr, val<A>{q_addr[i]});
                q_wmask[i] = select(fill_slot, wmask,
                             select(merge_slot, val<L>{q_wmask[i]} | wmask, val<L>{q_wmask[i]}));
                for (u64 lane = 0; lane < L; lane++) {
                    val<1> lane_merge = merge_slot & wmask_bits[lane];
                    val<W,T> merged_lane = update_ctr(val<W,T>{q_data[lane][i]}, taken_bits[lane]);
                    q_data[lane][i] = select(fill_slot, req_data_updated[lane],
                                      select(lane_merge, merged_lane, val<W,T>{q_data[lane][i]}));
                }
            }

            q_head = q_head_next;
            q_tail = q_tail_next;
            q_count = q_count_next;
            return stall_inner;
        });
        write_stall = stall_now;
    }

    val<1> stalled() const
    {
        return write_stall;
    }

    val<1> full() const
    {
        return q_count == val<QC>{D};
    }

    val<1> drain_urgent()
    {
        return (q_count != val<QC>{0}) & (q_head_wait > hard<8>{});
    }

#ifdef PERF_COUNTERS
    u64 wb_rd_req() const { return perf_wb_rd_req; }
    u64 wb_rd_pending_hit() const { return perf_wb_rd_pending_hit; }
    u64 wb_wr_req() const { return perf_wb_wr_req; }
    u64 wb_wr_hit_merge() const { return perf_wb_wr_hit_merge; }
    u64 wb_wr_enq() const { return perf_wb_wr_enq; }
    u64 wb_wr_drop_full() const { return perf_wb_wr_drop_full; }
    u64 wb_drain_req() const { return perf_wb_drain_req; }
    u64 wb_drain_do() const { return perf_wb_drain_do; }
    u64 wb_drain_block_wr() const { return perf_wb_drain_block_wr; }
    u64 wb_depth_peak() const { return perf_wb_depth_peak; }
    u64 wb_block_cycles() const { return perf_wb_block_cycles; }
    u64 wb_block_depth(u64 d) const { return (d <= D) ? perf_wb_block_depth[d] : 0; }
    u64 wb_lifetime_samples() const { return perf_wb_lifetime_samples; }
    u64 wb_lifetime_sum() const { return perf_wb_lifetime_sum; }
    u64 wb_lifetime_max() const { return perf_wb_lifetime_max; }
    u64 wb_lifetime_bin(u64 b) const { return (b <= WB_LIFETIME_BINS) ? perf_wb_lifetime_hist[b] : 0; }
#endif

    void reset()
    {
        mem.reset();
    }
};


// write-buffered RWRAM wrapper
// read path always goes directly to backend RWRAM (no WB lookup).
// write path supports WAW merge (same-address write updates buffered entry).
// one write request per cycle per instance.
template<u64 N, u64 M, u64 B, u64 D = 4, arith T = u64>
struct wb_rwram {
    static_assert(std::has_single_bit(M));
    static constexpr u64 A = std::bit_width(M-1);
    static_assert(D >= 2);
    static_assert(std::has_single_bit(D));
    static constexpr u64 Q = std::bit_width(D-1);
    static constexpr u64 QC = std::bit_width(D);

    rwram<N,M,B> mem;

    arr<reg<1>,D> q_valid;
    arr<reg<A>,D> q_addr;
    arr<reg<N,T>,D> q_data;
    reg<Q> q_head;
    reg<Q> q_tail;
    reg<QC> q_count;
    reg<4> q_head_wait;

    reg<1> write_stall;

#ifdef PERF_COUNTERS
    static constexpr u64 WB_LIFETIME_BINS = 64;
    u64 perf_wb_rd_req = 0;
    u64 perf_wb_rd_pending_hit = 0;
    u64 perf_wb_wr_req = 0;
    u64 perf_wb_wr_hit_merge = 0;
    u64 perf_wb_wr_enq = 0;
    u64 perf_wb_wr_drop_full = 0;
    u64 perf_wb_drain_req = 0;
    u64 perf_wb_drain_do = 0;
    u64 perf_wb_drain_block_wr = 0;
    u64 perf_wb_depth_peak = 0;
    u64 perf_wb_lifetime_samples = 0;
    u64 perf_wb_lifetime_sum = 0;
    u64 perf_wb_lifetime_max = 0;
    std::array<u64,WB_LIFETIME_BINS+1> perf_wb_lifetime_hist = {};
#endif

    wb_rwram(const char *label="") : mem{label} {}

    val<N,T> read(val<A> addr)
    {
#ifdef PERF_COUNTERS
        perf_wb_rd_req++;
        bool pending_hit = false;
        for (u64 i = 0; i < D; i++) {
            bool v = static_cast<bool>(val<1>{q_valid[i]});
            bool h = static_cast<bool>(val<A>{q_addr[i]} == addr);
            pending_hit |= (v & h);
        }
        perf_wb_rd_pending_hit += static_cast<u64>(pending_hit);
#endif
        return val<N,T>{mem.read(addr)};
    }

    // noconflict is passed to backend rwram write when dequeuing.
    // drain is allowed whenever queue is non-empty and noconflict=1.
    void write(val<A> addr, val<N,T> data, val<1> we, val<1> noconflict)
    {
        addr.fanout(hard<2*D+2>{});
        data.fanout(hard<D+3>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<D+8>{});
        q_head.fanout(hard<5+D>{});
        q_tail.fanout(hard<2+D>{});
        q_count.fanout(hard<3>{});
        q_valid.fanout(hard<2>{});
        q_addr.fanout(hard<3>{});
        q_data.fanout(hard<2>{});
        val<1> can_drain = noconflict.fo1();
        can_drain.fanout(hard<8>{});

        val<1> has_pending = q_count != val<QC>{0};
        has_pending.fanout(hard<8>{});
        val<1> do_dequeue = has_pending & can_drain;
        do_dequeue.fanout(hard<D+10>{});
        val<1> direct_write = wr_en & ~has_pending & can_drain;
        direct_write.fanout(hard<4>{});
        val<4> q_head_wait_prev = q_head_wait;
        val<1> wait_active = has_pending & ~do_dequeue;
        q_head_wait.fanout(hard<4>{});
        val<4> wait_inc = select(q_head_wait == hard<15>{}, val<4>{q_head_wait}, val<4>{q_head_wait + hard<1>{}});
        q_head_wait = select(wait_active, wait_inc, val<4>{0});

        arr<val<1>,D> hit = arr<val<1>,D>{[&](u64 i){
            return wr_en & val<1>{q_valid[i]} & (val<A>{q_addr[i]} == addr);
        }};
        val<D> hit_mask = hit.concat();
        val<1> has_hit = hit_mask != val<D>{0};
        has_hit.fanout(hard<6>{});
        val<D> hit_oh = hit_mask.one_hot();
        arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});
        val<D> head_mask = val<Q>{q_head}.decode().concat();
        val<1> hit_head = (hit_mask & head_mask) != val<D>{0};

        val<A> deq_addr = q_addr.select(q_head);
        val<N,T> deq_data_old = q_data.select(q_head);
        val<N,T> deq_data = select(hit_head, data, deq_data_old);
        val<A> ram_write_addr = select(do_dequeue, deq_addr, addr);
        val<N,T> ram_write_data = select(do_dequeue, deq_data, data);
        execute_if(do_dequeue | direct_write, [&](){
            mem.write(ram_write_addr, val<N>{ram_write_data}, can_drain);
        });

        val<QC> count_after_deq = select(do_dequeue, val<QC>{q_count - val<QC>{1}}, val<QC>{q_count});
        val<1> full_after_deq = count_after_deq == val<QC>{D};
        val<1> do_enqueue = wr_en & ~has_hit & ~full_after_deq & ~direct_write;
        do_enqueue.fanout(hard<D+8>{});
        val<1> stall_now = wr_en & ~has_hit & full_after_deq;
        write_stall = stall_now;

#ifdef PERF_COUNTERS
        perf_wb_wr_req += static_cast<u64>(wr_en);
        perf_wb_wr_hit_merge += static_cast<u64>(has_hit);
        perf_wb_wr_enq += static_cast<u64>(do_enqueue);
        perf_wb_wr_drop_full += static_cast<u64>(stall_now);
        perf_wb_drain_req += static_cast<u64>(has_pending & can_drain);
        perf_wb_drain_do += static_cast<u64>(do_dequeue);
        perf_wb_drain_block_wr += static_cast<u64>(has_pending & ~can_drain);
        wb_record_head_lifetime<WB_LIFETIME_BINS>(static_cast<u64>(do_dequeue) != 0,
                                                  static_cast<u64>(q_head_wait_prev),
                                                  perf_wb_lifetime_samples,
                                                  perf_wb_lifetime_sum,
                                                  perf_wb_lifetime_max,
                                                  perf_wb_lifetime_hist);
#endif

        val<Q> q_head_next = select(do_dequeue, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
        val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
        val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_deq + val<QC>{1}}, count_after_deq);

#ifdef PERF_COUNTERS
        u64 depth_now = static_cast<u64>(q_count_next);
        perf_wb_depth_peak = std::max(perf_wb_depth_peak, depth_now);
#endif

        for (u64 i = 0; i < D; i++) {
            val<1> is_head = q_head == val<Q>{i};
            val<1> is_tail = q_tail == val<Q>{i};
            val<1> clear_slot = do_dequeue & is_head;
            val<1> fill_slot = do_enqueue & is_tail;

            q_valid[i] = select(fill_slot, val<1>{1},
                         select(clear_slot, val<1>{0}, val<1>{q_valid[i]}));
            q_addr[i] = select(fill_slot, addr, val<A>{q_addr[i]});
            q_data[i] = select(hit_oh_split[i] | fill_slot, data, val<N,T>{q_data[i]});
        }

        q_head = q_head_next;
        q_tail = q_tail_next;
        q_count = q_count_next;
    }

    val<1> stalled() const
    {
        return write_stall;
    }

    val<1> full() const
    {
        return q_count == val<QC>{D};
    }

#ifdef PERF_COUNTERS
    u64 wb_rd_req() const { return perf_wb_rd_req; }
    u64 wb_rd_pending_hit() const { return perf_wb_rd_pending_hit; }
    u64 wb_wr_req() const { return perf_wb_wr_req; }
    u64 wb_wr_hit_merge() const { return perf_wb_wr_hit_merge; }
    u64 wb_wr_enq() const { return perf_wb_wr_enq; }
    u64 wb_wr_drop_full() const { return perf_wb_wr_drop_full; }
    u64 wb_drain_req() const { return perf_wb_drain_req; }
    u64 wb_drain_do() const { return perf_wb_drain_do; }
    u64 wb_drain_block_wr() const { return perf_wb_drain_block_wr; }
    u64 wb_depth_peak() const { return perf_wb_depth_peak; }
    u64 wb_lifetime_samples() const { return perf_wb_lifetime_samples; }
    u64 wb_lifetime_sum() const { return perf_wb_lifetime_sum; }
    u64 wb_lifetime_max() const { return perf_wb_lifetime_max; }
    u64 wb_lifetime_bin(u64 b) const { return (b <= WB_LIFETIME_BINS) ? perf_wb_lifetime_hist[b] : 0; }
#endif

    void reset()
    {
        mem.reset();
    }
};


// banked RAM with parameterized buffered writes
// D controls the pending-write queue depth.
template<u64 N, u64 M, u64 B, u64 D = 2>
struct my_rwram {
    // M entries of N bits (total), B banks
    static_assert(std::has_single_bit(M)); // number of entries is power of two
    static constexpr u64 A = std::bit_width(M-1); // address bits
    static_assert(B>=2 && B<=64);
    static_assert(std::has_single_bit(B)); // number of banks is power of two
    static constexpr u64 E = M/B; // entries per bank
    static_assert(E>1);
    static constexpr u64 L = std::bit_width(E-1); // local address bits
    static constexpr u64 I = std::bit_width(B-1); // bank ID bits
    static_assert(A==L+I);

    static_assert(D>=2);
    static_assert(std::has_single_bit(D));
    static constexpr u64 Q = std::bit_width(D-1); // queue index bits

    ram<val<N>,E> bank[B];
    reg<B> read_bank;

    // pending-write queue
    arr<reg<1>,D> q_valid;
    arr<reg<B>,D> q_bank;
    arr<reg<L>,D> q_localaddr;
    arr<reg<N>,D> q_data;
    reg<Q> q_head;
    reg<Q> q_tail;
    reg<Q+1> q_count;

    my_rwram(const char *label="") : bank{label} {}

    val<N> read(val<A> addr)
    {
        auto [localaddr,bankid] = split<L,I>(addr.fo1());
        localaddr.fanout(hard<B>{});
        arr<val<1>,B> banksel = bankid.fo1().decode();
        banksel.fanout(hard<2>{});
        arr<val<N>,B> data = [&] (u64 i) -> val<N> {
            return execute_if(banksel[i],[&](){return bank[i].read(localaddr);});
        };
        read_bank = banksel.concat();
        return data.fo1().fold_or();
    }

    void write(val<A> addr, val<N> data, val<1> noconflict)
    {
        // noconflict=1: write-through for current write.
        // noconflict=0: queue current write and flush queued write if possible.
        auto [localaddr,bankid] = split<L,I>(addr.fo1());
        data.fanout(hard<B+2>{});
        noconflict.fanout(hard<B+D+6>{});

        val<B> banksel = bankid.fo1().decode().concat();
        banksel.fanout(hard<4>{});
        arr<val<1>,B> banksel_split = banksel.make_array(val<1>{});
        banksel_split.fanout(hard<2>{});
        arr<val<1>,B> read_bank_split = read_bank.make_array(val<1>{});
        read_bank_split.fanout(hard<2>{});

        val<1> has_pending = (q_count != hard<0>{});
        has_pending.fanout(hard<4>{});
        val<B> head_bank = q_bank.select(q_head);
        head_bank.fanout(hard<2>{});
        arr<val<1>,B> head_bank_split = head_bank.make_array(val<1>{});
        head_bank_split.fanout(hard<2>{});
        val<L> head_localaddr = q_localaddr.select(q_head);
        val<N> head_data = q_data.select(q_head);
        val<1> head_conflict = ((head_bank & read_bank) != hard<0>{});
        val<1> dequeue = has_pending & ~head_conflict;
        dequeue.fanout(hard<D+5>{});

        // Write current request (if conflict-free) or flush queue head (if possible).
        val<1> direct_write = noconflict;
        direct_write.fanout(hard<B+2>{});
        for (u64 i=0; i<B; i++) {
            val<1> do_direct = direct_write & banksel_split[i];
            val<1> do_flush = dequeue & head_bank_split[i];
            do_direct.fanout(hard<3>{});
            do_flush.fanout(hard<2>{});
            execute_if(do_direct | do_flush, [&](){
                val<L> waddr = select(do_direct, localaddr, head_localaddr);
                val<N> wdata = select(do_direct, data, head_data);
                bank[i].write(waddr.fo1(), wdata.fo1());
            });
        }

        // Queue enqueue policy: if full and no dequeue, current write is dropped.
        val<1> enqueue_req = ~noconflict;
        val<1> queue_full = (q_count == hard<D>{});
        val<1> enqueue = enqueue_req & (~queue_full | dequeue);
        enqueue.fanout(hard<D+4>{});

        // Queue metadata updates.
        val<Q> qhead_inc = q_head + hard<1>{};
        val<Q> qtail_inc = q_tail + hard<1>{};
        val<Q+1> count_after_deq = select(dequeue, val<Q+1>{q_count - hard<1>{}}, val<Q+1>{q_count});
        q_count = select(enqueue, val<Q+1>{count_after_deq + hard<1>{}}, count_after_deq);
        q_head = select(dequeue, qhead_inc, val<Q>{q_head});
        q_tail = select(enqueue, qtail_inc, val<Q>{q_tail});

        // Per-slot state updates.
        arr<val<1>,D> next_valid = [&](u64 i) -> val<1> {
            val<1> is_head = has_pending & (q_head == i);
            val<1> is_tail = enqueue & (q_tail == i);
            is_head.fanout(hard<2>{});
            is_tail.fanout(hard<2>{});
            val<1> clear_slot = dequeue & is_head & ~is_tail;
            return select(is_tail, val<1>{1},
                   select(clear_slot, val<1>{0}, val<1>{q_valid[i]}));
        };
        arr<val<B>,D> next_bank = [&](u64 i) -> val<B> {
            val<1> is_tail = enqueue & (q_tail == i);
            return select(is_tail, banksel, val<B>{q_bank[i]});
        };
        arr<val<L>,D> next_localaddr = [&](u64 i) -> val<L> {
            val<1> is_tail = enqueue & (q_tail == i);
            return select(is_tail, localaddr, val<L>{q_localaddr[i]});
        };
        arr<val<N>,D> next_data = [&](u64 i) -> val<N> {
            val<1> is_tail = enqueue & (q_tail == i);
            return select(is_tail, data, val<N>{q_data[i]});
        };

        q_valid = next_valid;
        q_bank = next_bank;
        q_localaddr = next_localaddr;
        q_data = next_data;
    }

    void reset()
    {
        for (u64 i=0; i<B; i++) {
            bank[i].reset();
        }
        // Keep reset behavior aligned with rwram: clear SRAM contents only.
        // Writing queue registers here can conflict with write() in the same cycle.
    }
};


// global history updated by shifting left by one bit and then xoring with some branch bits
// (branch direction, PC bits, target bits, whatever...)
// N is the history length in bits
template<u64 N>
struct global_history {

    arr<reg<1>,N> h; // initial value = 0 (consistent with folded history)

    void update(valtype auto in)
    {
        auto input = in.fo1().make_array(val<1>{});
        static_assert(input.size<=N);
        for (u64 i=N-1; i>=input.size; i--) h[i] = h[i-1];
        for (u64 i=input.size-1; i>=1; i--) h[i] = h[i-1] ^ input[i].fo1();
        h[0] = input[0].fo1();
    }

    val<1>& operator[] (u64 i)
    {
        return h[i];
    }

    void fanout(hardval auto fo)
    {
        h.fanout(fo);
    }
};


// folding a global history means splitting it into equal size chunks and XORing all chunks together
// folded_gh does folding incrementally with a circular shift register
// see P. Michaud, "A PPM-like, tag-based branch predictor", Journal of ILP, vol. 7, 2005
template<u64 F>
struct folded_gh {
    static_assert(F!=0);

    reg<F> folded; // initial value = 0 (consistent with global history)

    val<F> get()
    {
        return folded;
    }

    void fanout(hardval auto fo)
    {
        folded.fanout(fo);
    }

    template<u64 MAXL>
    void update(global_history<MAXL> &gh, hardval auto ghlen, valtype auto in)
    {
        // left shift of global history corresponds to left rotate of folded history
        // the bit that is pushed out of the global history is XORed out of the folded history
        constexpr u64 inbits = std::min(F,std::min(in.size,ghlen.value));
        val<inbits> input = in.fo1(); // truncate input if longer than global history
        auto f = folded.make_array(val<1>{});
        static_assert(f.size==F);
        val<1> outbit = gh[ghlen-1];
        u64 outpos = ghlen % F;
        arr<val<1>,F> ff = [&](u64 i){
            if (i==0) {
                return (outpos==0)? f[F-1].fo1()^outbit.fo1() : f[F-1].fo1();
            } else {
                return (outpos==i)? f[i-1].fo1()^outbit.fo1() : f[i-1].fo1();
            }
        };
        auto x = input.fo1().make_array(val<1>{});
        arr<val<1>,F> y = [&](u64 i){return (i<x.size)? x[i].fo1()^ff[i].fo1() : ff[i].fo1();};
        folded = y.fo1().concat();
    }
};


// geometrically increasing global history lengths folds
// NH = number of history lengths, MINH = shortest, MAXH = longest
// FOLDS = fold sizes (in bits)
template<u64 NH, u64 MINH, u64 MAXH, u64... FOLDS>
struct geometric_folds {
    static_assert(NH>=2);
    static constexpr u64 NF = sizeof...(FOLDS); // number of folds per history length

    static constexpr auto HLEN = [] () {
        std::array<u64,NH> hlen;
        u64 prevhl = 0;
        for (u64 i=0; i<NH; i++) {
            u64 hl = MINH * mypow(f64(MAXH)/MINH,f64(i)/(NH-1));
            hl = std::max(prevhl+1,hl);
            hlen[NH-1-i] = hl;
            prevhl = hl;
        }
        return hlen;
    }();

    static_assert(HLEN[0]==MAXH); // HLEN[0] is the longest history

    global_history<MAXH> gh;
    std::array<std::tuple<folded_gh<FOLDS>...>,NH> folds;

    template<u64 J=0>
    auto get(u64 i)
    {
        if (i>=NH) {
            std::cerr << "geometric folds: out of bound access\n";
            std::terminate();
        }
        return std::get<J>(folds[i]).get();
    }

    void fanout(hardval auto fo)
    {
        for (u64 i=0; i<NH; i++) {
            static_loop<NF>([&]<u64 J>(){
                std::get<J>(folds[i]).fanout(fo);
            });
        }
    }

    void update(valtype auto branchbits)
    {
        // update folds before global history
        branchbits.fanout(hard<NH*NF+1>{});
        gh.fanout(hard<std::max(u64(2),NF+1)>{});
        static_loop<NH>([&]<u64 I>(){
            static_loop<NF>([&]<u64 J>(){
                std::get<J>(folds[I]).update(gh,hard<HLEN[I]>{},branchbits);
            });
        });
        gh.update(branchbits);
    }
};


#endif // COMMON_H
