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

#ifdef PERF_COUNTERS
    u64 perf_rw_wr_req = 0;
    u64 perf_rw_wr_noconflict0 = 0;
    u64 perf_rw_buf_flush = 0;
    u64 perf_rw_buf_overwrite_drop = 0;
#endif

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
        val<1> has_buffered = write_bank != hard<0>{};
        val<1> buffered_done = (write_bank & (current_write | read_bank)) == hard<0>{};
#ifdef PERF_COUNTERS
        val<1> overwrite_drop = ~noconflict & has_buffered & ~buffered_done;
        val<1> flush_done = has_buffered & buffered_done;
        perf_rw_wr_req++;
        perf_rw_wr_noconflict0 += static_cast<u64>(~noconflict);
        perf_rw_buf_flush += static_cast<u64>(flush_done);
        perf_rw_buf_overwrite_drop += static_cast<u64>(overwrite_drop);
#endif
        execute_if(buffered_done.fo1() | ~noconflict, [&](){
            write_bank = banksel & ~noconflict_mask;
            execute_if(~noconflict,[&](){
                write_localaddr = localaddr;
                write_data = data;
            });
        });
    }

#ifdef PERF_COUNTERS
    u64 rw_wr_req() const { return perf_rw_wr_req; }
    u64 rw_wr_noconflict0() const { return perf_rw_wr_noconflict0; }
    u64 rw_buf_flush() const { return perf_rw_buf_flush; }
    u64 rw_buf_overwrite_drop() const { return perf_rw_buf_overwrite_drop; }
#endif

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
        q_valid.fanout(hard<D + 4>{});
        q_addr.fanout(hard<2*D + 12>{});
        q_data.fanout(hard<D + 6>{});
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
        hit_mask.fanout(hard<3>{});
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

        // Only build compare tree when there is an incoming write.
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
    static_assert(B >= 2 && B <= 64);
    static_assert(std::has_single_bit(B));
    static constexpr u64 E = M / B;
    static_assert(E > 1);
    static constexpr u64 L = std::bit_width(E - 1);
    static constexpr u64 I = std::bit_width(B - 1);
    static_assert(A == L + I);
    static_assert(D >= 2);
    static_assert(std::has_single_bit(D));
    static constexpr u64 Q = std::bit_width(D-1);
    static constexpr u64 QC = std::bit_width(D);

    // Native banked RAM backend (no rwram wrapper).
    ram<val<N,T>,E> bank[B];

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

    wb_rwram(const char *label="") : bank{label} {}

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
        auto [localaddr, bankid] = split<L, I>(addr.fo1());
        localaddr.fanout(hard<B>{});
        arr<val<1>,B> banksel = bankid.fo1().decode();
        banksel.fanout(hard<2>{});
        arr<val<N,T>,B> data = [&](u64 i) -> val<N,T> {
            return execute_if(banksel[i], [&](){ return bank[i].read(localaddr); });
        };
        return data.fo1().fold_or();
    }

    // Legacy API: keep backward compatibility.
    void write(val<A> addr, val<N,T> data, val<1> we, val<A> read_addr, val<1> ren)
    {
        write(addr, data, we, read_addr, ren, val<1>{0});
    }

    // noconflict is an incremental condition:
    //   - noconflict=1 guarantees no bank conflict this cycle
    //   - noconflict=0 falls back to read_addr/ren bank-conflict checks
    // Full queue policy:
    //   - if head can write RAM: flush oldest entry then enqueue new request
    //   - otherwise: drop oldest entry then enqueue new request
    void write(val<A> addr, val<N,T> data, val<1> we, val<A> read_addr, val<1> ren, val<1> noconflict)
    {
        addr.fanout(hard<D + 1>{});
        data.fanout(hard<D + 1>{});
        // read_addr.fanout(hard<2>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<5>{});
        val<1> ren_en = ren.fo1();
        val<1> no_conflict_force = noconflict.fo1();
        no_conflict_force.fanout(hard<2>{});
        q_head.fanout(hard<4 + D>{});
        q_tail.fanout(hard<2 + D>{});
        q_count.fanout(hard<4>{});
        q_valid.fanout(hard<2>{});
        q_addr.fanout(hard<3>{});
        q_data.fanout(hard<2>{});

        val<A> wr_addr_for_split = addr;
        val<A> rd_addr_for_split = read_addr;
        auto [wr_localaddr, wr_bankid] = split<L, I>(wr_addr_for_split.fo1());
        auto [rd_localaddr, rd_bankid] = split<L, I>(rd_addr_for_split.fo1());
        static_cast<void>(rd_localaddr);
        wr_localaddr.fanout(hard<D + 1>{});
        val<B> wr_bank_mask = wr_bankid.fo1().decode().concat();
        val<B> rd_bank_mask = rd_bankid.fo1().decode().concat();
        val<B> rd_bank_mask_gated = rd_bank_mask & ren_en.replicate(hard<B>{}).concat();
        wr_bank_mask.fanout(hard<D + 2>{});
        rd_bank_mask_gated.fanout(hard<2>{});
        arr<val<1>,B> wr_bank_split = wr_bank_mask.make_array(val<1>{});
        wr_bank_split.fanout(hard<2>{});

        val<1> has_pending = q_count != val<QC>{0};
        has_pending.fanout(hard<3>{});

        val<D> hit_mask = execute_if(wr_en, [&](){
            arr<val<1>,D> hit = arr<val<1>,D>{[&](u64 i){
                auto [q_local, q_bank] = split<L, I>(val<A>{q_addr[i]});
                val<B> q_bank_mask = q_bank.fo1().decode().concat();
                val<1> same_bank = (q_bank_mask & wr_bank_mask) != val<B>{0};
                val<1> same_local = q_local == wr_localaddr;
                return val<1>{q_valid[i]} & same_bank & same_local;
            }};
            return hit.concat();
        });
        val<1> has_hit = hit_mask != val<D>{0};
        has_hit.fanout(hard<2>{});
        val<D> hit_oh = hit_mask.one_hot();
        arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});

        val<1> full_now = q_count == val<QC>{D};
        full_now.fanout(hard<3>{});

        val<A> head_addr = q_addr.select(q_head);
        auto [head_localaddr, head_bankid] = split<L, I>(head_addr);

        val<B> head_bank_mask = head_bankid.fo1().decode().concat();
        head_bank_mask.fanout(hard<2>{});

        arr<val<1>,B> head_bank_split = head_bank_mask.make_array(val<1>{});

        val<1> head_conflict = (head_bank_mask & rd_bank_mask_gated) != val<B>{0};
        val<1> head_can_write = no_conflict_force | ~head_conflict;
        head_can_write.fanout(hard<2>{});

        // Full queue path:
        // 1) try writing oldest entry to RAM (if no conflict)
        // 2) overwrite oldest slot with incoming request.
        val<1> oldest_action = wr_en & ~has_hit & full_now & has_pending;
        oldest_action.fanout(hard<3>{});
        val<1> do_flush = oldest_action & head_can_write;
        do_flush.fanout(hard<1+B>{});
        val<1> drop_oldest = oldest_action & ~head_can_write;

        val<1> pop_oldest = oldest_action;
        pop_oldest.fanout(hard<D + 5>{});

        val<1> current_conflict = (wr_bank_mask & rd_bank_mask_gated) != val<B>{0};
        val<1> direct_can_write = no_conflict_force | ~current_conflict;

        val<1> direct_write = wr_en & ~has_pending & direct_can_write;
        direct_write.fanout(hard<B + 2>{});

        val<N,T> head_data = q_data.select(q_head);

        val<L> wr_local_sel = wr_localaddr;
        wr_local_sel.fanout(hard<B>{});
        val<L> head_local_sel = head_localaddr;
        head_local_sel.fanout(hard<B>{});
        val<N,T> data_sel = data;
        data_sel.fanout(hard<B>{});
        val<N,T> head_data_sel = head_data;
        head_data_sel.fanout(hard<B>{});
        execute_if(do_flush | direct_write, [&](){
            for (u64 i = 0; i < B; i++) {
                val<1> do_head = do_flush & head_bank_split[i];
                val<1> do_direct = direct_write & wr_bank_split[i];
                do_head.fanout(hard<3>{});
                execute_if(do_head | do_direct, [&](){
                    val<L> waddr = select(do_head, head_local_sel, wr_local_sel);
                    val<N,T> wdata = select(do_head, head_data_sel, data_sel);
                    bank[i].write(waddr, wdata);
                });
            }
        });
        q_head_wait.fanout(hard<4>{});
        val<4> q_head_wait_prev = q_head_wait;
        val<1> wait_active = has_pending & ~pop_oldest;
       
        val<4> wait_inc = select(q_head_wait == hard<15>{}, val<4>{q_head_wait}, val<4>{q_head_wait + hard<1>{}});
        q_head_wait = select(wait_active, wait_inc, val<4>{0});

        val<QC> count_after_pop = select(pop_oldest, val<QC>{q_count - val<QC>{1}}, val<QC>{q_count});
        count_after_pop.fanout(hard<2>{});
        val<1> do_enqueue = wr_en & ~has_hit & (~full_now | pop_oldest) & ~direct_write;
        do_enqueue.fanout(hard<D+2>{});
        val<1> stall_now = wr_en & ~has_hit & full_now & ~pop_oldest;
        write_stall = stall_now;

#ifdef PERF_COUNTERS
        perf_wb_wr_req += static_cast<u64>(wr_en);
        perf_wb_wr_hit_merge += static_cast<u64>(has_hit);
        perf_wb_wr_enq += static_cast<u64>(do_enqueue);
        // New semantics: on full queue, oldest is flushed if possible else dropped.
        perf_wb_wr_drop_full += static_cast<u64>(drop_oldest);
        perf_wb_drain_req += static_cast<u64>(oldest_action);
        perf_wb_drain_do += static_cast<u64>(do_flush);
        perf_wb_drain_block_wr += static_cast<u64>(drop_oldest);
        wb_record_head_lifetime<WB_LIFETIME_BINS>(static_cast<u64>(do_flush) != 0,
                                                  static_cast<u64>(q_head_wait_prev),
                                                  perf_wb_lifetime_samples,
                                                  perf_wb_lifetime_sum,
                                                  perf_wb_lifetime_max,
                                                  perf_wb_lifetime_hist);
#endif

        val<Q> q_head_next = select(pop_oldest, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
        val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
        val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_pop + val<QC>{1}}, count_after_pop);

#ifdef PERF_COUNTERS
        u64 depth_now = static_cast<u64>(q_count_next);
        perf_wb_depth_peak = std::max(perf_wb_depth_peak, depth_now);
#endif

        for (u64 i = 0; i < D; i++) {
            val<1> is_head = q_head == val<Q>{i};
            val<1> is_tail = q_tail == val<Q>{i};
            val<1> clear_slot = pop_oldest & is_head;
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
        for (u64 i = 0; i < B; i++) {
            bank[i].reset();
        }
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
