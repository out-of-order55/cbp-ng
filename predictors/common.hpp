#ifndef COMMON_H
#define COMMON_H

#include "../harcom.hpp"

#ifndef WB_RWRAM_DROP_OLDEST_ON_CONFLICT
#define WB_RWRAM_DROP_OLDEST_ON_CONFLICT 1
#endif

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
    u64 perf_reads = 0;
    u64 perf_writes = 0;
    u64 perf_stale_reads = 0;
    u64 perf_drop_writes = 0;
    u64 perf_write_merges = 0;
    u64 perf_shadow_data[M] = {};
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
#ifdef PERF_COUNTERS
        perf_reads++;
        u64 full_addr = static_cast<u64>(addr);
        u64 bank_mask = u64{1} << static_cast<u64>(bankid);
        bool pending_same_addr =
            (static_cast<u64>(write_bank) == bank_mask) &&
            (static_cast<u64>(write_localaddr) == static_cast<u64>(localaddr)) &&
            (static_cast<u64>(write_data) != perf_shadow_data[full_addr]);
        if (pending_same_addr)
            perf_stale_reads++;
#endif
        return data.fo1().fold_or();
    }

    void write(val<A> addr, val<N> data, val<1> noconflict)
    {
        write(addr, data, val<1>{1}, noconflict);
    }

    void write(val<A> addr, val<N> data, val<1> we, val<1> noconflict)
    {
        // if noconflict is set, there was no read in this cycle: we do the write immediately;
        // we do the buffered write if no conflict with the read or the current write
        auto [localaddr,bankid] = split<L,I>(addr.fo1());
        data.fanout(hard<B+1>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<B+3>{});
        noconflict.fanout(hard<B+2>{});
        val<B> banksel = bankid.fo1().decode().concat();
        banksel.fanout(hard<2>{});
        val<B> wr_en_mask = wr_en.replicate(hard<B>{}).concat();
        val<B> noconflict_mask = noconflict.replicate(hard<B>{}).concat();
        noconflict_mask.fanout(hard<2>{});
        val<B> current_write = banksel & wr_en_mask & noconflict_mask;
        current_write.fanout(hard<3>{});
        arr<val<1>,B> current_write_split = current_write.make_array(val<1>{});
        current_write_split.fanout(hard<3>{});
        arr<val<1>,B> write_bank_split = write_bank.make_array(val<1>{});
        arr<val<1>,B> read_bank_split = read_bank.make_array(val<1>{});
#ifdef PERF_COUNTERS
        bool wr_req = static_cast<bool>(we);
        bool no_conflict = static_cast<bool>(noconflict);
        u64 req_bank = static_cast<u64>(bankid);
        if (wr_req) {
            perf_writes++;
            bool dropped_pending = !no_conflict && (static_cast<u64>(write_bank & read_bank) != 0);
            if (dropped_pending)
                perf_drop_writes++;
        }
#endif
        for (u64 i=0; i<B; i++) {
            execute_if((current_write_split[i] | (wr_en & write_bank_split[i].fo1() & ~read_bank_split[i].fo1())), [&](){
                val<L> a = select(current_write_split[i],localaddr,write_localaddr);
                val<N> d = select(current_write_split[i],data,write_data);
                bank[i].write(a.fo1(),d.fo1());
            });
#ifdef PERF_COUNTERS
            if (wr_req) {
                bool do_bank_write =
                    (no_conflict && req_bank == i) ||
                    (((static_cast<u64>(write_bank) >> i) & 1) != 0 &&
                     ((static_cast<u64>(read_bank) >> i) & 1) == 0);
                if (do_bank_write) {
                    bool write_current = no_conflict && req_bank == i;
                    u64 full_addr = (i << L) | (write_current ? static_cast<u64>(localaddr)
                                                               : static_cast<u64>(write_localaddr));
                    perf_shadow_data[full_addr] = write_current ? static_cast<u64>(data)
                                                                 : static_cast<u64>(write_data);
                }
            }
#endif
        }
        // buffer the current write if not done
        // keep the previous write if not done and the current write is done
        // otherwise invalidate buffered write
        val<1> buffered_done = (write_bank & (current_write | read_bank)) == hard<0>{};
        execute_if(wr_en & (buffered_done.fo1() | ~noconflict), [&](){
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
#ifdef PERF_COUNTERS
        for (u64 i=0; i<M; i++) {
            perf_shadow_data[i] = u64{0};
        }
#endif
    }

#ifdef PERF_COUNTERS
    u64 perf_total_reads() const
    {
        return perf_reads;
    }

    u64 perf_total_writes() const
    {
        return perf_writes;
    }

    u64 perf_total_stale_reads() const
    {
        return perf_stale_reads;
    }

    u64 perf_total_drop_writes() const
    {
        return perf_drop_writes;
    }

    u64 perf_total_write_merges() const
    {
        return perf_write_merges;
    }

    u64 perf_total_write_updates() const
    {
        return perf_write_merges;
    }
#endif
};


// write-buffered RAM wrapper
// read path always goes directly to backend RAM (no WB lookup).
// write path supports WAW update (same-address write updates buffered entry).
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

    wb_ram(const char *label="") : mem{label} {}

    val<N,T> read(val<A> addr)
    {
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
        has_hit.fanout(hard<2>{});
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

            val<Q> q_head_next = select(do_dequeue, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
            val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
            val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_deq + val<QC>{1}}, count_after_deq);

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

    void reset()
    {
        mem.reset();
    }
};

// write-buffered RAM wrapper with masked vector write.
// Storage is one aggregated row (L lanes, each W bits) per address.
// read path always goes directly to backend RAM (no WB lookup).
// write path supports WAW update on same address; lane masks accumulate in queue.
template<u64 W, u64 L, u64 M, u64 D = 16, arith T = i64>
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

    wb_mask_ram(const char *label="") : mem{label} {}

    arr<val<W,T>,L> read(val<A> addr)
    {
        val<A> ram_addr = addr;
        return mem.read(ram_addr);
    }

    // noconflict=1 means no backend read conflict this cycle.
    // wmask bit=1 means corresponding lane should be updated.
    // taken_mask bit controls up/down direction for update_ctr on each masked lane.
    void write(val<A> addr, arr<val<W,T>,L> data, val<L> wmask, val<L> taken_mask, val<1> we, val<1> noconflict)
    {
        addr.fanout(hard<2*D+2>{});
        val<1> wr_en = we.fo1() & (wmask != val<L>{0});
        val<1> has_pending = q_count != val<QC>{0};
        val<1> can_drain = noconflict.fo1();
        val<1> pending_ready = has_pending & can_drain;
        pending_ready.fanout(hard<2>{});
        val<1> stall_now = execute_if(wr_en | pending_ready, [&](){
            data.fanout(hard<2>{});
            wmask.fanout(hard<D+3>{});
            wr_en.fanout(hard<5>{});
            q_head.fanout(hard<4+D+L>{});
            q_tail.fanout(hard<2+D>{});
            q_count.fanout(hard<3>{});
            q_valid.fanout(hard<2>{});
            q_addr.fanout(hard<3>{});
            q_wmask.fanout(hard<2>{});
            for (u64 lane = 0; lane < L; lane++) {
                q_data[lane].fanout(hard<2>{});
            }
            can_drain.fanout(hard<8>{});
            has_pending.fanout(hard<8>{});
            // pending_ready.fanout(hard<L+D+12>{});
            val<1> do_dequeue = pending_ready;
            do_dequeue.fanout(hard<7+D>{});
            val<1> direct_write = wr_en & ~has_pending & can_drain;
            direct_write.fanout(hard<2>{});
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
            hit_mask.fanout(hard<3>{});
            val<1> has_hit = hit_mask != val<D>{0};
            has_hit.fanout(hard<3>{});
            val<D> hit_oh = hit_mask.one_hot();
            arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});
            val<D> head_mask = val<Q>{q_head}.decode().concat();
            val<1> hit_head = (hit_mask & head_mask) != val<D>{0};
            arr<val<1>,L> wmask_bits = wmask.make_array(val<1>{});
            arr<val<1>,L> taken_bits = taken_mask.make_array(val<1>{});
            wmask_bits.fanout(hard<D+2>{});
            taken_bits.fanout(hard<D+2>{});

            val<A> deq_addr = q_addr.select(q_head);
            arr<val<W,T>,L> deq_data_old = arr<val<W,T>,L>{[&](u64 lane){
                return q_data[lane].select(q_head);
            }};
            deq_data_old.fanout(hard<2>{});
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

            return execute_if(wr_en | do_dequeue, [&](){
                val<QC> count_after_deq = select(do_dequeue, val<QC>{q_count - val<QC>{1}}, val<QC>{q_count});
                val<1> full_after_deq = count_after_deq == val<QC>{D};
                val<1> do_enqueue = wr_en & ~has_hit & ~full_after_deq & ~direct_write;
                do_enqueue.fanout(hard<D+2>{});
                val<1> stall_inner = wr_en & ~has_hit & full_after_deq;

                val<Q> q_head_next = select(do_dequeue, val<Q>{q_head + val<Q>{1}}, val<Q>{q_head});
                val<Q> q_tail_next = select(do_enqueue, val<Q>{q_tail + val<Q>{1}}, val<Q>{q_tail});
                val<QC> q_count_next = select(do_enqueue, val<QC>{count_after_deq + val<QC>{1}}, count_after_deq);

                for (u64 i = 0; i < D; i++) {
                    val<1> is_head = q_head == val<Q>{i};
                    val<1> is_tail = q_tail == val<Q>{i};
                    val<1> clear_slot = do_dequeue & is_head;
                    val<1> fill_slot = do_enqueue & is_tail;
                    fill_slot.fanout(hard<L+3>{});
                    val<1> merge_slot = hit_oh_split[i];
                    merge_slot.fanout(hard<L+1>{});
                    val<1> valid_next = select(fill_slot, val<1>{1},
                                        select(clear_slot, val<1>{0}, val<1>{q_valid[i]}));
                    val<A> addr_next = select(fill_slot, addr, val<A>{q_addr[i]});
                    val<L> wmask_next = select(fill_slot, wmask,
                                        select(merge_slot, val<L>{q_wmask[i]} | wmask, val<L>{q_wmask[i]}));
                    q_valid[i] = valid_next;
                    q_addr[i] = addr_next;
                    q_wmask[i] = wmask_next;
                    for (u64 lane = 0; lane < L; lane++) {
                        val<1> lane_merge = merge_slot & wmask_bits[lane];
                        val<W,T> merged_lane = update_ctr(val<W,T>{q_data[lane][i]}, taken_bits[lane]);
                        val<W,T> lane_next = select(fill_slot, req_data_updated[lane],
                                             select(lane_merge, merged_lane, val<W,T>{q_data[lane][i]}));
                        q_data[lane][i] = lane_next;
                    }
                }

                q_head = q_head_next;
                q_tail = q_tail_next;
                q_count = q_count_next;
                return stall_inner;
            });
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

    void reset()
    {
        mem.reset();
    }
};


// write-buffered RWRAM wrapper
// read path always goes directly to backend RWRAM (no WB lookup).
// write path supports WAW update (same-address write updates buffered entry).
// one write request per cycle per instance.
// Each bank owns one pending-write slot.
template<u64 N, u64 M, u64 B, arith T = u64>
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
    static constexpr u64 DEPTH = 1;

    // Native banked RAM backend (no rwram wrapper).
    ram<val<N,T>,E> bank[B];

    arr<reg<1>,B> q_valid;
    arr<reg<A>,B> q_addr;
    arr<reg<N,T>,B> q_data;
    arr<reg<4>,B> pending_wait;

    reg<1> write_stall;

#ifdef PERF_COUNTERS
    u64 perf_reads = 0;
    u64 perf_writes = 0;
    u64 perf_stale_reads = 0;
    u64 perf_drop_writes = 0;
    u64 perf_write_merges = 0;
    u64 perf_read_pending_hits = 0;
    u64 perf_read_pending_hit_depth[DEPTH] = {};
    T perf_shadow_data[M] = {};
#endif

    wb_rwram(const char *label="") : bank{label} {}

    val<N,T> read(val<A> addr)
    {
        auto [localaddr, bankid] = split<L, I>(addr.fo1());
        localaddr.fanout(hard<B>{});
        arr<val<1>,B> banksel = bankid.fo1().decode();
        banksel.fanout(hard<2>{});
        arr<val<N,T>,B> data = [&](u64 i) -> val<N,T> {
            return execute_if(banksel[i], [&](){ return bank[i].read(localaddr); });
        };
#ifdef PERF_COUNTERS
        perf_reads++;
        u64 full_addr = static_cast<u64>(addr);
        u64 bank_idx = static_cast<u64>(bankid);
        bool stale = false;
        if (static_cast<bool>(q_valid[bank_idx]) &&
            static_cast<u64>(q_addr[bank_idx]) == full_addr) {
            perf_read_pending_hits++;
            perf_read_pending_hit_depth[0]++;
            stale = static_cast<T>(q_data[bank_idx]) != perf_shadow_data[full_addr];
        }
        if (stale)
            perf_stale_reads++;
#endif
        return data.fo1().fold_or();
    }

    void write(val<A> addr, val<N,T> data, val<1> we, val<1> noconflict)
    {
        write_impl(addr, data, we, addr, val<1>{0}, noconflict, val<1>{0}, val<1>{0});
    }

    // Legacy API: keep backward compatibility.
    void write(val<A> addr, val<N,T> data, val<1> we, val<A> read_addr, val<1> ren)
    {
        write_impl(addr, data, we, read_addr, ren, val<1>{0}, val<1>{0}, val<1>{0});
    }

    // noconflict is an incremental condition:
    //   - noconflict=1 guarantees no bank conflict this cycle
    //   - noconflict=0 falls back to read_addr/ren bank-conflict checks
    // Single-slot policy:
    //   - if the pending slot can drain, write it back first and optionally refill it
    //   - otherwise, keep the pending slot unless DROP_OLDEST is enabled
    void write(val<A> addr, val<N,T> data, val<1> we, val<A> read_addr, val<1> ren, val<1> noconflict)
    {
        write_impl(addr, data, we, read_addr, ren, noconflict, val<1>{0}, val<1>{0});
    }

    void write_update(val<A> addr, val<N,T> data, val<1> we, val<1> update_hit, val<1> update_incr,
                      val<A> read_addr, val<1> ren, val<1> noconflict)
    {
        write_impl(addr, data, we, read_addr, ren, noconflict, update_hit, update_incr);
    }

    void write_impl(val<A> addr, val<N,T> data, val<1> we, val<A> read_addr, val<1> ren,
                    val<1> noconflict, val<1> update_hit, val<1> update_incr)
    {
        addr.fanout(hard<2 * B>{});
        data.fanout(hard<2 * B>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<2 * B>{});
        val<1> update_hit_en = update_hit.fo1() & wr_en;
        update_hit_en.fanout(hard<2 * B>{});
        val<1> update_incr_en = update_incr.fo1();
        update_incr_en.fanout(hard<2 * B>{});
        val<1> ren_en = ren.fo1();
        val<1> no_conflict_force = noconflict.fo1();
        no_conflict_force.fanout(hard<B>{});
#ifdef PERF_COUNTERS
        if (static_cast<bool>(wr_en)) {
            perf_writes++;
        }
#endif
        write_stall = val<1>{0};
        q_valid.fanout(hard<4>{});
        q_addr.fanout(hard<5>{});
        q_data.fanout(hard<4>{});
        pending_wait.fanout(hard<4>{});

        val<A> wr_addr_for_split = addr;
        val<A> rd_addr_for_split = read_addr;
        auto [wr_localaddr, wr_bankid] = split<L, I>(wr_addr_for_split.fo1());
        auto [rd_localaddr, rd_bankid] = split<L, I>(rd_addr_for_split.fo1());
        static_cast<void>(rd_localaddr);
        wr_localaddr.fanout(hard<B>{});
        arr<val<1>,B> wr_bank_split = wr_bankid.fo1().decode();
        arr<val<1>,B> rd_bank_split = rd_bankid.fo1().decode();
        // rd_bank_split.fanout(hard<4>{});

        for (u64 b = 0; b < B; b++) {
            val<1> bank_sel = wr_bank_split[b];
            bank_sel.fanout(hard<2>{});
            val<1> bank_wr_en = wr_en & bank_sel;
            bank_wr_en.fanout(hard<5>{});
            val<1> bank_update_hit = update_hit_en & bank_sel;
            bank_update_hit.fanout(hard<4>{});
            val<1> bank_read_conflict = ren_en & rd_bank_split[b];

            val<1> has_pending = q_valid[b];
            has_pending.fanout(hard<8>{});
            val<1> head_can_write = no_conflict_force | ~bank_read_conflict;
            head_can_write.fanout(hard<4>{});
            val<1> bank_pending_ready = has_pending & head_can_write;
            bank_pending_ready.fanout(hard<2>{});
            val<1> bank_active = bank_wr_en | bank_pending_ready;

            // Gate the per-bank logic unless this bank is the target
            // of the current write or has a pending slot that can drain this cycle.
            execute_if(bank_active, [&](){
                val<1> slot_valid = val<1>{q_valid[b]};
                val<A> slot_addr = val<A>{q_addr[b]};
                val<N,T> slot_data = val<N,T>{q_data[b]};
                slot_addr.fanout(hard<2>{});

                val<1> has_hit = execute_if(bank_wr_en & has_pending, [&](){
                    return slot_valid & (slot_addr == addr);
                });
                has_hit.fanout(hard<2>{});
                val<N,T> head_data = slot_data;
                val<N,T> resolved_hit_data = select(bank_update_hit,
                                                    update_ctr(slot_data, update_incr_en),
                                                    data);
                resolved_hit_data.fanout(hard<3>{});
                val<L> head_localaddr = execute_if(has_pending, [&](){
                    auto [localaddr, head_bankid] = split<L, I>(slot_addr);
                    static_cast<void>(head_bankid);
                    return localaddr;
                });
                val<1> direct_write = bank_wr_en & ~has_pending & head_can_write;
                direct_write.fanout(hard<4>{});
                val<1> hit_head = has_hit;
                hit_head.fanout(hard<2>{});
                val<1> do_drain_ready = bank_pending_ready;
                do_drain_ready.fanout(hard<2>{});
                val<1> full_conflict = bank_wr_en & ~has_hit & has_pending;
                full_conflict.fanout(hard<3>{});
                val<1> blocked_full = full_conflict & ~head_can_write;
                blocked_full.fanout(hard<2>{});
#if WB_RWRAM_DROP_OLDEST_ON_CONFLICT
                val<1> pop_entry = do_drain_ready | blocked_full;
#else
                val<1> pop_entry = do_drain_ready;
#endif
                pop_entry.fanout(hard<6>{});

                execute_if(do_drain_ready | direct_write, [&](){
                    val<L> waddr = select(direct_write, wr_localaddr, head_localaddr);
                    val<N,T> pending_wdata = select(hit_head, resolved_hit_data, head_data);
                    val<N,T> wdata = select(direct_write, data, pending_wdata);
                    bank[b].write(waddr, wdata);
                });
#ifdef PERF_COUNTERS
                if (static_cast<bool>(bank_wr_en) || static_cast<bool>(bank_pending_ready)) {
                    if (static_cast<bool>(bank_wr_en) && static_cast<bool>(has_hit))
                        perf_write_merges++;
                    if (static_cast<bool>(blocked_full))
                        perf_drop_writes++;
                    bool do_backend_write =
                        static_cast<bool>(do_drain_ready) ||
                        static_cast<bool>(direct_write);
                    if (do_backend_write) {
                        bool write_direct = static_cast<bool>(direct_write);
                        bool write_updated_head = static_cast<bool>(do_drain_ready) && static_cast<bool>(hit_head);
                        T head_write_data = static_cast<T>(resolved_hit_data);
                        u64 full_addr = (b << L) | (write_direct ? static_cast<u64>(wr_localaddr)
                                                                  : static_cast<u64>(head_localaddr));
                        perf_shadow_data[full_addr] = write_direct ? static_cast<T>(data)
                                                                    : (write_updated_head ? head_write_data
                                                                                          : static_cast<T>(head_data));
                    }
                }
#endif

                execute_if(has_pending, [&](){
                    val<4> bank_wait_prev = pending_wait[b];
                    val<1> wait_active = ~pop_entry;
                    val<4> wait_inc = select(bank_wait_prev == hard<15>{},
                                             bank_wait_prev,
                                             val<4>{bank_wait_prev + hard<1>{}});
                    pending_wait[b] = select(wait_active, wait_inc, val<4>{0});
                });

                val<1> do_enqueue = bank_wr_en & ~has_hit & (~has_pending | pop_entry) & ~direct_write;
                do_enqueue.fanout(hard<4>{});
                val<1> clear_slot = pop_entry;
                val<1> fill_slot = do_enqueue;
                fill_slot.fanout(hard<4>{});
                val<1> update_slot = has_hit;
                update_slot.fanout(hard<2>{});
                val<1> touch_meta = clear_slot | fill_slot;
                val<1> touch_data = update_slot | fill_slot;

                execute_if(touch_meta, [&](){
                    q_valid[b] = select(fill_slot, val<1>{1}, val<1>{0});
                    q_addr[b] = select(fill_slot, addr, val<A>{q_addr[b]});
                });
                execute_if(touch_data, [&](){
                    q_data[b] = select(fill_slot, data, resolved_hit_data);
                });
            });
        }
    }

    val<1> stalled() const
    {
        return write_stall;
    }

    val<1> full() const
    {
        arr<val<1>,B> bank_full = [&](u64 i) -> val<1> {
            return val<1>{q_valid[i]};
        };
        return bank_full.fo1().fold_or();
    }

    void reset()
    {
        for (u64 i = 0; i < B; i++) {
            bank[i].reset();
        }
#ifdef PERF_COUNTERS
        for (u64 i = 0; i < M; i++) {
            perf_shadow_data[i] = T{};
        }
#endif
    }

#ifdef PERF_COUNTERS
    u64 perf_total_reads() const
    {
        return perf_reads;
    }

    u64 perf_total_writes() const
    {
        return perf_writes;
    }

    u64 perf_total_stale_reads() const
    {
        return perf_stale_reads;
    }

    u64 perf_total_drop_writes() const
    {
        return perf_drop_writes;
    }

    u64 perf_total_write_merges() const
    {
        return perf_write_merges;
    }

    u64 perf_total_write_updates() const
    {
        return perf_write_merges;
    }

    u64 perf_total_read_pending_hits() const
    {
        return perf_read_pending_hits;
    }

    u64 perf_total_read_pending_hit_depth(u64 depth) const
    {
        return depth == 0 ? perf_read_pending_hit_depth[0] : 0;
    }
#endif
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
