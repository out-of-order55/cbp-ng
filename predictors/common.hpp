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
        val<1> has_buffered = write_bank != hard<0>{};
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
// write path supports WAW merge (same-address write updates buffered entry).
// one write request per cycle per instance.
// Each bank owns an independent pending-write queue of depth D.
template<u64 N, u64 M, u64 B, u64 D = 2, arith T = u64>
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
    static constexpr u64 S = B * D;

    static constexpr u64 slot_index(u64 bank_id, u64 slot_id)
    {
        return bank_id * D + slot_id;
    }

    // Native banked RAM backend (no rwram wrapper).
    ram<val<N,T>,E> bank[B];

    arr<reg<1>,S> q_valid;
    arr<reg<A>,S> q_addr;
    arr<reg<N,T>,S> q_data;
    arr<reg<Q>,B> q_head;
    arr<reg<Q>,B> q_tail;
    arr<reg<QC>,B> q_count;
    arr<reg<4>,B> q_head_wait;

    reg<1> write_stall;

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
        return data.fo1().fold_or();
    }

    void write(val<A> addr, val<N,T> data, val<1> we, val<1> noconflict)
    {
        write(addr, data, we, addr, val<1>{0}, noconflict);
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
        addr.fanout(hard<B+B*D>{});
        data.fanout(hard<B+B*D>{});
        val<1> wr_en = we.fo1();
        wr_en.fanout(hard<2 * B>{});
        val<1> ren_en = ren.fo1();
        val<1> no_conflict_force = noconflict.fo1();
        no_conflict_force.fanout(hard<B>{});
        write_stall = val<1>{0};
        q_head.fanout(hard<2>{});
        q_tail.fanout(hard<2>{});
        q_count.fanout(hard<3>{});
        q_valid.fanout(hard<4>{});
        q_addr.fanout(hard<5>{});
        q_data.fanout(hard<4>{});
        q_head_wait.fanout(hard<4>{});

        val<A> wr_addr_for_split = addr;
        val<A> rd_addr_for_split = read_addr;
        auto [wr_localaddr, wr_bankid] = split<L, I>(wr_addr_for_split.fo1());
        auto [rd_localaddr, rd_bankid] = split<L, I>(rd_addr_for_split.fo1());
        static_cast<void>(rd_localaddr);
        wr_localaddr.fanout(hard<B>{});
        arr<val<1>,B> wr_bank_split = wr_bankid.fo1().decode();
        arr<val<1>,B> rd_bank_split = rd_bankid.fo1().decode();
        // wr_bank_split.fanout(hard<D + 12>{});
        // rd_bank_split.fanout(hard<4>{});

        for (u64 b = 0; b < B; b++) {
            val<1> bank_sel = wr_bank_split[b];
            bank_sel.fanout(hard<2>{});
            val<1> bank_wr_en = wr_en & bank_sel;
            bank_wr_en.fanout(hard<5>{});
            val<1> bank_read_conflict = ren_en & rd_bank_split[b];

            val<Q> bank_head = q_head[b];
            val<Q> bank_tail = q_tail[b];
            val<QC> bank_count = q_count[b];
            bank_head.fanout(hard<3 + D>{});
            bank_tail.fanout(hard<1 + D>{});
            bank_count.fanout(hard<4>{});

            val<1> has_pending = bank_count != val<QC>{0};
            has_pending.fanout(hard<8>{});
            val<1> head_can_write = no_conflict_force | ~bank_read_conflict;
            head_can_write.fanout(hard<4>{});
            val<1> bank_pending_ready = bank_sel & has_pending & head_can_write;
            bank_pending_ready.fanout(hard<2>{});
            val<1> bank_active = bank_wr_en | bank_pending_ready;

            // Gate the heavy per-bank compare/select logic unless this bank is the target
            // of the current write or is the selected pending-drain bank with no conflict.
            execute_if(bank_active, [&](){
                arr<val<1>,D> bank_valid = [&](u64 s) -> val<1> {
                    return val<1>{q_valid[slot_index(b, s)]};
                };
                arr<val<A>,D> bank_addr = [&](u64 s) -> val<A> {
                    return val<A>{q_addr[slot_index(b, s)]};
                };
                bank_addr.fanout(hard<2>{});
                arr<val<N,T>,D> bank_data = [&](u64 s) -> val<N,T> {
                    return val<N,T>{q_data[slot_index(b, s)]};
                };

                val<1> full_now = bank_count == val<QC>{D};
                full_now.fanout(hard<2>{});

                val<D> hit_mask = execute_if(bank_wr_en & has_pending, [&](){
                    arr<val<1>,D> hit = [&](u64 s) -> val<1> {
                        return bank_valid[s] & (bank_addr[s] == addr);
                    };
                    return hit.concat();
                });
                hit_mask.fanout(hard<2>{});
                val<1> has_hit = hit_mask != val<D>{0};
                has_hit.fanout(hard<2>{});
                val<D> hit_oh = hit_mask.one_hot();
                arr<val<1>,D> hit_oh_split = hit_oh.make_array(val<1>{});
                // hit_oh_split.fanout(hard<2>{});

                val<A> head_addr = execute_if(has_pending, [&](){
                    return bank_addr.select(bank_head);
                });
                val<N,T> head_data = execute_if(has_pending, [&](){
                    return bank_data.select(bank_head);
                });
                val<L> head_localaddr = execute_if(has_pending, [&](){
                    auto [localaddr, head_bankid] = split<L, I>(head_addr);
                    static_cast<void>(head_bankid);
                    return localaddr;
                });


                val<1> direct_write = bank_wr_en & ~has_pending & head_can_write;
                direct_write.fanout(hard<4>{});

                val<1> do_drain_only = ~wr_en & bank_pending_ready;
                do_drain_only.fanout(hard<2>{});
                val<1> oldest_action = bank_wr_en & ~has_hit & full_now & has_pending;
                oldest_action.fanout(hard<3>{});
                val<1> do_flush = oldest_action & head_can_write;

                val<1> drop_oldest = oldest_action & ~head_can_write;

                val<1> pop_entry = do_drain_only | oldest_action;
                pop_entry.fanout(hard<D + 5>{});

                execute_if(do_drain_only | do_flush | direct_write, [&](){
                    val<L> waddr = select(direct_write, wr_localaddr, head_localaddr);
                    val<N,T> wdata = select(direct_write, data, head_data);
                    bank[b].write(waddr, wdata);
                });

                execute_if(has_pending, [&](){
                    val<4> bank_wait_prev = q_head_wait[b];
                    val<1> wait_active = ~pop_entry;
                    val<4> wait_inc = select(bank_wait_prev == hard<15>{},
                                             bank_wait_prev,
                                             val<4>{bank_wait_prev + hard<1>{}});
                    q_head_wait[b] = select(wait_active, wait_inc, val<4>{0});
                });

                val<QC> count_after_pop = select(pop_entry,
                                                 val<QC>{bank_count - val<QC>{1}},
                                                 bank_count);
                count_after_pop.fanout(hard<2>{});
                val<1> do_enqueue = bank_wr_en & ~has_hit & (~full_now | pop_entry) & ~direct_write;
                do_enqueue.fanout(hard<D + 3>{});
                val<QC> q_count_next = select(do_enqueue,
                                              val<QC>{count_after_pop + val<QC>{1}},
                                              count_after_pop);

                execute_if(pop_entry, [&](){
                    q_head[b] = val<Q>{bank_head + val<Q>{1}};
                });
                execute_if(do_enqueue, [&](){
                    q_tail[b] = val<Q>{bank_tail + val<Q>{1}};
                });
                execute_if(pop_entry | do_enqueue, [&](){
                    q_count[b] = q_count_next;
                });

                for (u64 s = 0; s < D; s++) {
                    const u64 idx = slot_index(b, s);
                    val<1> is_head = bank_head == val<Q>{s};
                    val<1> is_tail = bank_tail == val<Q>{s};
                    val<1> clear_slot = pop_entry & is_head;
                    val<1> fill_slot = do_enqueue & is_tail;
                    fill_slot.fanout(hard<4>{});
                    val<1> touch_meta = clear_slot | fill_slot;
                    val<1> touch_data = hit_oh_split[s] | fill_slot;

                    execute_if(touch_meta, [&](){
                        q_valid[idx] = select(fill_slot, val<1>{1}, val<1>{0});
                        q_addr[idx] = select(fill_slot, addr, val<A>{q_addr[idx]});
                    });
                    execute_if(touch_data, [&](){
                        q_data[idx] = data;
                    });
                }
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
            return val<QC>{q_count[i]} == val<QC>{D};
        };
        return bank_full.fo1().fold_or();
    }

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
