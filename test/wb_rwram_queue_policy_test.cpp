#include "harcom.hpp"
#include "predictors/common.hpp"

#include <iostream>
#include <string>

using namespace hcm;

class harcom_superuser {
public:
    void next_cycle() { panel.next_cycle(); }
    auto get(valtype auto x) { return x.get(); }
} hsu;

struct TestSummary {
    u64 passed = 0;
    u64 failed = 0;

    void expect(bool cond, const std::string &name)
    {
        if (cond) {
            passed++;
            std::cout << "[PASS] " << name << "\n";
        } else {
            failed++;
            std::cout << "[FAIL] " << name << "\n";
        }
    }
};

using WBRW = wb_rwram<8,8,2,4>;
static constexpr u64 QDEPTH = 4;

static constexpr u64 ADDR_SPACE = (1ull << WBRW::A);

template<u64 DEPTH>
static val<WBRW::A> delay_addr(val<WBRW::A> x)
{
    if constexpr (DEPTH == 0) {
        return x;
    } else {
        return delay_addr<DEPTH - 1>(x ^ val<WBRW::A>{0});
    }
}

template<u64 N, u64 M, u64 B, u64 D>
static void clear_wb_state(wb_rwram<N,M,B,D> &wb)
{
    for (u64 i = 0; i < D; i++) {
        wb.q_valid[i] = val<1>{0};
        wb.q_addr[i] = val<wb_rwram<N,M,B,D>::A>{0};
        wb.q_data[i] = val<N>{0};
    }
    wb.q_head = val<wb_rwram<N,M,B,D>::Q>{0};
    wb.q_tail = val<wb_rwram<N,M,B,D>::Q>{0};
    wb.q_count = val<wb_rwram<N,M,B,D>::QC>{0};
    wb.q_head_wait = val<4>{0};
    wb.write_stall = val<1>{0};
}

template<typename WB>
static u64 qcount_of(WB &wb)
{
    return static_cast<u64>(hsu.get(val<WB::QC>{wb.q_count}));
}

static u64 read_mem(WBRW &wb, u64 addr)
{
    val<WBRW::A> probe_addr = delay_addr<512>(val<WBRW::A>{addr});
    return static_cast<u64>(hsu.get(wb.read(probe_addr)));
}

static u64 bank_of(u64 addr)
{
    auto [local, bank] = split<WBRW::L, WBRW::I>(val<WBRW::A>{addr});
    [[maybe_unused]] auto keep_local = local;
    return static_cast<u64>(hsu.get(bank));
}

static u64 find_addr_with_bank(u64 bank, bool equal)
{
    for (u64 a = 0; a < ADDR_SPACE; a++) {
        bool cond = (bank_of(a) == bank);
        if (cond == equal) {
            return a;
        }
    }
    return 0;
}

static u64 head_addr_of(WBRW &wb)
{
    return static_cast<u64>(hsu.get(wb.q_addr.select(val<WBRW::Q>{wb.q_head})));
}

static void reset_case(WBRW &wb)
{
    wb.reset();
    clear_wb_state(wb);
    hsu.next_cycle();
}

static void fill_full_with_enqueued_writes(WBRW &wb, u64 base_data)
{
    // Use read_addr=addr and ren=1 to force conflict so writes enqueue.
    for (u64 i = 0; i < QDEPTH; i++) {
        wb.write(val<WBRW::A>{i}, val<8>{base_data + i}, val<1>{1},
                 val<WBRW::A>{i}, val<1>{1}, val<1>{0});
        hsu.next_cycle();
    }
}

int main()
{
    TestSummary t;
    WBRW wb{"wb_rwram_queue_policy"};
    panel.make_floorplan();

    // Case 1: direct write when queue empty and no read-bank conflict.
    reset_case(wb);
    wb.write(val<WBRW::A>{1}, val<8>{11}, val<1>{1}, val<WBRW::A>{0}, val<1>{0}, val<1>{0});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == 0, "direct write keeps queue empty");
    t.expect(read_mem(wb, 1) == 11, "direct write commits to RAM");

    // Case 2: queue empty but read-bank conflict blocks direct write -> enqueue.
    reset_case(wb);
    wb.write(val<WBRW::A>{2}, val<8>{22}, val<1>{1}, val<WBRW::A>{2}, val<1>{1}, val<1>{0});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == 1, "read-bank conflict enqueues write");
    t.expect(read_mem(wb, 2) == 0, "enqueued write not yet in RAM");

    // Case 2b: queue empty + read-bank conflict + noconflict=1 => force direct write.
    reset_case(wb);
    wb.write(val<WBRW::A>{3}, val<8>{33}, val<1>{1}, val<WBRW::A>{3}, val<1>{1}, val<1>{1});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == 0, "noconflict forces direct write under conflict");
    t.expect(read_mem(wb, 3) == 33, "forced direct write commits to RAM");

    // Case 3: full queue + non-conflicting read bank => flush oldest then overwrite oldest slot.
    reset_case(wb);
    fill_full_with_enqueued_writes(wb, 10);
    t.expect(qcount_of(wb) == QDEPTH, "queue reaches full depth");
    u64 old_head_addr = head_addr_of(wb);
    u64 old_head_bank = bank_of(old_head_addr);
    u64 rd_nonconf = find_addr_with_bank(old_head_bank, false);
    wb.write(val<WBRW::A>{6}, val<8>{99}, val<1>{1}, val<WBRW::A>{rd_nonconf}, val<1>{1}, val<1>{0});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == QDEPTH, "full queue stays full after oldest replacement");
    t.expect(read_mem(wb, old_head_addr) == 10, "oldest flushed to RAM when non-conflicting");

    // Case 4: full queue + conflicting read bank => drop oldest (no RAM commit) then overwrite slot.
    reset_case(wb);
    fill_full_with_enqueued_writes(wb, 30);
    t.expect(qcount_of(wb) == QDEPTH, "queue full again for conflict-drop case");
    old_head_addr = head_addr_of(wb);
    u64 rd_conf = find_addr_with_bank(bank_of(old_head_addr), true);
    wb.write(val<WBRW::A>{7}, val<8>{77}, val<1>{1}, val<WBRW::A>{rd_conf}, val<1>{1}, val<1>{0});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == QDEPTH, "full queue remains full after drop-oldest path");
    t.expect(read_mem(wb, old_head_addr) == 0, "oldest dropped on conflict (not committed to RAM)");

    // Case 5: full queue + conflicting read bank + noconflict=1 => force flush oldest.
    reset_case(wb);
    fill_full_with_enqueued_writes(wb, 40);
    t.expect(qcount_of(wb) == QDEPTH, "queue full again for forced-flush case");
    old_head_addr = head_addr_of(wb);
    rd_conf = find_addr_with_bank(bank_of(old_head_addr), true);
    wb.write(val<WBRW::A>{5}, val<8>{55}, val<1>{1}, val<WBRW::A>{rd_conf}, val<1>{1}, val<1>{1});
    hsu.next_cycle();
    t.expect(qcount_of(wb) == QDEPTH, "full queue stays full after forced flush replacement");
    t.expect(read_mem(wb, old_head_addr) == 40, "noconflict forces oldest flush under conflict");

    std::cout << "\n==== wb_rwram queue-policy test summary ====\n";
    std::cout << "passed: " << t.passed << "\n";
    std::cout << "failed: " << t.failed << "\n";

    return (t.failed == 0) ? 0 : 1;
}
