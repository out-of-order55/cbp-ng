# CBP-NG Championship Tutorial

Welcome to the tutorial for the Next-Generation Championship in Branch
Prediction! This document aims to provide a basic introduction to the
prediction interface and the HARCOM library. It does not cover all of HARCOM's
features - for that you are referred to [the HARCOM
manual](https://github.com/AmpereComputing/cbp-ng/raw/refs/heads/main/docs/harcom.pdf).

In this tutorial, we will introduce the interface and basic HARCOM topics by
considering a number of working CBP-NG predictors. The source code for each
predictor can be found in the `./predictors/tutorial` directory, and they may
be compiled and run like any other predictor, for example:

```console
./compile cbp -DPREDICTOR="tutorial_00"
./cbp ./gcc_test_trace.gz gcc_test_trace 0 40000
```

## Predictor 00 ([tutorial_00.hpp](/predictors/tutorial/tutorial_00.hpp))

The first predictor we'll consider is the simplest - it predicts not-taken for
every instruction. Though the predictor itself isn't very interesting, we can
use it to study the interface between the simulator and the predictor. Here is
a flowchart of the calls the simulator makes into predictors:

![predictor interface flowchart](/docs/interface-flowchart.svg)

Note that `predict1()` and `predict2()` are each called once for the first
instruction in each prediction block. Because our predictor never predicts more
than one instruction per block, `reuse_predict1()` and `reuse_predict2()` are
never called (more on this later).

The arguments and return types of all methods involved in the predictor
interface are `val` - a basic HARCOM type. A `val` represents a value with a
particular type and bit width. They can be acted upon by other HARCOM C++
operations (via operator overloading or named methods). In addition to its
value, each `val` carries a notion of its 'time' - this represents when the
value was first available to be used. When an operation is performed on a
HARCOM `val`, the result is assigned the time of its input(s) plus the time it
took to perform the operation. In addition to latency (time), HARCOM also
records the energy and area used by each operation.

For example, the time of the `val<64> inst_pc` argument to `predict1()` will be
when that prediction could begin (i.e. when the target of the previous
prediction was known). The latency of the level-1 predictor is the difference
between the return value's time and the time of its input, `inst_pc`.

Our first predictor is a little "odd" - it returns `hard<0>{}` for all
predictions instead of a `val`. The `hard` HARCOM type acts as a hardware
parameter/constant, and the template parameter (0 here) is its value. They
always have time 0 because they don't depend on other values or logic. You use
`hard` to introduce constant values into your predictor's logic.

Back to our predictor's return value: returning 0 signifies "not-taken" (1 is
naturally "taken").

## Predictor 01 ([tutorial_01.hpp](/predictors/tutorial/tutorial_01.hpp))

The next predictor stores information about past branches to try to improve
prediction. In particular, it stores the PC (program counter, a.k.a.
instruction pointer) of the last taken taken. If it finds a branch with a PC
matching our saved PC, it predicts it taken.

This predictor introduces a new type, `reg` (register). `reg`s are unlike
`val`s in that they may be updated (although they may only be written once per
cycle). They can be used to persist information between cycles/calls. You can
see that in both `predict1()` and `predict2()`, the value of `reg<64>
last_taken` is compared to the current instruction's PC, `val<64> inst_pc` with
`==`. The result of this comparison is a `val<1>` set to 1 if equal and 0
otherwise. It can directly return this result, because the equality test will
return 1 when the current PC matches the last taken branch.

This predictor must now also include some update logic. It stores the PC of
taken branches, and overwrite `last_taken` whenever it sees a new taken branch.
It accomplishes the first part by saving the branch PC and direction to two
other registers in `update_condbr()`, which is called after prediction for each
conditional branch.

Then in `update_cycle()`, called once for each prediction block, it
conditionally updates `last_taken`. If it mispredicted the just-predicted
branch it updates the PC (otherwise there is no need). And if it mispredicted
and the current branch was not taken, it resets `last_taken` to be 0 (this
means the current branch was not taken but its PC is stored in `last_taken`).
It accomplishes both of these conditional operations using the `select()`
method. If the first argument to `select()` is 1, it returns the second
argument; otherwise, it returns the third. This is similar to the ternary
operator in C++.

## Predictor 02 ([tutorial_02.hpp](/predictors/tutorial/tutorial_02.hpp))

The third predictor is a simple bimodal predictor. Instead of storing one piece
of global information like the previous example, it uses the current PC to
index into a RAM array (a HARCOM `ram` type) of 2-bit saturating counters.

The prediction logic (`predict1()`) is straightforward - it hashes the PC to
form an index, uses that index to access the RAM array, and returns the top bit
of the saturating counter as the direction prediction. The calculation of the
index introduces two new HARCOM features: `make_array(val<6>{})` chunks the
64-bit `inst_pc` into a HARCOM `arr` made up of 6-bit `val`s. `fold_xor()`,
then uses a series of exclusive-or operations to combine all of the 6-bit
`val`s in the `arr` into a single 6-bit index. `counters.read(index)` returns
the current value at the specified index.

To update, it increases the read counter value of the branch was executed
taken, and decrease it if not. The catch is that HARCOM `ram`s can only perform
one read OR one write per cycle, and it was already read at prediction-time
(you will get a runtime error if you violate this rule). One way to work around
this "hardware" limitation is to inform the simulator that the predictor needs
an additional cycle by calling `need_extra_cycle` with a 1. Rather than always
incur a 1-cycle update penalty, this predictor conditionally calculates whether
an update is needed by checking if the counter was already saturated, and uses
`execute_if` to conditionally update the RAM.

## Predictor 03 ([tutorial_03.hpp](/predictors/tutorial/tutorial_03.hpp))

The main novelty of this predictor is that unlike previous tutorial predictors
it predicts more than one instruction per prediction block. Instead of storing
a single counter per RAM entry like predictor_02.hpp, it stores an array of
counters per entry (one per possible instruction in the block). It calls
`reuse_prediction(1)` to indicate the predictor is ready to predict an
additional instruction in this block (traditionally to re-use the same RAM
lookup(s)). `predict1()` and `predict2()` are called to predict the first
instruction in the block; `reuse_predict1()` and `reuse_predict2()` are called
for the remaining instructions (until a branch is mispredicted, executed taken,
or the predictor calls `reuse_prediction(0)`).

You can also see that handling updates for more than one branch at a time has
added some complexity. As conditional branches are observed in
`update_condbr()`, their PCs and directions are recorded into arrays of
registers. Note that those arrays are indexed using a C++ integer
(`num_branches`)! Though it is expected that most of the types in your
predictors will be HARCOM types, this is an exception - it is allowed to count
calls and use those counts in this way. In fact, it would be difficult to
structure this update logic if it were disallowed.

The logic in `update_cycle()` reads the accumulated branch information and uses
it to update the counters for all conditional branches executed in the block,
before saving them back to the RAM (if needed).

Note: Because this predictor is the first in the tutorial to be a template, you
must either add "<>" to the end of the class name to accept the default
template parameters or something like "<3>" to build it with 2^3 instructions
per block:

```console
./compile cbp -DPREDICTOR="tutorial_03<>"
```
