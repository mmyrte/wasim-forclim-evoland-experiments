# CPU spike: hand-over, fusion and column loops

Status 2026-10-19. Code, specification and raw output: [`gt4py/spike/`](../gt4py/spike/) (`results-2026-10-19.txt`, `wasim_timing/`). Machine: this container, x86-64, 4 cores, no GPU; gt4py 1.2.3, DaCe 2.0.0a10, clang/gcc `-O3` without `-march=native`. All times are ms per hourly step, warm (without compilation).

The question: can a port reach WaSiM's speed, and which formulation should the final model use? Three parts:

1. the 2-D processes, called one by one as today versus with data kept in the backend and one call per process or one program per step;
2. WaSiM's per-cell loops in the Richards zone (sub-steps, iteration to convergence) in C++, DaCe, gt4py.cartesian and gt4py.next, including fixed loop counts that gt4py.next can run without host decisions;
3. WaSiM itself: one thread, OpenMP on 4 cores, and 4 runs side by side.

## Results

### WaSiM

| | Doubs (6,022 cells) | Rhone (21,584 cells) |
| --- | --- | --- |
| whole step, 1 thread | 76.7 | 184.4 |
| whole step, OpenMP 4 threads | 39.0 | 94.9 |
| of which the Richards zone, 1 / 4 threads | 59.9 / 24.5 | 115.9 / 49.4 |
| potential ET, interception, snow part 1, snow evaporation, 1 / 4 threads | 3.87 / 2.37 | 13.6 / 7.85 |

Four single-thread Doubs runs side by side take 22.1 s against 19.4 s for one run alone, so they have 3.5 times the throughput of one run; one run with 4 OpenMP threads takes 10.4 s. For ensembles, separate processes beat threads by a factor of 1.9.

The Richards zone uses 1 sub-step in 99.9 % of the Doubs cell-steps (maximum 2): its cost is the iteration to convergence.

### Richards zone: the proxy

The proxy (spec in [`gt4py/spike/README.md`](../gt4py/spike/README.md)) has the loop structure and cost of WaSiM's: single-threaded on Doubs 80 ms per step against WaSiM's 60, with 4 threads 22 against 24.5. Per cell-step it needs 5.9 iterations on average, the median is 4, the 99th percentile 33, the maximum 100 (`max_iter`).

| Formulation (adaptive: same algorithm) | Doubs, 1 thread | 4 threads | Rhone, 1 thread | 4 threads | Same result as C++ |
| --- | --- | --- | --- | --- | --- |
| C++, loops per cell (reference) | 80 | 22 | 284 | 73 | — |
| DaCe, `dace.map` over cells with the loops inside | 83 | 25 | 283 | 82 | bit-identical |
| DaCe, also the time loop inside the program | 81 | 25 | 295 | 82 | bit-identical |
| gt4py.cartesian, lockstep ticks from Python | 1,169 | 335 | 5,348 | 1,467 | bit-identical |
| gt4py.next, lockstep, one call per stage | 1,538 | 659 | 5,970 | 2,534 | bit-identical |

**Lockstep** means one stencil call per iteration ("tick") for all cells, with masks for the cells that are done, because neither gt4py flavour can put a column sweep inside a per-cell loop. A step then needs as many ticks as its slowest cell has iterations: 36.7 (Doubs) and 37.9 (Rhone) ticks per step, against 5.9 and 5.7 iterations per cell on average. So **6.2–6.6 times the work**, and each tick also costs 2.3–3 times an average C++ iteration over all cells (single thread: cartesian 32 ms, next 42 ms per tick, C++ 80 ms / 5.9 ≈ 14 ms). gt4py ends up 15–35 times slower than C++ and WaSiM.

gt4py.next called as one program per tick (all stages in one `@gtx.program`): 54 ms per tick on Doubs, 1.3 times the 42 ms of the stages called one by one; on gt4py's DaCe backend 48 ms (results differ in the last bit: that backend compiles with `-march=native`, so multiply-adds are fused).

**Fixed loop counts** (gt4py.next, no host decision during a step; N sub-steps, I iterations per sub-step, deviations of the 48-step totals from the adaptive result):

| N, I | ms per step, Doubs / Rhone (1 thread) | drainage | infiltration excess |
| --- | --- | --- | --- |
| 1, 4 / 6 / 8 | 172–354 / 695–1,278 | unstable (totals of 10^13–10^15 mm) | unstable |
| 1, 12 | 505 / 1,887 | +2.6 % / −3.0 % | +11.5 % / +10.7 % |
| 1, 16 | 687 / 2,618 | −10.1 % / −0.6 % | +8.2 % / +7.2 % |
| 2, 8 | 709 / 2,602 | +8.4 % / +1.0 % | +7.4 % / +7.5 % |

Cutting WaSiM's iteration short does not give a slightly less accurate result: below about 12 iterations the scheme diverges, and at 12–16 iterations the fluxes are several percent off, while each step still costs 6–9 times the C++ reference. Fixed counts are no way to keep this scheme and make it fast. A next-friendly column would need a different numerical scheme (for example an implicit step with a fixed number of Newton iterations, each a tridiagonal solve as two scans), which is a model change.

**On a GPU** the per-cell loops of the DaCe form would run in groups of 32 cells (a warp), each as long as its slowest cell: from the reference's iteration counts that is 20.4 iterations per cell-step against 5.9, about 3.5 times the work. That is about half of lockstep's waste of 6.2–6.6 times, and could shrink further by ordering cells by soil type. It is a GPU measurement still to make (wasim-gt#17).

### 2-D processes: hand-over and fusion

Potential ET (2 layers), interception (2 layers), snow evaporation and T-index snow for one daylight step, three ways, all bit-identical to each other. Members are copies of the grid stacked along I (an ensemble dimension):

| | Doubs, 1 thread | 4 threads | Rhone, 1 thread | 4 threads |
| --- | --- | --- | --- | --- |
| current process functions (host glue per call) | 53.5 | 46.0 | 57.8 | 39.0 |
| the same field operators called directly, fields kept in the backend | 21.3 | 11.7 | 29.1 | 11.1 |
| one `@gtx.program` for the whole chain | 27.9 | 18.4 | 36.1 | 22.2 |
| per member with 4 / 16 members (direct calls) | 23.4 / 27.8 | 11.5 / 10.7 | 32.3 | 11.9 |
| WaSiM, same modules (snow method 7, more expensive) | 3.87 | 2.37 | 13.6 | 7.85 |

- **Keeping the data in the backend** (no masks rebuilt, no `asnumpy`/`as_field` in a step) halves the time, or better with threads.
- **One gt4py.next program per chain is slower than direct calls** (by 1.2–2 times) and took 400 s to compile on first use. Composing programs does not pay on CPU with gt4py 1.2.3.
- **The time per member does not fall with more members.** This is not per-call overhead: it is arithmetic on the whole bounding rectangle of the grid, 44,421 points for Doubs (of which 6,022 model cells, 14 %) and 47,965 for Rhone (21,584, 45 %). Per grid point the direct calls take 0.48 µs (Doubs) and 0.61 µs (Rhone), WaSiM 0.64 µs per model cell. That is also why the earlier per-module times in `performance.md` hardly grew from Doubs to Rhone: the rectangles have the same size.
- **Computing only the model cells** would make the chain as fast as WaSiM: 0.48 µs × 6,022 cells ≈ 2.9 ms on Doubs against WaSiM's 3.9 ms.

## Conclusions

1. **Equal performance with WaSiM on CPU is in reach, but not with gt4py for the columns.** DaCe with WaSiM's per-cell loops ran at C++ speed and bit-identical; both gt4py flavours have to run those loops in lockstep from the host and were 15–35 times slower. Fixed loop counts are unstable or several percent off for WaSiM's scheme.
2. **For the 2-D processes the lever is the layout, not fusion:** compute model cells only (a compressed cell dimension, with neighbour tables where a process needs neighbours) and keep fields in the backend between processes. Fusing calls into one gt4py program made things slower.
3. **Ensembles:** WaSiM already gets 3.5 times the throughput from 4 processes on 4 cores. A port has to beat that, not WaSiM's single-thread time.
4. **Way forward** (wasim-gt `docs/architecture.md` §9): keep replicating module by module with gt4py, bit for bit; then build the fast model directly in DaCe, with a compressed cell layout and the per-cell loops as they are, and measure on a GPU.

## Caveats

- The column result comes from a proxy, not WaSiM's Richards zone; the proxy's cost matches WaSiM's on Doubs (1 thread 80 against 60 ms, 4 threads 22 against 24.5 ms) and its iteration counts are what drives the comparison.
- One machine, 4 cores, no GPU. The fused-chain inputs are one real daylight step with made-up precipitation and snow (the dump has none).
- DaCe 2.0.0a10 miscompiled one access pattern (`gt4py/spike/dace_offset_bug.py`: an index offset dropped inside a map); its default flags include `-ffast-math`. Both need regression tests before the model depends on DaCe.
