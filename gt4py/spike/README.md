# CPU spike: fusion, column loops and coordination cost

Three questions, answered on CPU (this container: 4 cores, no GPU) before more modules are ported:

1. What does it cost to call the 2-D processes one by one, and what does one fused call per step
   cost instead? (`fused_chain.py`)
2. How can WaSiM's per-cell loops in the Richards zone (adaptive sub-steps, iteration to
   convergence) be written in gt4py.cartesian, gt4py.next and DaCe, how fast is each compared to
   C++, and how close do fixed, next-compatible loop counts get to the adaptive result?
   (`richards_*`)
3. How do these compare with WaSiM itself, single-threaded, with OpenMP and as several
   processes at once? (`WASIM_TIMING`, see `../../validation/performance.md`)

Results and conclusions: [`../../validation/spike-cpu.md`](../../validation/spike-cpu.md).

## The Richards proxy

The Richards zone is not ported yet (wasim-gt#7), so the column question is answered with a proxy
that has the structure and cost of WaSiM's loop (`unsatzon.cc`, `run`, method 3) without its many
couplings. Per cell and time step `DT`:

1. **Sub-steps** (Courant criterion, as WaSiM): `dtmin = min(DT, min over levels 1..L-1 of
   dz / (k_rel ks))`, with `k_rel = 1` when it rains and `k_rel(θ)` otherwise; at least
   `min_step`; `n = ceil(DT / dtmin)`, `dt = DT / n`; rain and ET are split evenly.
2. Per sub-step: ET taken from the top `root_levels` levels (not below θr), then `old = avg = θ`,
   `h1 = h(old)`, and iterate (`z` = iteration count):
   - heads `hh = zc − (h1 + h(avg)) / 2`, conductivity `K = ks k_rel(avg)` (van Genuchten–Mualem);
   - Darcy fluxes between levels `q = sqrt(K_k K_k+1) (hh_k − hh_k+1) / ((dz_k + dz_k+1) / 2)`,
     rain at the top, free drainage `K_L` at the bottom (down positive);
   - restrictions: top-down, no level below θr (limits outflow); bottom-up, no level above θs
     (limits inflow; at the top this is infiltration excess);
   - new water contents `θ = old + (q_top − q_bottom) dt / dz` (clipped to 0..1) and WaSiM's
     moving average of `new_thetas` (`avg = θ` in the first 3 iterations, then weighted by `z − 3`);
   - until the largest change is ≤ `max_err` (5 × after 20 iterations) or `z = max_iter`.
3. Outputs: θ, drainage, infiltration excess, ET, sub-steps and iterations per cell and step.

Not in the proxy: groundwater coupling, interflow, drainage, macropores, frost, heat transport,
snow layers, WaSiM's lookup tables for h(θ) and k(θ). Case (`proxy_case.py`): the 7 soils of the
Doubs soil table (20 levels each), 6,022 (Doubs) or 21,584 (Rhone) columns, 48 hourly steps with
two rain events per day and daytime ET; `max_err = 1e-5`, `max_iter = 100` as in the golden cases.
The proxy is calibrated by its cost: single-threaded it takes 80 ms per Doubs step, WaSiM's
Richards zone 60 ms; with 4 threads 22 and 24.5 ms.

## Formulations

| File | Loops over sub-steps and iterations | Columns | Same result as C++ |
| --- | --- | --- | --- |
| `richards_ref.cpp` | per cell, in C++ (OpenMP over cells) | per cell, in C++ | reference |
| `richards_dace.py` | per cell, Python loops in a `dace.map` (compiled) | per cell | bit-identical |
| `richards_cartesian.py` | lockstep "ticks" from Python, masks per cell | `FORWARD`/`BACKWARD` computations | bit-identical (adaptive) |
| `richards_next.py` | lockstep ticks from Python, masks per cell | scan operators | bit-identical (adaptive) |

Lockstep: a tick is one iteration for every cell that is not done; a cell that converges
finishes its sub-step and starts the next in the next tick. `adaptive` asks the host after every
tick whether all cells are done; `fixed` uses N sub-steps for every cell and at most I iterations
per sub-step and runs exactly N·I ticks, with no host check.

## Running

```sh
uv run --project ../../../wasim-gt python proxy_case.py WORK/case_doubs          # or use run_all.sh
c++ -O3 -fopenmp -o richards_ref richards_ref.cpp && OMP_NUM_THREADS=1 ./richards_ref WORK/case_doubs OUT 3
python richards_dace.py WORK/case_doubs OUT [--loop]
python richards_cartesian.py WORK/case_doubs OUT [--mode fixed --substeps 1 --iters 8]
python richards_next.py WORK/case_doubs OUT [--mode fixed --substeps 1 --iters 8]
python compare.py WORK/case_doubs REF_OUT OUT
python fused_chain.py DUMP CONTROL_INTERN 2021-06-13T01 --members 1 4 16
./run_all.sh WORK ../../../wasim-gt        # everything, sequentially
```

DaCe needs `DACE_compiler_cpu_args="-std=c++17 -fPIC -O3 -fopenmp"`: its default adds
`-march=native -ffast-math`, which changes results.

## Found on the way

- **DaCe 2.0.0a10 drops an index offset** inside a `dace.map`: `p[c, 1, k] + p[c, 1, k + 1]` in a
  sequential loop reads `p[c, 1, k]` twice. Workaround in `richards_dace.py`: copy the column into a
  local array first. Fixed in DaCe main, not yet released; reproducer and details in
  [`../../dace/known-issues/`](../../dace/known-issues/README.md).
- **gt4py.next on gtfn:** a scan cannot take another scan's output within one field operator
  ("Scans are only supported at the top level of a stencil closure"), so each scan is its own
  stage of a program; one large field operator for a whole tick did not finish compiling within
  15 minutes (type inference on the inlined tree), and a program of many stages takes minutes to
  compile (400 s for `fused_chain.py`'s) and runs slower than calling its stages one by one.
- **gt4py.next's build cache** is per process by default (`GT4PY_BUILD_CACHE_LIFETIME=session`):
  every new process recompiles. `run_all.sh` sets it to `persistent`.
- **WaSiM defect** (now in wasim-eth `TODO.md`): `if (gwin[k] += qu_sl[k] >= 0)` in
  `flux_restrictions(_opt)`.
