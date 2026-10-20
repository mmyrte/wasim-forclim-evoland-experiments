# DaCe: known issues

Defects of DaCe met while building wasim-gt's DaCe form (wasim-gt `docs/architecture.md` §2, §9). Each
has a script that prints OK or WRONG per case and exits with the number of wrong cases, so it can run
as a regression test in wasim-gt's CI against the pinned DaCe version.

## Loop indices inside a map get the wrong offset (`loop_in_map_offset.py`)

**Status:** wrong in DaCe **2.0.0a10 from PyPI** (the version gt4py 1.2.3 resolves to); **fixed in DaCe
main** (correct at `be19bf4a98f5161ed258b25ad01fcc469083af08`, 2026-10-20). Not reported upstream, as
main already has it right. Pin a DaCe version that contains the fix, and keep the script as a test.

**What happens.** Inside `for c in dace.map[0:N]`, a sequential `for k in range(...)` that indexes an
argument array with an offset of the loop variable loses the offset:

```python
@dace.program
def f(p: dace.float64[N, L], out: dace.float64[N, L]):
    for c in dace.map[0:N]:
        for k in range(L - 1):
            out[c, k] = p[c, k + 1]      # reads p[c, k]
```

| Case | 2.0.0a10 (PyPI) | main be19bf4 |
| --- | --- | --- |
| map + loop, `p[c, k + 1]` | wrong: reads `p[c, k]` | OK |
| map + loop over `range(1, L)`, `p[c, k - 1]` | wrong: writes `out[c, k - 1] = p[c, k - 1]` | OK |
| map + loop, `p[c, 1, k] + p[c, 1, k + 1]` | wrong: `2 · p[c, 1, k]` | OK |
| two loops, no map | OK | OK |
| 2-D map over `c, k` | OK | OK |
| map + loop, column copied into a local array first | OK | OK |

Same results on Python 3.13.16 and 3.14.6 (numpy 2.5.4, gcc 13.3, x86-64), and with
`to_sdfg(simplify=False)`. The generated C++ passes the column's base pointer and indexes it with `k`.

**Where it bit.** The Richards proxy of the CPU spike (`gt4py/spike/richards_dace.py`) read
`p[c, DZ, k + 1]` in the flux between two levels, so its results differed from the C++ reference only at
horizon boundaries, where the layer thickness changes. The workaround there is the last case of the
table: copy the column into a local array, then index the local array.

**Run:**

```sh
python dace/known-issues/loop_in_map_offset.py    # exit code = number of wrong cases
```
