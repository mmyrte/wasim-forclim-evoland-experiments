#!/usr/bin/env python3
"""DaCe: array indices that depend on a sequential loop nested in a map get the wrong offset.

Inside `for c in dace.map[...]`, a sequential `for k in range(...)` that reads `p[c, k + 1]` reads
`p[c, k]`; with `range(1, L)` and `p[c, k - 1]` the write index shifts as well. The same code without
the map, or as a 2-D map, is correct. Wrong in DaCe 2.0.0a10 from PyPI (Python 3.13.16 and 3.14.6,
numpy 2.5.4, gcc 13.3, x86-64), also with `to_sdfg(simplify=False)`; correct in DaCe main at
be19bf4a98f5161ed258b25ad01fcc469083af08 (2026-10-20) on both Python versions.

Each case compares with numpy and prints OK or WRONG; the exit code is the number of wrong cases, so
the script can serve as a regression test (see README.md).

usage: python loop_in_map_offset.py
"""
import sys

import dace
import numpy as np

N = dace.symbol("N")
L = dace.symbol("L")


@dace.program
def map_loop_plus(p: dace.float64[N, L], out: dace.float64[N, L]):
    for c in dace.map[0:N]:
        for k in range(L - 1):
            out[c, k] = p[c, k + 1]


@dace.program
def map_loop_minus(p: dace.float64[N, L], out: dace.float64[N, L]):
    for c in dace.map[0:N]:
        for k in range(1, L):
            out[c, k] = p[c, k - 1]


@dace.program
def map_loop_sum_3d(p: dace.float64[N, 2, L], out: dace.float64[N, L]):
    for c in dace.map[0:N]:
        for k in range(L - 1):
            out[c, k] = p[c, 1, k] + p[c, 1, k + 1]


@dace.program
def loop_loop_plus(p: dace.float64[N, L], out: dace.float64[N, L]):
    for c in range(N):
        for k in range(L - 1):
            out[c, k] = p[c, k + 1]


@dace.program
def map2d_plus(p: dace.float64[N, L], out: dace.float64[N, L]):
    for c, k in dace.map[0:N, 0:L - 1]:
        out[c, k] = p[c, k + 1]


@dace.program
def map_loop_local_copy(p: dace.float64[N, L], out: dace.float64[N, L]):
    """The workaround used in gt4py/spike/richards_dace.py: copy the column into a local array first."""
    for c in dace.map[0:N]:
        col = np.empty(L, dtype=np.float64)
        for k in range(L):
            col[k] = p[c, k]
        for k in range(L - 1):
            out[c, k] = col[k + 1]


def main():
    p2 = np.arange(8, dtype=np.float64).reshape(2, 4)
    p3 = np.arange(16, dtype=np.float64).reshape(2, 2, 4)
    plus = np.zeros((2, 4))
    plus[:, :3] = p2[:, 1:]
    minus = np.zeros((2, 4))
    minus[:, 1:] = p2[:, :3]
    sum3 = np.zeros((2, 4))
    sum3[:, :3] = p3[:, 1, :3] + p3[:, 1, 1:]
    cases = [
        ("map + loop, p[c, k + 1]", map_loop_plus, p2, plus),
        ("map + loop, p[c, k - 1]", map_loop_minus, p2, minus),
        ("map + loop, p[c, 1, k] + p[c, 1, k + 1]", map_loop_sum_3d, p3, sum3),
        ("loop + loop (no map), p[c, k + 1]", loop_loop_plus, p2, plus),
        ("2-D map, p[c, k + 1]", map2d_plus, p2, plus),
        ("map + loop, local copy first", map_loop_local_copy, p2, plus),
    ]
    wrong = 0
    for name, prog, p, expected in cases:
        out = np.zeros((2, 4))
        prog(p.copy(), out, N=2, L=4)
        ok = np.array_equal(out, expected)
        wrong += not ok
        print(f"{'OK   ' if ok else 'WRONG'} {name}: got {out[0].tolist()}, expected {expected[0].tolist()}")
    print(f"dace {dace.__version__}: {wrong} wrong case(s)")
    return wrong


if __name__ == "__main__":
    sys.exit(main())
