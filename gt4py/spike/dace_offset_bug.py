"""DaCe drops the +1 in p[c, 1, k + 1] read inside a loop nested in a map (2.0.0a10)."""
import dace
import numpy as np

N = dace.symbol("N")
L = dace.symbol("L")


@dace.program
def f(p: dace.float64[N, 2, L], out: dace.float64[N, L]):
    for c in dace.map[0:N]:
        for k in range(L - 1):
            out[c, k] = p[c, 1, k] + p[c, 1, k + 1]


p = np.arange(2 * 2 * 4, dtype=np.float64).reshape(2, 2, 4).copy()
out = np.zeros((2, 4))
f(p, out, N=2, L=4)
ref = np.zeros((2, 4))
ref[:, :3] = p[:, 1, :3] + p[:, 1, 1:]
print("dace:", out[:, :3].tolist())
print("expected:", ref[:, :3].tolist())
print("OK" if np.array_equal(out, ref) else "WRONG")
