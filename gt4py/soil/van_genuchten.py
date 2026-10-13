"""WaSiM's soil hydraulic functions as GT4Py (gt4py.cartesian) functions for column stencils.

WaSiM (src/defltuzr.cc, SOILtabletype) does not evaluate van Genuchten directly: at start-up it
tabulates theta(h) and k_rel(theta) at 30 fixed suctions per soil horizon (LookUpTableGeneration),
and h_calc, k_calc and get_theta interpolate linearly in these tables after a binary search.
These three calls take about 36% of a Rhone run (gt4py/README.md).

Storing the tables per cell and level would need 90 values per soil cell-level (several GB for
Switzerland). Instead these functions recompute the few table entries a lookup needs from the
horizon's van Genuchten parameters, with the same floating-point expressions as the table
generation, so the results are bit-identical to WaSiM (checked by test_van_genuchten.py against
the original C++ functions). The search is a fixed five-step bisection (30 entries), which keeps
GPU threads in lockstep.

Inputs per cell and level: thetas, thetar, alpha, n, m (= 1 - 1/n), tau. The fixed suctions are
passed as a 2-D field with a data dimension of 31 (index 0 unused), because gtscript cannot index
a constant tuple at run time; runtime data-dimension indexing needs the gt:* or dace:* backends.
"""
import numpy as np
import gt4py.cartesian.gtscript as gs

# h_tab of SOILtabletype::LookUpTableGeneration, suction in m, index 1..30 (index 0 never set there)
H_TAB = np.array([0.0, 0.0, 0.0005, 0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1,
                  0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6, 2.6, 3.45, 6.9, 20.0, 35.0, 50.0, 75.0,
                  100.0, 125.0, 150.0])
NTAB = 30


@gs.function
def theta_entry(htab, idx, ts, tr, alpha, n, m):
    """theta_h[idx] of the lookup table."""
    th = tr
    if idx == 0:
        th = ts
    elif idx < 30:
        th = (1.0 / (1.0 + (alpha * htab[0, 0][idx]) ** n)) ** m
        th = th * (ts - tr)
        th = th + tr
    return th


@gs.function
def k_entry(htab, idx, alpha, n, m, tau):
    """k_tab[idx] of the lookup table (relative conductivity)."""
    k = 0.0
    if idx < 30:
        temp1 = (1.0 / (1.0 + (alpha * htab[0, 0][idx]) ** n)) ** m
        temp2 = temp1 ** tau
        temp3 = (1.0 - (1.0 - temp1 ** (1.0 / m)) ** m) ** 2.0
        k = temp2 * temp3
    return k


@gs.function
def theta_probe(htab, theta, idx, step, ts, tr, alpha, n, m):
    nxt = idx + step
    if nxt <= 30:
        th = theta_entry(htab, nxt, ts, tr, alpha, n, m)
        if th > theta:
            idx = nxt
    return idx


@gs.function
def theta_search(htab, theta, ts, tr, alpha, n, m):
    """Last table index in 1..30 with theta_h > theta (the binary search of h_calc and k_calc)."""
    idx = 1
    idx = theta_probe(htab, theta, idx, 16, ts, tr, alpha, n, m)
    idx = theta_probe(htab, theta, idx, 8, ts, tr, alpha, n, m)
    idx = theta_probe(htab, theta, idx, 4, ts, tr, alpha, n, m)
    idx = theta_probe(htab, theta, idx, 2, ts, tr, alpha, n, m)
    idx = theta_probe(htab, theta, idx, 1, ts, tr, alpha, n, m)
    return idx


@gs.function
def h_probe(htab, h, idx, step):
    nxt = idx + step
    if nxt <= 30:
        if htab[0, 0][nxt] < h:
            idx = nxt
    return idx


@gs.function
def h_calc(htab, theta, ts, tr, alpha, n, m, thetarex, thetasex):
    """Suction in m for water content theta (SOILtabletype::h_calc)."""
    if thetasex > 0.0 and thetarex > 0.0:
        rel_theta = (theta - thetarex) / (thetasex - thetarex)
        theta = tr + (ts - tr) * rel_theta
    th1 = theta_entry(htab, 1, ts, tr, alpha, n, m)
    h = 0.0
    if theta <= tr:
        h = htab[0, 0][30]
    elif theta < th1:
        off = theta_search(htab, theta, ts, tr, alpha, n, m)
        t0 = theta_entry(htab, off, ts, tr, alpha, n, m)
        t1 = theta_entry(htab, off + 1, ts, tr, alpha, n, m)
        a = (theta - t1) / (t0 - t1)
        h = htab[0, 0][off] * a + htab[0, 0][off + 1] * (1.0 - a)
    return h


@gs.function
def k_calc(htab, theta, se, ts, tr, alpha, n, m, tau, thetarex, thetasex):
    """Relative hydraulic conductivity (0..1) for water content theta and liquid saturation se
    (SOILtabletype::k_calc)."""
    if thetasex > 0.0 and thetarex > 0.0:
        rel_theta = (theta - thetarex) / (thetasex - thetarex)
        theta = tr + (ts - tr) * rel_theta * se
    else:
        theta = tr + (theta - tr) * se
    if theta < tr:
        theta = tr
    th1 = theta_entry(htab, 1, ts, tr, alpha, n, m)
    k = 1.0
    if theta < th1:
        off = theta_search(htab, theta, ts, tr, alpha, n, m)
        t0 = theta_entry(htab, off, ts, tr, alpha, n, m)
        t1 = theta_entry(htab, off + 1, ts, tr, alpha, n, m)
        a = (theta - t1) / (t0 - t1)
        k = k_entry(htab, off, alpha, n, m, tau) * a + k_entry(htab, off + 1, alpha, n, m, tau) * (1.0 - a)
    return k


@gs.function
def get_theta(htab, h, ts, tr, alpha, n, m, thetarex, thetasex):
    """Water content for suction h in m (SOILtabletype::get_theta)."""
    res = tr
    if h < htab[0, 0][30]:
        off = 0
        off = h_probe(htab, h, off, 16)
        off = h_probe(htab, h, off, 8)
        off = h_probe(htab, h, off, 4)
        off = h_probe(htab, h, off, 2)
        off = h_probe(htab, h, off, 1)
        a = 0.0
        if off > 0:  # WaSiM reads the never-initialized h_tab[0] for h <= 0; with 0 there, a = 0
            if htab[0, 0][off + 1] - htab[0, 0][off] >= 0.001:
                a = (h - htab[0, 0][off]) / (htab[0, 0][off + 1] - htab[0, 0][off])
        res = theta_entry(htab, off + 1, ts, tr, alpha, n, m) * a + theta_entry(htab, off, ts, tr, alpha, n, m) * (1.0 - a)
        # as in WaSiM, the conversion to other thetar/thetas is skipped for h beyond the table
        if thetasex > 0.0 and thetarex > 0.0:
            rel_theta = (res - tr) / (ts - tr)
            res = thetarex + (thetasex - thetarex) * rel_theta
    return res
