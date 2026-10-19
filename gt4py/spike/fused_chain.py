#!/usr/bin/env python3
"""Fused 2-D chain: potential ET, interception, snow evaporation, T-index snow as one program.

The wasim-gt processes are called one by one today (potential_et, interception,
snow_evaporation.conventional, snow_tindex), each through Python glue that builds masks on the
host and calls one field operator after the other. Here the same field operators of wasim-gt are
called from one gt4py.next program, so a step is one call with all fields staying in the
backend's memory, and the per-call overhead is paid once. Both paths get the same inputs and
must give the same outputs (bit-identical).

Inputs: one daylight step of a WaSiM grid dump (as benchmark_et.py); precipitation, interception
capacities and storages, and snow storages are made up (the dump has none), the same for both
paths. Only steps without a dark part are handled (the driver's night calls are left out).

Ensemble: --members E stacks E copies of the grids along I, as E parameter sets would be.

usage: fused_chain.py DUMP_DIR CONTROL_INTERN START [--step T] [--repeat N] [--members E ...]
                      [--backend gtfn_cpu|dace_cpu]
"""
import argparse
import datetime
import os
import sys
import time

import numpy as np
import gt4py.next as gtx
from gt4py.next import float64

from wasim_gt.fields import F, I, J, M, with_backend
from wasim_gt.processes import interception as icp
from wasim_gt.processes import potential_et as pet
from wasim_gt.processes import snow_evaporation as sev
from wasim_gt.processes.interception import glacier_mix, layer_step
from wasim_gt.processes.potential_et import add, combine_layer, get_etp, lai_sum, pick, shortwave
from wasim_gt.processes.snow_evaporation import take
from wasim_gt.processes.snow_tindex import snow_step
from wasim_gt.processes.sun import Site, daylight

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from potential_et import multilayer_table, numbers, section  # noqa: E402
from wasim_dump import Dump  # noqa: E402

SNOW = dict(t0r=1.0, transition=1.0, t0=0.0, c0=2.5, cwh=0.1)


@gtx.program
def hydro_step(
    rg: F, ssd: F, ssd_mean: F, exkorr: F, albedo: F, lake: F, pond: F, use_lake: bool,
    cr0: float64, cr1: float64, cr2: float64, cr3: float64, rgmax: float64, bright: bool, albedo_ow: float64,
    tm: F, dhm: F, rlf: F, u: F, rse: F, k_lambert: F, lai_scale: F, pond_cell: M,
    lai1: F, vbg1: F, rsc1: F, rsi1: F, z01: F, veg1: M, bottom1: M,
    lai2: F, vbg2: F, rsc2: F, rsi2: F, z02: F, veg2: M, bottom2: M,
    t_first: bool, f_first: bool, night: bool, anteil: float64, sigma: float64, vv_factor: float64, z0_ow: float64,
    prec: F, cap1: F, cap2: F, si1: F, si2: F, ice: F,
    ssnow: F, sliq: F, no_snow: M, snow_veg1: M, snow_veg2: M, snow: M,
    t0r: float64, transition: float64, t0: float64, c0: float64, cwh: float64, steps_per_day: float64,
    zero: F,
    # temporaries
    rk: F, relrg: F, ssd_used: F, e1: F, ev1: F, ei1: F, ra1: F, rr1: F, ra_l1: F, rr_l1: F, etp0a: F, eip0a: F,
    lai_kum2: F, e2: F, ev2: F, ei2: F, ra2: F, rr2: F, etp0b: F, eip0b: F,
    iei1: F, iq1: F, isi1: F, ietp1: F, ietp0: F, iei0: F, isi0: F, iq0: F,
    iei2: F, iq2: F, isi2: F, ietp2: F, ietp0b: F, iei0b: F, isi0b: F, iq0b: F,
    s1: F, d1: F, v1: F, s2: F, d2: F, v2: F,
    # outputs
    rg_o: F, rgex: F, cloudiness: F, etp1: F, eip1: F, evap1: F, etp2: F, eip2: F, evap: F, etp_total: F,
    throughfall: F, ssnow_o: F, sliq_o: F, outflow: F, evap_o: F, evrs: F,
):
    shortwave(rg, ssd, ssd_mean, exkorr, albedo, lake, pond, use_lake, cr0, cr1, cr2, cr3, rgmax, bright,
                  albedo_ow, out=(rg_o, rgex, rk, relrg, ssd_used, cloudiness))
    get_etp(tm, dhm, rlf, u, rse, ssd_used, rk, lai1, vbg1, rsc1, rsi1, z01, k_lambert, lai_scale, pond_cell,
                bottom1, zero, zero, zero, t_first, night, anteil, sigma, vv_factor, cr0, z0_ow,
                out=(e1, ev1, ei1, ra1, rr1))
    pick(veg1, ra1, zero, out=ra_l1)
    pick(veg1, rr1, zero, out=rr_l1)
    combine_layer(veg1, zero, zero, zero, e1, ei1, ev1, zero, zero, zero, out=(etp1, eip1, evap1, etp0a, eip0a))
    lai_sum(veg2, zero, lai1, out=lai_kum2)
    get_etp(tm, dhm, rlf, u, rse, ssd_used, rk, lai2, vbg2, rsc2, rsi2, z02, k_lambert, lai_scale, pond_cell,
                bottom2, lai_kum2, ra_l1, rr_l1, f_first, night, anteil, sigma, vv_factor, cr0, z0_ow,
                out=(e2, ev2, ei2, ra2, rr2))
    combine_layer(veg2, zero, zero, zero, e2, ei2, ev2, evap1, etp0a, eip0a, out=(etp2, eip2, evap, etp0b, eip0b))
    add(etp0b, evap, out=etp_total)
    layer_step(prec, lai1, vbg1, cap1, etp1, si1, veg1, etp_total, zero, zero, zero,
                   out=(iei1, iq1, isi1, ietp1, ietp0, iei0, isi0, iq0))
    layer_step(iq1, lai2, vbg2, cap2, etp2, si2, veg2, ietp0, iei0, isi0, iq0,
                   out=(iei2, iq2, isi2, ietp2, ietp0b, iei0b, isi0b, iq0b))
    glacier_mix(iq0b, prec, ice, out=throughfall)
    take(ssnow, ietp1, snow_veg1, zero, out=(s1, d1, v1))
    take(s1, ietp2, snow_veg2, v1, out=(s2, d2, v2))
    take(s2, evap, snow, v2, out=(ssnow_o, evap_o, evrs))
    snow_step(tm, throughfall, evrs, ssnow, sliq, no_snow, t0r, transition, t0, c0, cwh, steps_per_day,
              out=(ssnow_o, sliq_o, outflow))


def hydro_step_lean(b, a):
    """The program hydro_step as direct field-operator calls from Python (fields preallocated, no host glue)."""
    w = lambda op: with_backend(op, b)  # noqa: E731  (cached: binding on every call is slow)
    o = {}
    w(shortwave)(a["rg"], a["ssd"], a["ssd_mean"], a["exkorr"], a["albedo"], a["lake"], a["pond"], a["use_lake"],
                 a["cr0"], a["cr1"], a["cr2"], a["cr3"], a["rgmax"], a["bright"], a["albedo_ow"],
                 out=(a["rg_o"], a["rgex"], a["rk"], a["relrg"], a["ssd_used"], a["cloudiness"]), offset_provider=o)
    w(get_etp)(a["tm"], a["dhm"], a["rlf"], a["u"], a["rse"], a["ssd_used"], a["rk"], a["lai1"], a["vbg1"], a["rsc1"],
               a["rsi1"], a["z01"], a["k_lambert"], a["lai_scale"], a["pond_cell"], a["bottom1"], a["zero"], a["zero"],
               a["zero"], True, False, a["anteil"], a["sigma"], a["vv_factor"], a["cr0"], a["z0_ow"],
               out=(a["e1"], a["ev1"], a["ei1"], a["ra1"], a["rr1"]), offset_provider=o)
    w(pick)(a["veg1"], a["ra1"], a["zero"], out=a["ra_l1"], offset_provider=o)
    w(pick)(a["veg1"], a["rr1"], a["zero"], out=a["rr_l1"], offset_provider=o)
    w(combine_layer)(a["veg1"], a["zero"], a["zero"], a["zero"], a["e1"], a["ei1"], a["ev1"], a["zero"], a["zero"],
                     a["zero"], out=(a["etp1"], a["eip1"], a["evap1"], a["etp0a"], a["eip0a"]), offset_provider=o)
    w(lai_sum)(a["veg2"], a["zero"], a["lai1"], out=a["lai_kum2"], offset_provider=o)
    w(get_etp)(a["tm"], a["dhm"], a["rlf"], a["u"], a["rse"], a["ssd_used"], a["rk"], a["lai2"], a["vbg2"], a["rsc2"],
               a["rsi2"], a["z02"], a["k_lambert"], a["lai_scale"], a["pond_cell"], a["bottom2"], a["lai_kum2"],
               a["ra_l1"], a["rr_l1"], False, False, a["anteil"], a["sigma"], a["vv_factor"], a["cr0"], a["z0_ow"],
               out=(a["e2"], a["ev2"], a["ei2"], a["ra2"], a["rr2"]), offset_provider=o)
    w(combine_layer)(a["veg2"], a["zero"], a["zero"], a["zero"], a["e2"], a["ei2"], a["ev2"], a["evap1"], a["etp0a"],
                     a["eip0a"], out=(a["etp2"], a["eip2"], a["evap"], a["etp0b"], a["eip0b"]), offset_provider=o)
    w(add)(a["etp0b"], a["evap"], out=a["etp_total"], offset_provider=o)
    w(layer_step)(a["prec"], a["lai1"], a["vbg1"], a["cap1"], a["etp1"], a["si1"], a["veg1"], a["etp_total"], a["zero"],
                  a["zero"], a["zero"], out=(a["iei1"], a["iq1"], a["isi1"], a["ietp1"], a["ietp0"], a["iei0"], a["isi0"],
                                              a["iq0"]), offset_provider=o)
    w(layer_step)(a["iq1"], a["lai2"], a["vbg2"], a["cap2"], a["etp2"], a["si2"], a["veg2"], a["ietp0"], a["iei0"],
                  a["isi0"], a["iq0"], out=(a["iei2"], a["iq2"], a["isi2"], a["ietp2"], a["ietp0b"], a["iei0b"],
                                            a["isi0b"], a["iq0b"]), offset_provider=o)
    w(glacier_mix)(a["iq0b"], a["prec"], a["ice"], out=a["throughfall"], offset_provider=o)
    w(take)(a["ssnow"], a["ietp1"], a["snow_veg1"], a["zero"], out=(a["s1"], a["d1"], a["v1"]), offset_provider=o)
    w(take)(a["s1"], a["ietp2"], a["snow_veg2"], a["v1"], out=(a["s2"], a["d2"], a["v2"]), offset_provider=o)
    w(take)(a["s2"], a["evap"], a["snow"], a["v2"], out=(a["ssnow_o"], a["evap_o"], a["evrs"]), offset_provider=o)
    w(snow_step)(a["tm"], a["throughfall"], a["evrs"], a["ssnow"], a["sliq"], a["no_snow"], a["t0r"], a["transition"],
                 a["t0"], a["c0"], a["cwh"], a["steps_per_day"], out=(a["ssnow_o"], a["sliq_o"], a["outflow"]),
                 offset_provider=o)


def backend_of(name):
    if name == "gtfn_cpu":
        return gtx.gtfn_cpu
    from gt4py.next.program_processors.runners.dace import run_dace_cpu
    return run_dace_cpu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_intern")
    ap.add_argument("start")
    ap.add_argument("--step", type=int, default=11)
    ap.add_argument("--repeat", type=int, default=20)
    ap.add_argument("--members", type=int, nargs="+", default=[1])
    ap.add_argument("--backend", default="gtfn_cpu")
    ap.add_argument("--no-program", action="store_true", help="skip the fused program (its compilation is slow)")
    args = ap.parse_args()

    text = open(args.control_intern, encoding="latin-1").read()
    coords, et, lake = (section(text, s) for s in ("coordinates", "evapotranspiration", "lake_model"))
    step = float(et[1])
    site = Site(float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3]),
                float(section(text, "output_interval")[3].split()[0]) / 60.0)
    amp = numbers(et[30])
    lakes_on = int(lake[0]) != 0
    params = pet.Parameters(cr=tuple(numbers(et[25])), dt_koeff=numbers(et[28])[0], dtmax=tuple(numbers(et[29])[:12]),
                            ampl_day=amp[0], ampl_annual=amp[1], albedo_open_water=float(lake[2]) if lakes_on else 0.0,
                            z0_open_water=float(lake[3]) if lakes_on else 0.0, lakes=lakes_on)
    d = Dump(args.dump_dir)
    p1, t = "after_expokorr", args.step
    date, hour = args.start.split("T")
    now = datetime.datetime.strptime(date, "%Y-%m-%d") + datetime.timedelta(hours=int(hour), minutes=step * t)
    dl = daylight(site, now.day, now.month, int(now.hour + site.time_shift), step)
    if dl.dark >= 0.01 or dl.no_radiation:
        raise SystemExit(f"step {t} has a dark part ({dl.dark}); pick a full daylight step")
    b = backend_of(args.backend)
    rng = np.random.default_rng(0)

    for members in args.members:
        base = {g: np.array(d(p1, g)[t]) for g in (
            "[zone_grid]", "[elevation_model]", "landuse", "landuse1", "landuse2", "lake_codes", "ponding_storage",
            "albedo", "SurfaceEvaporationResistance", "corr_temperature", "global_radiation", "sunshine_duration",
            "radiation_correction", "air_humidity", "wind_speed", "SSDTemp")
            + tuple(f"{g}{n}" for g in ("vegetation_coverage_degree", "leaf_area_index", "SurfaceCanopyResistance",
                                        "SurfaceIntercepResistance", "RoughnessLength") for n in (1, 2))}
        tile = lambda a: np.concatenate([a] * members, axis=0)  # noqa: E731
        g = {k: tile(v) for k, v in base.items()}
        active = g["[zone_grid]"] > -9999
        nrows, ncols = active.shape
        dom = gtx.domain({I: nrows, J: ncols})
        fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))  # noqa: E731
        raw = lambda a: gtx.as_field(dom, np.array(a))  # noqa: E731
        msk = lambda a: gtx.as_field(dom, np.array(a, dtype=bool))  # noqa: E731
        ml = multilayer_table(text)
        kl, ls = np.full(active.shape, 0.3), np.full(active.shape, 10.0)
        for code, (k, s) in ml.items():
            kl[g["landuse"] == code], ls[g["landuse"] == code] = k, s
        prec = np.where(active, rng.uniform(0.0, 3.0, active.shape), 0.0)
        ssnow0 = np.where(active & (rng.uniform(size=active.shape) < 0.3), rng.uniform(0, 50, active.shape), 0.0)
        sliq0 = 0.05 * ssnow0
        cap = [np.full(active.shape, 0.35), np.full(active.shape, 0.2)]
        si = [np.where(active, rng.uniform(0, 0.3, active.shape), 0.0) for _ in range(2)]
        tk = fld(g["corr_temperature"])

        # --- current path: process functions one after the other
        layers = [pet.Layer(landuse=raw(g[f"landuse{n}"]), lai=raw(g[f"leaf_area_index{n}"]),
                            coverage=fld(g[f"vegetation_coverage_degree{n}"]), rsc=fld(g[f"SurfaceCanopyResistance{n}"]),
                            rsi=fld(g[f"SurfaceIntercepResistance{n}"]), z0=fld(g[f"RoughnessLength{n}"])) for n in (1, 2)]
        ins = [fld(g[k]) for k in ("global_radiation", "sunshine_duration", "radiation_correction", "air_humidity",
                                   "wind_speed", "albedo")]
        rest = dict(dhm=fld(g["[elevation_model]"]), rse=fld(g["SurfaceEvaporationResistance"]), layers=layers,
                    k_lambert=fld(kl), lai_scale=fld(ls), lake=raw(g["lake_codes"]), pond=raw(g["ponding_storage"]))
        veg = [(g[f"landuse{n}"] != pet.NODATA) & active for n in (1, 2)]
        snow_mask = ssnow0 > 0.0
        no_snow = msk(lakes_on & (g["lake_codes"] > 0) & (g["ponding_storage"] > 0))
        ice = gtx.zeros(dom, dtype=float64)
        steps_per_day = 1440.0 / step

        def current():
            state = pet.State(ssd_mean=raw(g["SSDTemp"]))
            r = pet.potential_et(site, now.day, now.month, now.hour, step, params, tk, tk, *ins, state=state,
                                 active=active, domain=dom, backend=b, **rest)
            ic = icp.interception(fld(prec), [icp.Layer(lai=raw(g[f"leaf_area_index{n}"]),
                                                        coverage=fld(g[f"vegetation_coverage_degree{n}"]),
                                                        capacity=fld(cap[n - 1]), etp=r["etp"][n], storage=fld(si[n - 1]),
                                                        has_vegetation=msk(veg[n - 1])) for n in (1, 2)],
                                  r["etp"][0], ice, dom, backend=b)
            ss, _, ev, evrs = sev.conventional(fld(ssnow0), [ic[1]["etp"], ic[2]["etp"]],
                                               [msk(veg[0]), msk(veg[1])], r["evap"], msk(snow_mask), dom, backend=b)
            out = tuple(gtx.zeros(dom, dtype=float64) for _ in range(3))
            with_backend(snow_step, b)(tk, ic[0]["q"], evrs, fld(ssnow0), fld(sliq0), no_snow, SNOW["t0r"],
                                      SNOW["transition"], SNOW["t0"], SNOW["c0"], SNOW["cwh"], steps_per_day,
                                      out=out, offset_provider={})
            return dict(etp=r["etp"][0].asnumpy(), throughfall=ic[0]["q"].asnumpy(), evrs=evrs.asnumpy(),
                        evap=ev.asnumpy(), ssnow=out[0].asnumpy(), sliq=out[1].asnumpy(), outflow=out[2].asnumpy())

        # --- fused path: one program call; all fields allocated once
        z = lambda: gtx.zeros(dom, dtype=float64)  # noqa: E731
        bottom = [msk(np.full(active.shape, False) | (g["leaf_area_index2"] == pet.NODATA)), msk(np.full(active.shape, True))]
        fins = dict(
            rg=ins[0], ssd=ins[1], ssd_mean=raw(g["SSDTemp"]), exkorr=ins[2], albedo=ins[5], lake=rest["lake"],
            pond=rest["pond"], use_lake=params.lakes, cr0=params.cr[0], cr1=params.cr[1], cr2=params.cr[2],
            cr3=params.cr[3], rgmax=1260.0 * step / 60.0 * dl.bright, bright=True, albedo_ow=params.albedo_open_water,
            tm=tk, dhm=rest["dhm"], rlf=ins[3], u=ins[4], rse=rest["rse"], k_lambert=rest["k_lambert"],
            lai_scale=rest["lai_scale"],
            pond_cell=msk(params.lakes & (g["lake_codes"] > 0) & (g["ponding_storage"] > 0)),
            **{f"{k}{n}": v for n, lay in ((1, layers[0]), (2, layers[1]))
               for k, v in (("lai", lay.lai), ("vbg", lay.coverage), ("rsc", lay.rsc), ("rsi", lay.rsi), ("z0", lay.z0))},
            veg1=msk(veg[0]), bottom1=bottom[0], veg2=msk(veg[1]), bottom2=bottom[1],
            t_first=True, f_first=False, night=False, anteil=dl.bright, sigma=step / 60 * 0.0000000567,
            vv_factor=dl.bright * step * 60.0 * 0.216876, z0_ow=params.z0_open_water,
            prec=fld(prec), cap1=fld(cap[0]), cap2=fld(cap[1]), si1=fld(si[0]), si2=fld(si[1]), ice=ice,
            ssnow=fld(ssnow0), sliq=fld(sliq0), no_snow=no_snow, snow_veg1=msk(snow_mask & veg[0]),
            snow_veg2=msk(snow_mask & veg[1]), snow=msk(snow_mask), steps_per_day=steps_per_day, zero=z(), **SNOW)
        names = ("rk relrg ssd_used e1 ev1 ei1 ra1 rr1 ra_l1 rr_l1 etp0a eip0a lai_kum2 e2 ev2 ei2 ra2 rr2 etp0b eip0b "
                 "iei1 iq1 isi1 ietp1 ietp0 iei0 isi0 iq0 iei2 iq2 isi2 ietp2 ietp0b iei0b isi0b iq0b s1 d1 v1 s2 d2 v2 "
                 "rg_o rgex cloudiness etp1 eip1 evap1 etp2 eip2 evap etp_total throughfall ssnow_o sliq_o outflow "
                 "evap_o evrs").split()
        bufs = {k: z() for k in names}
        prog = hydro_step.with_backend(b)
        allargs = {**fins, **bufs}

        def fused():
            prog(**allargs, offset_provider={})

        def lean():
            hydro_step_lean(b, allargs)

        ref = current()

        def check():
            got = dict(etp=bufs["etp_total"], throughfall=bufs["throughfall"], evrs=bufs["evrs"], evap=bufs["evap_o"],
                       ssnow=bufs["ssnow_o"], sliq=bufs["sliq_o"], outflow=bufs["outflow"])
            return max(float(np.abs(ref[k] - got[k].asnumpy())[active].max()) for k in ref)

        lean()
        worst_lean = check()
        tc = timed(current, args.repeat)
        tl = timed(lean, args.repeat)
        line = (f"{members} member(s), {active.sum()} active cells ({nrows} x {ncols}), {args.backend}: "
                f"current {tc:7.2f} ms, lean calls {tl:7.2f} ms (max |diff| {worst_lean:.3g})")
        if not args.no_program:
            for k in bufs:
                bufs[k] = z()
            allargs.update(bufs)
            t0 = time.perf_counter()
            fused()
            t_compile = time.perf_counter() - t0
            worst = check()
            tf = timed(fused, args.repeat)
            line += f", one program {tf:7.2f} ms (max |diff| {worst:.3g}; compiled in {t_compile:.0f} s)"
        print(line + f"; per member: current {tc / members:.2f}, lean {tl / members:.2f} ms", flush=True)


def timed(fn, repeat, warmup=2):
    for _ in range(warmup):
        fn()
    t = time.perf_counter()
    for _ in range(repeat):
        fn()
    return 1000.0 * (time.perf_counter() - t) / repeat


if __name__ == "__main__":
    main()
