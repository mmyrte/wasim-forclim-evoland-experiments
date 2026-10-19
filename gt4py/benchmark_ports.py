#!/usr/bin/env python3
"""Warm per-step timing of wasim-gt ports, measured inside their validation scripts.

Runs a validation script of this directory unchanged, with the port's functions wrapped by a timer.
Only time spent inside the port functions counts (outermost call, so nested port calls are not counted
twice): Python glue and host-device copies of the port are included, dump I/O and the comparison with
WaSiM are not. The first step's calls are warm-up (compilation) and are left out, as are later calls that
compile a new operator variant (longer than 0.5 s and 50 times the median). This complements
benchmark_et.py; see ../validation/performance.md.

usage: benchmark_ports.py PROCESS --steps N -- SCRIPT_ARGS...
  PROCESS: meteo | snow | glacier_melt | permafrost | groundwater
  N: the number of model steps (or groundwater sub steps) the script processes, to group calls by step
Example:
  OMP_NUM_THREADS=1 python3 gt4py/benchmark_ports.py snow --steps 48 -- <dump dir> <run dir>/control.wasim
"""
import argparse
import importlib
import pathlib
import runpy
import statistics
import sys
import time

HERE = pathlib.Path(__file__).parent

# process -> (validation script, {group: (module, [function names])})
PROCESSES = {
    "meteo": ("meteo_interpolation.py", {
        "meteo interpolation (host fits and grids)": ("wasim_gt.processes.meteo_interpolation", [
            "correct_precipitation_profile", "distance_weights", "fit_regression", "interpolate_point",
            "meteo_interpolation", "precipitation_factor", "rotate_coordinates"])}),
    "snow": ("snow_energy_balance.py", {
        "snow part 1, method 7": ("wasim_gt.processes.snow_energy_balance", ["snow_energy_balance"]),
        "snow part 2, slides (host)": ("wasim_gt.processes.snow_slides", ["snow_slides"])}),
    "glacier_melt": ("glacier_melt.py", {
        "snow part 1 on cells and glacier fraction": (
            "wasim_gt.processes.snow_energy_balance", ["snow_energy_balance"]),
        "glacier melt, method 13": ("wasim_gt.processes.glacier_melt", ["glacier_melt"])}),
    "permafrost": ("permafrost.py", {
        "permafrost": ("wasim_gt.processes.permafrost", ["permafrost"])}),
    "groundwater": ("groundwater.py", {
        "groundwater (per sub step)": ("wasim_gt.processes.groundwater", ["groundwater"])}),
}


class Timer:
    def __init__(self):
        self.depth = 0
        self.calls = {}  # group -> list of durations of outermost calls

    def wrap(self, group, fn):
        self.calls.setdefault(group, [])

        def timed(*args, **kwargs):
            outer = self.depth == 0
            self.depth += 1
            start = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                self.depth -= 1
                if outer:
                    self.calls[group].append(time.perf_counter() - start)
        return timed


def main():
    argv = sys.argv[1:]
    own, script_args = (argv[:argv.index("--")], argv[argv.index("--") + 1:]) if "--" in argv else (argv, [])
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("process", choices=sorted(PROCESSES))
    ap.add_argument("--steps", type=int, required=True)
    args = ap.parse_args(own)
    script, groups = PROCESSES[args.process]
    timer = Timer()
    for group, (module, names) in groups.items():
        mod = importlib.import_module(module)
        for name in names:
            setattr(mod, name, timer.wrap(group, getattr(mod, name)))
    sys.path.insert(0, str(HERE))
    sys.argv = [str(HERE / script), *script_args]
    try:
        runpy.run_path(str(HERE / script), run_name="__main__")
    except SystemExit as exc:
        if exc.code not in (0, None):
            print(f"note: the validation script exited with {exc.code} (differences to WaSiM)")
    print(f"\n{args.process}: ms per step, warm (first step left out), {args.steps} steps")
    for group, calls in timer.calls.items():
        if not calls:
            print(f"  {group}: no calls")
            continue
        per_step = len(calls) / args.steps
        warm = calls[int(round(per_step)):] if len(calls) > per_step else calls
        # a call that compiles a new operator variant later in the run (e.g. another station count) is
        # left out: longer than 0.5 s and 50 times the median call
        limit = max(0.5, 50 * statistics.median(warm))
        compiled = [c for c in warm if c > limit]
        warm = [c for c in warm if c <= limit]
        warm_steps = max(1, args.steps - 1)
        print(f"  {group}: {1000 * sum(warm) / warm_steps:.2f} ms per step  "
              f"({len(calls)} calls, {per_step:.1f} per step, median call {1000 * statistics.median(warm):.2f} ms, "
              f"first step {1000 * sum(calls[:int(round(per_step))]):.0f} ms, "
              f"{len(compiled)} later compiling calls left out)")


if __name__ == "__main__":
    main()
