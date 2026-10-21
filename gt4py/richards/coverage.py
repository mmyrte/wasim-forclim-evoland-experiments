"""Measure an explicit, conservative Richards source-line set for the two setups.

python coverage.py WASIM_SRC DOUBS_OBJECTS RHONE_OBJECTS OUTPUT_DIR
Outputs contain source locations and counts, never case values. The scope is a
reviewable line set, not an assertion that every branch in a retained line runs.
Allocation and I/O are excluded; common guards and shared numerical statements
remain. Graph exchange and tracer/heat-only bodies have explicit exclusions.
"""

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parents[2] / "validation"))
from wasim_coverage import function_range, gcov_lines  # noqa: E402

FUNCTIONS = [
    "pre_run",
    "run",
    "post_run",
    "quflow",
    "gradient_restriction",
    "GetMeanTheta",
    "evapotranspiration",
    "hydraulic_heads_iterative",
    "calc_surface_interaction",
    "specific_fluxes",
    "flux_restrictions_opt",
    "new_thetas",
    "groundwater_recharge",
    "CalculateMacroporeInfiltration",
]
# Inclusive exclusions; all other lines within each named function are retained.
EXCLUSIONS = {
    "pre_run": [
        (7767, 7848, "#9 heat/layered snow"),
        (7897, 7897, "tracer statistics"),
        (7899, 7902, "#9 heat statistics"),
        (7905, 7910, "#8 river table and debug"),
    ],
    "run": [
        (8077, 8092, "#8 ordered river-link gather"),
        (8111, 8247, "row allocation; audited separately"),
        (8357, 8472, "ArcVeg and #9 layered snow"),
        (8476, 8716, "tracers, layered snow, heat and debug"),
        (8818, 8988, "legacy diagnostic dump/replay"),
        (9047, 9095, "tracer groundwater mixing"),
        (9130, 9141, "tracer decay/mixing"),
        (9171, 9189, "tracer mixing"),
        (9359, 9367, "debug tracer balance"),
        (9375, 9380, "debug warning"),
        (9393, 9401, "debug tracer balance"),
        (9416, 9422, "debug-only correction"),
        (9458, 9463, "legacy uncoupled baseflow"),
        (9527, 9531, "#13 irrigation"),
        (9674, 9855, "#9 heat/layered snow"),
        (9897, 9901, "unoptimized flux restrictions: >=200 layers or layered snow"),
        (9914, 9920, "debug output"),
        (9924, 9949, "#9 layered snow"),
        (9955, 10430, "diagnostic dump output and #8 river-table update"),
        (10445, 10462, "#9 layered snow"),
        (10579, 10582, "#13 irrigation"),
        (10625, 10666, "tracer transport"),
        (10674, 10686, "tracer mixing"),
        (10711, 10728, "tracer mixing"),
        (10736, 10742, "tracer clearing"),
        (10801, 10888, "#9 layered snow and ArcVeg"),
        (10904, 11158, "#9 heat-derived grids"),
        (11170, 11199, "#13 irrigation diagnostics"),
        (11226, 11318, "debug balance"),
        (11326, 11381, "#9 heat and layered snow writeback"),
        (11411, 11554, "#9 glacier temperature neighbour search"),
        (11557, 11740, "control-cell diagnostics, deallocation and #8 scatter"),
        (11743, 11890, "HRU spin-up"),
        (11937, 11973, "tracer statistics"),
        (12000, 12016, "uncoupled baseflow"),
        (12033, 12115, "tracer statistics"),
        (12117, 12124, "#9 heat statistics"),
        (12206, 12296, "buffered diagnostic output"),
    ],
    "post_run": [
        (13110, 13116, "ArcVeg and #8 river table"),
        (13149, 13152, "tracers"),
        (13158, 13162, "tracers"),
        (13176, 13179, "tracers"),
        (13185, 13189, "tracers"),
        (13225, 13402, "errors, statistics and legacy output"),
    ],
    "gradient_restriction": [(n, n, "#9 layered-snow lookup") for n in (4877, 4902, 4905, 4934, 4937)],
    "hydraulic_heads_iterative": [(5241, 5241, "#9 layered-snow lookup"), (5245, 5245, "#9 layered-snow lookup")],
    "evapotranspiration": [(5179, 5191, "#9 frozen-soil ET")],
    "specific_fluxes": [
        (5480, 5484, "#9 snow conductivity"),
        (5536, 5551, "ArcVeg"),
        (5660, 5660, "clay layer disabled"),
        (5694, 5698, "#9 snow maxratio"),
        (5784, 5796, "ArcVeg and #9 snow interflow"),
        (5801, 5814, "drainage disabled"),
    ],
    "flux_restrictions_opt": [
        (6222, 6222, "#9 SE array"),
        (6226, 6230, "#9 negative snow offset"),
        (6556, 6556, "#9 SE array"),
    ],
    "groundwater_recharge": [(6685, 6701, "uncoupled/external groundwater options")],
    "calc_surface_interaction": [(5437, 5438, "tracers"), (5447, 5448, "tracers")],
}


def regime(function, number):
    if function in ("evapotranspiration", "GetMeanTheta") or 9251 <= number <= 9452:
        return (
            "DroughtAndRoots",
            (
                "First vary rooting depth, dry/wet ET limits and vegetation coverage in copied "
                "controls; then use sustained dry forcing to wilting and a wet-root column."
            ),
        )
    if function == "CalculateMacroporeInfiltration":
        return (
            "MacroporePulse",
            (
                "First lower PMacroThreshold and vary capacity/depth in copied soil tables; then "
                "apply a positive rain pulse to a partly saturated column."
            ),
        )
    if function == "calc_surface_interaction":
        return (
            "RiverExchange",
            (
                "Coordinate #8: vary river depth/colmation and available river flow; use deep "
                "contact levels and limited reinfiltration. Keep river-link ordering explicit."
            ),
        )
    if (
        function in ("new_thetas", "hydraulic_heads_iterative", "groundwater_recharge")
        or 9538 <= number <= 9651
        or 10434 <= number <= 10578
    ):
        return (
            "RisingFallingGroundwater",
            (
                "First vary starting groundwater state and module substeps in copied controls; "
                "use synthetic upward/downward exchange across layer interfaces and the moving "
                "active-depth witness."
            ),
        )
    if function in ("specific_fluxes", "gradient_restriction", "flux_restrictions_opt"):
        return (
            "FluxDirectionAndBounds",
            (
                "First vary rainfall, min step and convergence settings; then synthetic "
                "saturation/wilting, upward and downward gradients, and thaw-depth input from #6."
                " Heat/SE-only cases remain #9."
            ),
        )
    if function in ("pre_run", "post_run") or number >= 11898:
        return (
            "ZonalAndDaily",
            (
                "Vary travel-time classes and reservoir constants in copied controls; add a daily"
                " interval with forcing aggregation defined explicitly. Retain deterministic "
                "zone/cell reduction order."
            ),
        )
    return (
        "SurfacePondAndAdaptiveStep",
        (
            "First vary maximum pond depth, minimum dynamic step and iteration limit in "
            "copied controls; then a saturated, dry, or high-rain column. Daily forcing needs"
            " an explicit aggregation case."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("source", "doubs_objects", "rhone_objects", "output"):
        parser.add_argument(key, type=pathlib.Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = {
        name: gcov_lines(str(args.source), str(obj), "unsatzon.cc")
        for name, obj in [("Doubs", args.doubs_objects), ("Rhone", args.rhone_objects)]
    }
    summary, partition, missing = {}, {}, []
    for name in FUNCTIONS:
        first, last = function_range(str(args.source), "unsatzon.cc", "unsatzonmodellklasse::" + name)
        all_lines = [n for n, (count, _) in data["Doubs"].items() if first <= n <= last and count is not None]
        excluded = {n: reason for lo, hi, reason in EXCLUSIONS.get(name, []) for n in range(lo, hi + 1)}
        selected = [n for n in all_lines if n not in excluded]
        summary[name] = dict(
            first=first,
            last=last,
            whole=len(all_lines),
            in_scope=len(selected),
            doubs=sum(data["Doubs"][n][0] > 0 for n in selected),
            rhone=sum(data["Rhone"][n][0] > 0 for n in selected),
            union=sum(any(data[case][n][0] > 0 for case in data) for n in selected),
        )
        partition[name] = dict(
            in_scope_lines=selected,
            exclusions=[dict(first=lo, last=hi, reason=reason) for lo, hi, reason in EXCLUSIONS.get(name, [])],
        )
        for n in selected:
            if not any(data[case][n][0] > 0 for case in data):
                group, proposal = regime(name, n)
                missing.append(dict(function=name, line=n, group=group, proposal=proposal))
    (args.output / "coverage.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output / "line-set.json").write_text(json.dumps(partition, indent=2) + "\n")
    (args.output / "missing-lines.json").write_text(json.dumps(missing, indent=2) + "\n")
    lines = [
        "# Richards in-scope line coverage",
        "",
        (
            "Conservative numerical/control-flow line set, single-thread no-OpenMP build. "
            "Allocation/I/O, tracers, heat, ArcVeg, unused methods and graph communication "
            "have explicit exclusions in line-set.json. Guards shared with these options "
            "remain. This is line coverage, not branch coverage."
        ),
        "",
        "| Function | Range | Whole executable | In scope | Doubs | Rhone | Union |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for name, s in summary.items():
        lines.append(
            f"| {name} | {s['first']}–{s['last']} | {s['whole']} | {s['in_scope']} "
            f"| {s['doubs']} | {s['rhone']} | {s['union']} |"
        )
    lines += [
        "",
        "## Every unexecuted in-scope line",
        "",
        "| Regime | Function | Source line | Reaching proposal |",
        "| --- | --- | --- | --- |",
    ]
    for item in sorted(missing, key=lambda x: (x["group"], x["line"])):
        lines.append(f"| {item['group']} | {item['function']} | unsatzon.cc:{item['line']} | {item['proposal']} |")
    (args.output / "coverage.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
