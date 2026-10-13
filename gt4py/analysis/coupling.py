#!/usr/bin/env python3
"""Static analysis of how WaSiM's modules share grids.

WaSiM keeps all grids in one table (modelobject::gridinfo). A module creates its own grids
with initialize_grid(..., identification, ...) and gets grids of other modules by name:
    X = gridinfo[getindex(_Identification_)].pointer
This script finds, per source file:
  - created grids (identification -> pointer variable)
  - acquired grids (identification -> pointer variable), and whether the module writes them
and combines both into a list of cross-module writes: grids written by a module other
than the one that created them. These are the "entanglements" a port has to make explicit.

The analysis is regex-based and approximate: identifiers built at run time (sprintf into a
name variable) are resolved when the sprintf is in the same file, pointer aliases and
writes through helper functions are not followed. Check findings in the code.

usage: coupling.py WASIM_SRC_DIR [--json out.json] [--markdown out.md]
"""
import argparse
import json
import os
import re
from collections import defaultdict

MODULE_FILES = {
    "meteoint.cc": "meteo interpolation + precipitation correction",
    "expokorr.cc": "radiation/temperature correction",
    "etpetr.cc": "evapotranspiration",
    "snowmod.cc": "snow and glaciers",
    "snow_canopy.cc": "snow canopy interception",
    "interzep.cc": "interception",
    "irrigate.cc": "irrigation",
    "siltingup.cc": "silting up",
    "surfrout.cc": "surface routing",
    "unsatzon.cc": "unsaturated zone (Richards)",
    "gwflow.cc": "groundwater",
    "abfluss.cc": "routing, lakes, reservoirs",
    "dataassim.cc": "data assimilation",
    "defltuzr.cc": "initialization (zones, tables, MPI)",
    "modapuzr.cc": "application (time loop, wiring)",
    "teilmod.cc": "submodel base, statistics",
}


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def read_defines(src):
    defines = {}
    for name in ("modconst.h", "wasim.h"):
        p = os.path.join(src, name)
        if os.path.exists(p):
            for m in re.finditer(r'^#define\s+(_\w+_)\s+"([^"]*)"', open(p, encoding="utf-8", errors="replace").read(), re.M):
                defines[m.group(1)] = m.group(2)
    return defines


def resolve_idents(arg, text, defines):
    """Identifiers (#define names) mentioned in a getindex/initialize_grid argument."""
    found = re.findall(r"\b(_\w+_)\b", arg)
    found = [f for f in found if f in defines and f != "_ungueltiger_Name_"]
    if found:
        return found
    var = arg.strip().split("[")[0].strip()
    if re.fullmatch(r"\w+", var):
        # name variable filled by sprintf(var, "...", _Define_, ...) in the same file
        for m in re.finditer(r"sprintf\s*\(\s*" + re.escape(var) + r"\b[^;]*?\b(_\w+_)\b", text):
            if m.group(1) in defines:
                return [m.group(1)]
    return []


def write_count(text, var):
    pat = re.compile(r"(?<![\w.>])" + re.escape(var) + r"\s*((\[[^\]\[]*\])+)\s*(=|\+=|-=|\*=|/=)(?!=)")
    return sum(1 for m in pat.finditer(text) if m.group(1).count("[") >= 2)


def element_access_count(text, var):
    pat = re.compile(r"(?<![\w.>])" + re.escape(var) + r"\s*(\[[^\]\[]*\]){2,}")
    return len(pat.findall(text))


def analyse(src):
    defines = read_defines(src)
    result = {}
    acq_re = re.compile(
        r"(\w+)\s*(?:\[[^\]=;]*\])?\s*=\s*gridinfo\s*\[\s*getindex(?:_first_of_\d+)?\s*\(([^;]*?)\)\s*\]\s*\.\s*(pointer|stackpointer)")
    idx_re = re.compile(r"(\w+)\s*=\s*getindex(?:_first_of_\d+)?\s*\(([^;]*?)\)\s*;")
    via_idx_re = re.compile(r"(\w+)\s*(?:\[[^\]=;]*\])?\s*=\s*gridinfo\s*\[\s*(\w+)\s*\]\s*\.\s*(pointer|stackpointer)")
    init_re = re.compile(r"initialize_grid\s*\(\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),")
    for fname, label in MODULE_FILES.items():
        path = os.path.join(src, fname)
        if not os.path.exists(path):
            continue
        text = strip_comments(open(path, encoding="utf-8", errors="replace").read())
        created = defaultdict(set)
        for m in init_re.finditer(text):
            grid_ptr, stack_ptr, ident = m.group(3).strip(), m.group(4).strip(), m.group(6)
            ptr = grid_ptr if grid_ptr != "NULL" else stack_ptr
            ptr = re.sub(r"\[.*", "", ptr).strip("&* ")
            for d in resolve_idents(ident, text, defines):
                created[d].add(ptr)
        acquired = defaultdict(set)
        for m in acq_re.finditer(text):
            for d in resolve_idents(m.group(2), text, defines):
                acquired[d].add(m.group(1))
        index_vars = defaultdict(list)
        for m in idx_re.finditer(text):
            index_vars[m.group(1)].append(m.group(2))
        for m in via_idx_re.finditer(text):
            for arg in index_vars.get(m.group(2), []):
                for d in resolve_idents(arg, text, defines):
                    acquired[d].add(m.group(1))
        entry = {"label": label, "created": {}, "acquired": {}}
        for d, ptrs in sorted(created.items()):
            entry["created"][d] = {"name": defines[d], "pointers": sorted(ptrs)}
        for d, ptrs in sorted(acquired.items()):
            ptrs = sorted(p for p in ptrs if p not in ("pointer",))
            writes = sum(write_count(text, p) for p in ptrs)
            uses = sum(element_access_count(text, p) for p in ptrs)
            entry["acquired"][d] = {"name": defines[d], "pointers": ptrs, "element_writes": writes, "element_uses": uses}
        result[fname] = entry
    return result, defines


def cross_writes(result):
    creator = defaultdict(set)
    for f, e in result.items():
        for d in e["created"]:
            creator[d].add(f)
    rows = []
    for f, e in result.items():
        for d, a in e["acquired"].items():
            if a["element_writes"] > 0 and f not in creator.get(d, set()):
                rows.append({"grid": d, "name": a["name"], "created_by": sorted(creator.get(d, [])) or ["(input grid)"],
                             "written_by": f, "writes": a["element_writes"], "pointers": a["pointers"]})
    return sorted(rows, key=lambda r: (r["written_by"], r["grid"]))


def readers(result):
    used = defaultdict(set)
    for f, e in result.items():
        for d, a in e["acquired"].items():
            if a["element_uses"] > 0:
                used[d].add(f)
    return used


def to_markdown(result, rows):
    out = ["# Grids shared between WaSiM modules (generated by coupling.py)", ""]
    out.append("## Cross-module writes")
    out.append("")
    out.append("Grids a module writes although another module (or the input) owns them.")
    out.append("")
    out.append("| Grid | Internal name | Created by | Written by | Element writes | Pointer(s) |")
    out.append("| --- | --- | --- | --- | --- | --- |")
    for r in rows:
        out.append(f"| `{r['grid']}` | {r['name']} | {', '.join(r['created_by'])} | {r['written_by']} | {r['writes']} | {', '.join(r['pointers'])} |")
    out.append("")
    out.append("## Per module")
    for f, e in result.items():
        out.append("")
        out.append(f"### {f} ({e['label']})")
        out.append("")
        out.append("Creates: " + (", ".join(f"`{d}`" for d in e["created"]) or "-"))
        out.append("")
        out.append("Uses from the grid table: " + (", ".join(
            f"`{d}`{' (writes)' if a['element_writes'] else ''}" for d, a in e["acquired"].items()) or "-"))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--json")
    ap.add_argument("--markdown")
    args = ap.parse_args()
    result, defines = analyse(args.src)
    rows = cross_writes(result)
    if args.json:
        json.dump({"modules": result, "cross_writes": rows}, open(args.json, "w"), indent=1)
    if args.markdown:
        open(args.markdown, "w").write(to_markdown(result, rows))
    for f, e in result.items():
        print(f"{f:16s} creates {len(e['created']):3d}  acquires {len(e['acquired']):3d}  writes {sum(1 for a in e['acquired'].values() if a['element_writes'])}")
    print(f"{len(rows)} cross-module writes")


if __name__ == "__main__":
    main()
