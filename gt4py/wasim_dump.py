"""Read grid dumps written by WaSiM's debug hook (WASIM_DUMP_GRIDS / WASIM_DUMP_DIR).

The hook (wasim-eth src/modapuzr.cc, griddump) appends one record per call to
<dir>/<point>__<grid>.f64: float64, rows 1..nrows, columns 1..ncols, row-major, and for stacks
levels+1 values per cell. <dir>/manifest.txt lists "point grid nrows ncols levels".
Points: start, after_meteo, after_expokorr, after_etp, after_irrigation, after_interception,
after_snow1, after_snowevap, after_siltingup, after_unsatzon, after_surfrout, after_etreal,
after_gwflow (these four once per sub step), after_snow2, after_post_run, after_routing.
"""
import os
import re

import numpy as np


class Dump:
    def __init__(self, directory):
        self.dir = directory
        self.shape = {}
        self._cache = {}
        with open(os.path.join(directory, "manifest.txt")) as f:
            for line in f:
                point, grid, nrows, ncols, levels = line.split()
                self.shape[(point, grid)] = (int(nrows), int(ncols), int(levels))

    def __call__(self, point, grid):
        """All records of a grid at a point: (steps, rows, cols) or (steps, rows, cols, levels+1)."""
        if (point, grid) in self._cache:
            return self._cache[(point, grid)]
        nrows, ncols, levels = self.shape[(point, grid)]
        data = np.memmap(os.path.join(self.dir, f"{point}__{grid}.f64"), dtype="f8", mode="r")
        data = data.reshape(-1, nrows, ncols, levels + 1) if levels else data.reshape(-1, nrows, ncols)
        self._cache[(point, grid)] = data
        return data

    def has(self, point, grid):
        return (point, grid) in self.shape


def landuse_table(control_file, key):
    """{landuse code: value} for one parameter of [landuse_table] in a (filled) control file."""
    text = open(control_file, encoding="utf-8", errors="replace").read()
    section = text.split("[landuse_table]", 1)[1].split("\n[", 1)[0]
    values = {}
    for m in re.finditer(r"^\s*(\d+)\s+\S+\s*\{(.*?)\}", section, re.M | re.S):
        v = re.search(r"\b" + key + r"\s*=\s*([-+0-9.eE]+)", m.group(2))
        if v:
            values[int(m.group(1))] = float(v.group(1))
    return values
