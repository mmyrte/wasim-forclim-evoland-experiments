#!/bin/bash
# All measurements of the CPU spike, one after the other (no other load on the machine).
# usage: run_all.sh WORK_DIR WASIM_GT_DIR   (WORK_DIR: cases and outputs; dumps of the golden runs below)
# Compiles first (warm caches), then times. Results: WORK_DIR/results.txt
set -u
W=$1
PY=$2/.venv/bin/python
S=$(cd "$(dirname "$0")" && pwd)
R=$W/results.txt
export DACE_compiler_cpu_args="-std=c++17 -fPIC -O3 -fopenmp"  # DaCe's default adds -march=native -ffast-math
export GT4PY_BUILD_CACHE_LIFETIME=persistent GT4PY_BUILD_CACHE_DIR=$W/gt4py_cache  # default: recompile per process
cd "$W" || exit 1
: > "$R"
log() { echo "$*" | tee -a "$R"; }

c++ -O3 -fopenmp -o richards_ref "$S/richards_ref.cpp"
for case in doubs rhone; do
  [ -d case_$case ] || "$PY" "$S/proxy_case.py" case_$case --cells $([ $case = doubs ] && echo 6022 || echo 21584) \
      --seed $([ $case = doubs ] && echo 1 || echo 2)
  log "== Richards proxy, $case"
  for th in 1 4; do
    mkdir -p out_${case}_ref$th
    log "C++ reference, $th thread(s): $(OMP_NUM_THREADS=$th ./richards_ref case_$case out_${case}_ref$th 3)"
  done
  ref=out_${case}_ref1
  for th in 1 4; do
    log "DaCe, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/richards_dace.py" case_$case out_${case}_dace$th 2>/dev/null | tail -1)"
    log "  $("$PY" "$S/compare.py" case_$case $ref out_${case}_dace$th | head -1)"
    log "DaCe, time loop inside, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/richards_dace.py" case_$case out_${case}_daceL$th --loop 2>/dev/null | tail -1)"
  done
  for th in 1 4; do
    log "cartesian lockstep, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/richards_cartesian.py" case_$case out_${case}_cart$th 2>/dev/null | tail -1)"
    log "  $("$PY" "$S/compare.py" case_$case $ref out_${case}_cart$th | head -1)"
    log "next lockstep, stages called one by one, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/richards_next.py" case_$case out_${case}_next$th --split 2>/dev/null | tail -1)"
    log "  $("$PY" "$S/compare.py" case_$case $ref out_${case}_next$th | head -1)"
  done
  # the same stages composed into one program: much slower per tick, so 4 steps only
  mkdir -p out_${case}_ref1_4 && OMP_NUM_THREADS=1 ./richards_ref case_$case out_${case}_ref1_4 1 4 > /dev/null
  log "next lockstep, one program per tick, 4 steps, 1 thread: $(OMP_NUM_THREADS=1 "$PY" "$S/richards_next.py" case_$case out_${case}_nextP --steps 4 2>/dev/null | tail -1)"
  log "  $("$PY" "$S/compare.py" case_$case out_${case}_ref1_4 out_${case}_nextP | head -1)"
  log "cartesian fixed N=1 I=8, 4 threads: $(OMP_NUM_THREADS=4 "$PY" "$S/richards_cartesian.py" case_$case out_${case}_cartF --mode fixed --substeps 1 --iters 8 2>/dev/null | tail -1)"
  for cfg in "1 4" "1 6" "1 8" "1 12" "1 16" "2 8"; do
    set -- $cfg
    o=out_${case}_fixed_N$1_I$2
    log "next fixed N=$1 I=$2: $(OMP_NUM_THREADS=1 "$PY" "$S/richards_next.py" case_$case $o --mode fixed --substeps $1 --iters $2 --split 2>/dev/null | tail -1)"
    "$PY" "$S/compare.py" case_$case $ref $o | sed 's/^/  /' | tee -a "$R"
  done
done

log "== fused 2-D chain (potential ET, interception, snow evaporation, T-index snow)"
G=$W/../dumps
for th in 1 4; do
  log "Doubs, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/fused_chain.py" $G/bench_doubs $W/../runs/bench_doubs/control.wasim.intern 2021-06-13T01 --members 1 4 16 2>/dev/null | tail -3 | tr '\n' ' ')"
  log "Rhone, $th thread(s): $(OMP_NUM_THREADS=$th "$PY" "$S/fused_chain.py" $G/bench_rhone $W/../runs/bench_rhone/control.wasim.intern 2021-06-13T24 --step 12 --members 1 4 2>/dev/null | tail -2 | tr '\n' ' ')"
done
