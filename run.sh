#!/usr/bin/env bash
# ============================================================================
# run.sh — Nihar's one-shot runner for the 5 "Sums of Three Positive Cubes"
# research extensions.
#
# What it does, in order, for each ext1..ext5:
#   1. compiles it fresh with g++ -O3 -std=c++17 -pthread
#   2. runs it (each writes its own results_extN_*.txt)
#   3. deletes the compiled binary and any downloaded OEIS b-file scratch
#      data afterwards, leaving ONLY the .cpp sources and the results_*.txt
#      output files behind — no build artifacts, no clutter.
#
# Usage:
#   chmod +x run.sh
#   ./run.sh
#
# Works on Linux and macOS (Intel or Apple Silicon) out of the box. On
# Windows, run it from Git Bash / WSL, or just run the five g++ commands
# printed below manually from PowerShell/cmd.
# ============================================================================
set -e

SOURCES=(
  "ext1_extreme_value_cube_multiplicity.cpp"
  "ext2_equidistribution_explicit_bounds.cpp"
  "ext3_elliptic_descent_families.cpp"
  "ext4_chao_richness_estimation.cpp"
  "ext5_parallel_asymptotic_regression.cpp"
)

echo "=============================================================="
echo " Building and running all 5 extensions — results_*.txt only"
echo " will be left behind when this finishes."
echo "=============================================================="

for src in "${SOURCES[@]}"; do
    if [ ! -f "$src" ]; then
        echo "!! Skipping $src — file not found in this directory."
        continue
    fi
    bin="${src%.cpp}"
    echo ""
    echo "--- $src ---"
    echo "  compiling..."
    g++ -O3 -std=c++17 -pthread -o "$bin" "$src"

    echo "  running..."
    ./"$bin" "$@"

    echo "  cleaning up binary and any scratch download files..."
    rm -f "$bin"
    rm -f A003072_b_file*.txt A003325_b_file*.txt
done

echo ""
echo "=============================================================="
echo " Done. Results:"
echo "=============================================================="
ls -la results_*.txt 2>/dev/null || echo "(no results_*.txt found — check the compile/run logs above)"