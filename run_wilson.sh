#!/bin/bash

set -e

# Simulation parameters
L=${1:-32}
BETA=${2:-2.0}
N_SWEEP=${3:-1e6}

EXEC="./u1gauge2d"
OUTDIR="results"

mkdir -p "$OUTDIR"

BASE_SEED=12345
JOB=0

# Wilson loop dimensions to measure
R_VALUES=(1 2 3 4 5 6 7 8)
T_VALUES=(1 2 3 4 5 6 7 8)

for R in "${R_VALUES[@]}"; do
    for T in "${T_VALUES[@]}"; do

        JOB=$((JOB + 1))
        SEED=$((BASE_SEED + JOB))

        FILENAME="Wc_R${R}_T${T}"

        echo "Running: R=$R T=$T"

        "$EXEC" \
            "$L" \
            "$BETA" \
            "$N_SWEEP" \
            "$FILENAME" \
            "$R" \
            "$T" \
            "$SEED"

    done
done

echo "Campaign completed!"