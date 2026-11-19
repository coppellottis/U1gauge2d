
#!/bin/bash

SWEEPS=100000
MAX_JOBS=4

CONFIGS=(
    "16 2"
    "24 4"
    "28 6"
    "32 8"
    "40 12"
)

mkdir -p results logs

for CONFIG in "${CONFIGS[@]}"; do

    read -r L BETA <<< "$CONFIG"

    echo "Starting L=$L beta=$BETA"

    R=1
    T=1
    # for R in {1..6}; do
    #     for T in {1..6}; do

            FILE="Q_W_L${L}_beta${BETA}"
            LOG="logs/L${L}_b${BETA}.log"

            # Limita il numero di processi contemporanei
            while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
                wait -n
            done

            # Seed differente per ogni processo
            SEED=$(od -An -N4 -tu4 /dev/urandom | tr -d ' ')

            ./u1gauge2d "$L" "$BETA" "$SWEEPS" \
                "$FILE" "$R" "$T" "$SEED" > "$LOG" 2>&1 &

    #     done
    # done
done

# Attende il completamento di tutti i processi
wait

echo "All simulations completed!"


# #!/bin/bash

# set -e

# # Simulation parameters
# L=${1:-32}
# BETA=${2:-2.0}
# N_SWEEP=${3:-1e6}

# EXEC="./u1gauge2d"
# OUTDIR="results"

# mkdir -p "$OUTDIR"

# BASE_SEED=12345
# JOB=0

# # Wilson loop dimensions to measure
# R_VALUES=(1 2 3 4 5 6 7 8)
# T_VALUES=(1 2 3 4 5 6 7 8)

# for R in "${R_VALUES[@]}"; do
#     for T in "${T_VALUES[@]}"; do

#         JOB=$((JOB + 1))
#         SEED=$((BASE_SEED + JOB))

#         FILENAME="Wc_R${R}_T${T}"

#         echo "Running: R=$R T=$T"

#         "$EXEC" \
#             "$L" \
#             "$BETA" \
#             "$N_SWEEP" \
#             "$FILENAME" \
#             "$R" \
#             "$T" \
#             "$SEED"

#     done
# done

# echo "Campaign completed!"