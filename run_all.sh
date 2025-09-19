#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# optional: Caches in einen lokalen Ordner mounten
mkdir -p .cache/colabfold .cache/mpl

run_one() {
  local RUN_DIR="$1"
  echo "=== RUN START: ${RUN_DIR} ==="

  python3 scripts/universal/05_insert_peptides.py "${RUN_DIR}" \
    --n-peptides 1 --orientation parallel
  python3 scripts/universal/06_equilibrate.py "${RUN_DIR}"

  echo "=== RUN DONE : ${RUN_DIR} ==="
}

# Hier nur die RUN_DIRs auflisten – Reihenfolge == Ausführungsreihenfolge
RUNS=(
  simulations/solvia_8_run_2
  simulations/solvia_14_run_2
  simulations/solvia_32_run_2
  simulations/solvia_68_run_2
  simulations/solvia_126_run_2
  simulations/solvia_215_run_2
  simulations/solvia_482_run_2
  simulations/solvia_490_run_2
  simulations/solvia_515_run_2
  simulations/solvia_524_run_2
  simulations/solvia_527_run_2
  simulations/solvia_617_run_2
  simulations/solvia_624_run_2
  simulations/solvia_850_run_2
  simulations/solvia_858_run_2
  simulations/solvia_941_run_2
  simulations/solvia_974_run_2
  simulations/solvia_1023_run_2
  simulations/solvia_1045_run_2
)

for d in "${RUNS[@]}"; do
  run_one "$d"
done