#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# optional: Caches in einen lokalen Ordner mounten
mkdir -p .cache/colabfold .cache/mpl

run_one() {
  local RUN_DIR="$1"
  echo "=== RUN START: ${RUN_DIR} ==="

 python3 scripts/universal/06_equilibrate.py ${RUN_DIR}

  echo "=== RUN DONE : ${RUN_DIR} ==="
}

# Hier nur die RUN_DIRs auflisten – Reihenfolge == Ausführungsreihenfolge
RUNS=(
simulations/solvia_1363_run_1
simulations/solvia_1564_run_1
simulations/solvia_1587_run_1
simulations/solvia_1663_run_1
simulations/solvia_1680_run_1
simulations/solvia_1684_run_1
simulations/solvia_1743_run_1
simulations/solvia_1844_run_1
simulations/solvia_1941_run_1
simulations/solvia_1952_run_1
simulations/solvia_1962_run_1
simulations/solvia_2012_run_1
)

for d in "${RUNS[@]}"; do
  run_one "$d"
done