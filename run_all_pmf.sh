#!/bin/bash

# Liste deiner Simulationen
runs=(
  simulations/solvia_253_run_1
  simulations/solvia_302_run_1
  simulations/solvia_313_run_1
  simulations/solvia_677_run_1
  simulations/solvia_794_run_1
  simulations/solvia_1051_run_1
  simulations/solvia_1115_run_1
  simulations/solvia_1343_run_1
)

# Für jede GPU einen Run starten
for i in {0..7}; do
  export CUDA_VISIBLE_DEVICES=$i
  run_dir=${runs[$i]}
  echo "▶️  Starte Run $run_dir auf GPU $i"
  
  python3 scripts/universal/08_run_pmf.py "$run_dir" \
    --replicate 1 \
    --tag pmf_midplane_gpu$i > logs/run_gpu$i.log 2>&1 &
done

# Warten bis alle fertig sind
wait
echo "✅ Alle 8 PMF-Runs abgeschlossen!"