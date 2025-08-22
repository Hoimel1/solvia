# 🧬 SOLVIA Produktions-Pipeline für Einzelpeptide – Schritt 1: ColabFold (Container)


Dieses Dokument erklärt Schritt für Schritt, wie Du ColabFold im Container ausführst. Die Pipeline ist modular: 1) ColabFold → 2) martinize2 → 3) INSANE → 4) GROMACS.

## Voraussetzungen
- Docker und Docker Compose sind installiert.
- NVIDIA GPU + Treiber verfügbar (für GPU-Beschleunigung).
- Dieses Repository ist lokal geklont.

Schnelltests:
```bash
docker --version
docker compose version
docker run --rm --gpus all nvidia/cuda:12.2.0-base nvidia-smi
```

## Einmalige Vorbereitung
1) Optional: Docker ohne sudo verwenden
```bash
sudo usermod -aG docker $USER
newgrp docker
docker ps
```

2) Eigenes ColabFold-Image bauen (aus diesem Repo)
```bash
docker compose build colabfold
```

3) Cache-Verzeichnis für ColabFold-Datenbanken (auf dem Host) anlegen
```bash
mkdir -p colabfold_cache/colabfold_envdb_202308
```

## 1) Run-Verzeichnis anlegen
```bash
python scripts/universal/01_setup_run.py data/raw/fasta_split/SOLVIA_1.fasta
```
Beispielstruktur:
```
simulations/solvia_1_run_1/
├── input/peptide.fasta
├── metadata.yaml
└── colabfold/
```

## 2) ColabFold im Container starten
```bash
RUN_DIR="simulations/solvia_1_run_1"

docker compose run --rm --gpus all \
  colabfold \
  /work/${RUN_DIR}/input/peptide.fasta /work/${RUN_DIR}/colabfold \
  --num-seeds 3 \
  --num-models 5 \
  --msa-mode mmseqs2_uniref_env \
  --pair-mode unpaired_paired
```
Hinweise:
- Ergebnisse liegen in `${RUN_DIR}/colabfold/` (im Host-Repo sichtbar).
- Die Variable `COLABFOLD_DB` ist im Compose-Service gesetzt und zeigt im Container auf `/cache/colabfold_envdb_202308` (gemappt auf `./colabfold_cache`).

## 3) Bestes Modell auswählen und protokollieren
```bash
python scripts/universal/02_run_colabfold.py ${RUN_DIR}
```
Ergebnis:
```
${RUN_DIR}/colabfold/model_selection.yaml
```
Dort steht, welches Modell gewählt wurde (inkl. pLDDT).

---
Nächster Schritt: martinize2 im Container, basierend auf dem ausgewählten PDB.

# 🧬 SOLVIA-Produktions-Pipeline für Batch-Peptide – Schritt 1: ColabFold (Container)

## Ziel
Mehrere FASTA-Dateien automatisch durch ColabFold schicken und die besten Modelle auswählen.

## Voraussetzung
- Deine FASTA-Dateien liegen z. B. unter `data/raw/fasta_split/`.
- Dateinamen-Endung: `.fasta` (ein Header pro Datei, Peptid-ID in `>`-Zeile).

## 1) Runs für alle FASTA-Dateien anlegen
```bash
for fasta in data/raw/fasta_split/*.fasta; do
  echo "Setup für $(basename "$fasta")"
  python scripts/universal/01_setup_run.py "$fasta"
done
```

## 2) ColabFold für alle Runs starten
Variante A – sequentiell (einfach, stabil):
```bash
for run_dir in simulations/*_run_*; do
  echo "ColabFold für $run_dir"
  docker compose run --rm --gpus all \
    colabfold \
    "/work/${run_dir}/input/peptide.fasta" "/work/${run_dir}/colabfold" \
    --num-seeds 3 --num-models 5 --msa-mode mmseqs2_uniref_env --pair-mode unpaired_paired
done
```

Variante B – parallel (z. B. 2 Jobs gleichzeitig), erfordert GNU parallel:
```bash
ls -d simulations/*_run_* | parallel -j 2 \
  'docker compose run --rm --gpus all colabfold \
   "/work/{}/input/peptide.fasta" "/work/{}/colabfold" \
   --num-seeds 3 --num-models 5 --msa-mode mmseqs2_uniref_env --pair-mode unpaired_paired'
```

## 3) Bestes Modell je Run auswählen
```bash
for run_dir in simulations/*_run_*; do
  echo "Select best model für $run_dir"
  python scripts/universal/02_run_colabfold.py "$run_dir"
done
```

## Ergebnisse
- Pro Run: `simulations/<id>_run_X/colabfold/model_selection.yaml` enthält pLDDT und Pfad zum gewählten PDB.
- Logs: `simulations/<id>_run_X/logs/colabfold.log` (falls vorhanden).