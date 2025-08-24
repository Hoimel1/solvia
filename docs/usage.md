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
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
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

4) .env für Dateirechte (empfohlen, damit Container-Dateien Dir gehören)
```bash
printf "UID=$(id -u)\nGID=$(id -g)\n" > .env
cat .env
# docker-compose.yml nutzt user: "${UID}:${GID}"
```

## 1) Run-Verzeichnis anlegen
```bash
python3 scripts/universal/01_setup_run.py data/raw/fasta_split/SOLVIA_1.fasta
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
RUN_DIR="simulations/solvia_12_run_1"

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

---

# 🧪 Schritt 2: Coarse-Graining mit martinize2 (vermouth) – Container

Ziel: Aus der besten ColabFold-Struktur (`model_selection.yaml`) eine Martini 3 CG-Struktur und Topologie erzeugen (`cg_pdb/*.pdb`, `cg_pdb/*.itp`).

Wichtige Voraussetzung:
- Martini 3 Force-Fields sind im Repo unter `force_fields/` vorhanden (z. B. `martini_v3.0.0.itp`).
- Passe einmalig `config/solvia_config.yaml` an, damit Pfade im Container stimmen:
```bash
# Empfohlen für Container: setze force_fields auf /work/force_fields
sed -i 's#^\s*force_fields:.*#  force_fields: "/work/force_fields"#' config/solvia_config.yaml
```

## Image bauen (einmalig oder nach Änderungen)
```bash
# Aus diesem Repo heraus bauen
docker compose build martinize

# Optional: Manuell taggen und pushen (falls kein CI genutzt wird)
# echo <GHCR_PAT> | docker login ghcr.io -u <GITHUB_USER> --password-stdin
# docker build -t ghcr.io/<OWNER>/<REPO>/martinize:latest -f docker/Dockerfile.martinize docker
# docker push ghcr.io/<OWNER>/<REPO>/martinize:latest

# Prüfen, dass DSSP im Container verfügbar ist
docker compose run --rm --entrypoint dssp martinize --version
```

## Pfad- und Input-Check 
Wichtig: Du übergibst nur den `RUN_DIR`. Das Skript liest automatisch `colabfold/model_selection.yaml` und wählt die dort eingetragene PDB (`best_pdb`).

1) Prüfen, ob Auswahl vorhanden ist
```bash
RUN_DIR="simulations/solvia_1_run_1"
[ -f "${RUN_DIR}/colabfold/model_selection.yaml" ] && echo "OK: Auswahl vorhanden" || echo "FEHLT: Bitte Schritt ColabFold/Selektion ausführen"
```

2) Anzeigen, welche PDB verwendet wird (ohne Zusatztools)
```bash
python - <<'PY'
import sys, yaml
run_dir = "simulations/solvia_1_run_1"  # ggf. anpassen
sel = yaml.safe_load(open(f"{run_dir}/colabfold/model_selection.yaml"))
print(sel.get("best_pdb", "kein best_pdb gefunden"))
PY
```

3) Prüfen, ob Force-Fields im Container-Pfad verfügbar sind
```bash
[ -f force_fields/martini_v3.0.0.itp ] && echo "OK: FF vorhanden" || echo "FEHLT: martini_v3.0.0.itp"
grep -n "force_fields:" -n config/solvia_config.yaml
```

## Einzelnes Peptid ausführen
```bash
RUN_DIR="simulations/solvia_1_run_1"

docker compose run --rm \
  --entrypoint python \
  martinize \
  scripts/universal/03_coarse_grain.py ${RUN_DIR}

# Validierung
ls -l "${RUN_DIR}"/cg_pdb
head -n 5 "${RUN_DIR}"/cg_pdb/*_cg.pdb
```
Ergebnis:
```
${RUN_DIR}/cg_pdb/<PEPTID>_cg.pdb
${RUN_DIR}/cg_pdb/<PEPTID>.itp
${RUN_DIR}/cg_pdb/<PEPTID>_complete.top
${RUN_DIR}/logs/martinize2.log
```

Hinweis: Für sehr niedrige pLDDT (<70) setzt das Skript automatisch `--martini3-idp`.

## Batches (mehrere Runs)
```bash
# Sequentiell über alle Runs
echo "Starte martinize2 für alle Runs..."
for run_dir in simulations/*_run_*; do
  if [ -f "$run_dir/colabfold/model_selection.yaml" ]; then
    echo "martinize2: $run_dir"
    # (optional) Rechte korrigieren & alte Outputs säubern
    sudo chown -R $(id -u):$(id -g) "$run_dir" || true
    rm -f "$run_dir"/cg_pdb/* || true
    docker compose run --rm --entrypoint python martinize \
      scripts/universal/03_coarse_grain.py "$run_dir"
    ls -l "$run_dir"/cg_pdb || true
  else
    echo "Überspringe (keine Auswahl): $run_dir"
  fi
done
```
Optional parallel (mit GNU parallel):
```bash
ls -d simulations/*_run_* | parallel -j 2 \
  'test -f {}/colabfold/model_selection.yaml && \
   docker compose run --rm --entrypoint python martinize \
   scripts/universal/03_coarse_grain.py {}'
```

## Troubleshooting (martinize)
- "DSSP executable not found":
  - Stelle sicher, dass Du das Image neu gebaut hast: `docker compose build martinize`
  - Prüfe DSSP im Container: `docker compose run --rm --entrypoint dssp martinize --version`
- "could not select device driver with capabilities: [[gpu]]":
  - Gilt nur für GPU-Container (z. B. ColabFold/GROMACS). Installiere NVIDIA Container Toolkit und starte Docker neu (siehe Hinweise im oberen Teil dieser Datei) oder teste ohne `--gpus all` für nicht-GPU-Container wie `martinize`.

---

# 🧪 Schritt 3: Membranaufbau mit INSANE – Container

Ziel: Eine asymmetrische RBC-Membran gemäß `config/solvia_config.yaml` erstellen. Ergebnis: `membrane_template.gro`, `membrane_template.top`, `composition.yaml`.

## Einmalig vorbereiten (INSANE)
```bash
# Image bauen
docker compose build insane

# Pfad für globale Templates im Container korrekt setzen (optional, nur für --global-template)
sed -i 's#^\s*templates:.*#  templates: "/work/data/templates"#' config/solvia_config.yaml

# Verzeichnis anlegen (nur falls global genutzt)
mkdir -p data/templates/membrane
```

## 1) Membran im Run-Verzeichnis erstellen (empfohlen)
```bash
RUN_DIR="simulations/solvia_1_run_1"

docker compose run --rm \
  --entrypoint python \
  insane \
  scripts/universal/04_build_membrane.py "${RUN_DIR}"

# Validierung
ls -l "${RUN_DIR}"/membrane_template
head -n 5 "${RUN_DIR}"/membrane_template/membrane_template.gro
grep -n "\\[ molecules \\]" "${RUN_DIR}"/membrane_template/membrane_template.top | cat
cat "${RUN_DIR}"/membrane_template/composition.yaml | cat
```

Ausgaben liegen in `${RUN_DIR}/membrane_template/`. Logs: `${RUN_DIR}/logs/insane_membrane.log`.

## 2) Globale Membranvorlage (optional)
```bash
docker compose run --rm --entrypoint python insane \
  scripts/universal/04_build_membrane.py --global-template

# Dateien unter: data/templates/membrane/
ls -l data/templates/membrane
```

## 3) Batches (mehrere Runs)
```bash
echo "Starte INSANE für alle Runs..."
for run_dir in simulations/*_run_*; do
  echo "INSANE: $run_dir"
  docker compose run --rm --entrypoint python insane \
    scripts/universal/04_build_membrane.py "$run_dir"
  ls -l "$run_dir"/membrane_template || true
done
```

## Troubleshooting (INSANE)
- "insane: command not found" oder Fehler beim Start:
  - Image neu bauen: `docker compose build insane`
  - Teste CLI direkt: `docker compose run --rm insane --help`
- Keine Dateien im Run-Ordner erzeugt:
  - Prüfe Logs: `${RUN_DIR}/logs/insane_membrane.log`
  - Stelle sicher, dass `${RUN_DIR}` existiert (siehe Schritt 1 der ColabFold-Anleitung)
- Globales Template schlägt fehl oder wird außerhalb des Repos geschrieben:
  - Setze `directories.templates` in `config/solvia_config.yaml` auf `"/work/data/templates"`
- Rechte-/Ownership-Probleme:
  - `.env` mit UID/GID anlegen (siehe ganz oben) oder korrigiere besitzende User: `sudo chown -R $(id -u):$(id -g) ${RUN_DIR}`

---

Nächster Schritt: GROMACS-Setup und -Simulation im Container (Equilibration, Production, Analyse).
\
\
# ⚙️ Schritt 4: GROMACS – Systemaufbau, Equilibration, Produktion (Container)

Ziel: Aus `system/*.gro/.top` ein laufendes MD-System erzeugen, dann Equilibration (EM→NVT→NPT) und eine Produktionssimulation starten.

## Einmalig vorbereiten (GROMACS)
```bash
# GPU-Image bauen (enthält gmx + Python libs)
docker compose build gromacs

# GPU testen (Host):
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
```

## 0) Peptide in die Membran einsetzen (wenn noch nicht erfolgt)
```bash
RUN_DIR="simulations/solvia_1_run_1"
docker compose run --rm --entrypoint python3 insane \
  scripts/universal/05_insert_peptides.py "${RUN_DIR}" --occupancy low

ls -l "${RUN_DIR}"/system
```

Einzelnes Peptid (Single-Peptide, zentriert) statt Occupancy:
```bash
RUN_DIR="simulations/solvia_1_run_1"
docker compose run --rm --entrypoint python3 insane \
  scripts/universal/05_insert_peptides.py "${RUN_DIR}" --n-peptides 1

# Ausgaben: system/system_n1.gro und system/system_n1.top
ls -l "${RUN_DIR}"/system
```

Mehrere Replikate (unabhängige Seeds/Platzierungen/Orientierungen):
```bash
RUN_DIR="simulations/solvia_1_run_1"
docker compose run --rm --entrypoint python3 insane \
  scripts/universal/05_insert_peptides.py "${RUN_DIR}" --n-peptides 1 --replicates 2

# Ausgaben: system/system_n1_rep1.* und system/system_n1_rep2.*
ls -l "${RUN_DIR}"/system
```

## 1) Equilibration (EM → NVT → NPT)
```bash
RUN_DIR="simulations/solvia_1_run_1"
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/06_equilibrate.py "${RUN_DIR}" --occupancy low

# Outputs prüfen
ls -l "${RUN_DIR}"/equilibration/em
ls -l "${RUN_DIR}"/equilibration/nvt
ls -l "${RUN_DIR}"/equilibration/npt
cat "${RUN_DIR}"/equilibration/summary.yaml | cat
```

## 2) Produktion starten
```bash
RUN_DIR="simulations/solvia_1_run_1"
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/07_run_production.py "${RUN_DIR}" --occupancy low --time 2

# Fortschritt/Logs überwachen
tail -f "${RUN_DIR}"/production/production.log
```
### Replikate sequentiell (kein Parallelbetrieb der GPU)
```bash
RUN_DIR="simulations/solvia_1_run_1"

# Equilibration Rep1
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/06_equilibrate.py "$RUN_DIR" --tag n1_rep1

# Produktion Rep1
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/07_run_production.py "$RUN_DIR" --tag n1_rep1 --time 2

# Equilibration Rep2 (erst wenn Rep1 fertig oder bewusst nacheinander absetzen)
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/06_equilibrate.py "$RUN_DIR" --tag n1_rep2

# Produktion Rep2
docker compose run --rm --entrypoint python3 gromacs \
  scripts/universal/07_run_production.py "$RUN_DIR" --tag n1_rep2 --time 2
```

Hinweis: Keine parallelen `docker compose run`-Aufrufe starten; die obigen Kommandos nacheinander absetzen, um GPU-Konkurrenz zu vermeiden.

Hinweise:
- `--time` in ns; bei >10 ns läuft der Prozess im Hintergrund und schreibt Log nach `${RUN_DIR}/production/production.log`.
- GPU wird automatisch genutzt (NVIDIA Container Toolkit erforderlich).

## 3) Batch-Modus (Equilibration + Produktion)
```bash
for run_dir in simulations/*_run_*; do
  echo "== Equilibrate: $run_dir =="
  docker compose run --rm --gpus all --entrypoint python gromacs \
    scripts/universal/06_equilibrate.py "$run_dir" --occupancy low || continue

  echo "== Production: $run_dir =="
  docker compose run --rm --gpus all --entrypoint python gromacs \
    scripts/universal/07_run_production.py "$run_dir" --occupancy low --time 2 || true
done
```

## Troubleshooting (GROMACS)
- "could not select device driver with capabilities: [[gpu]]": NVIDIA Container Toolkit installieren/konfigurieren; teste `docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi`.
- grompp-Fehler (`.mdp`, Topologie): Siehe Meldung in der Konsole; prüfe `system/*.top`, `equilibration/mdp/*.mdp`. Erhöhe notfalls `-maxwarn` im Skript minimal.
- Produktion startet nicht: Prüfe, ob `equilibration/summary.yaml` auf `npt_completed: true` steht und Dateien unter `equilibration/npt/` existieren.