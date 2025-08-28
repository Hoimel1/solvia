# 🧬 SOLVIA Pipeline - PMF-optimierte Toxizitätsvorhersage

## Übersicht

SOLVIA ist eine **vollständig containerisierte** Pipeline zur Vorhersage der hämolytischen Toxizität antimikrobieller Peptide mittels PMF-Berechnungen (Potential of Mean Force). Die Pipeline nutzt Coarse-Grained Molecular Dynamics mit dem MARTINI 3 Kraftfeld.

**🚀 Wichtig: Alle Tools laufen über Docker Container** - keine lokale Installation erforderlich!
- ColabFold, Martinize, INSANE und GROMACS werden automatisch via Docker aufgerufen
- GPU-Support für GROMACS ist integriert (NVIDIA Runtime)

## Pipeline-Architektur

```
ColabFold → Martinize2 → INSANE → GROMACS → PMF/Umbrella → MBAR → Toxizitätsklassifikation
```

## Voraussetzungen

- Docker und Docker Compose (für ColabFold und Martinize2)
- GROMACS 2023.3 mit GPU-Unterstützung
- Python 3.8+ mit pymbar
- NVIDIA GPU (L4 oder besser)
- ~50 GB Speicher pro Peptid

## Schritt 1: Setup und Strukturvorhersage mit ColabFold

### 1.1 Run-Verzeichnis anlegen

```bash
# Arbeitsverzeichnis
cd /home/michelhuller/solvia

# Run-Setup für einzelnes Peptid
python3 scripts/universal/01_setup_run.py data/raw/fasta/SOLVIA_1.fasta

# Struktur:
# simulations/solvia_1_run_1/
#   ├── input/peptide.fasta
#   └── metadata.yaml
```

### 1.2 ColabFold im Container ausführen

```bash
# Variable für Run-Verzeichnis
RUN_DIR="simulations/solvia_1_run_1"

# ColabFold für Strukturvorhersage
docker compose run --rm \
  -e XDG_CACHE_HOME=/cache \
  -e MPLCONFIGDIR=/cache/mpl \
  colabfold \
  /work/${RUN_DIR}/input/peptide.fasta /work/${RUN_DIR}/colabfold \
  --num-seeds 3 --num-models 5 \
  --msa-mode mmseqs2_uniref_env \
  --pair-mode unpaired_paired

# Bestes Modell auswählen
python3 scripts/universal/02_run_colabfold.py ${RUN_DIR}
```

**Output:** 5 Modelle mit pLDDT-Scores in `${RUN_DIR}/colabfold/`

## Schritt 2: Coarse-Graining mit Martinize2

### 2.1 Automatisiertes Coarse-Graining

```bash
# EMPFOHLEN: Python-Script (vollautomatisch)
python3 scripts/universal/03_coarse_grain.py ${RUN_DIR}

# Das Script handhabt automatisch:
# - Docker-Container-Aufruf von Martinize
# - ITP-Datei-Verwaltung (verhindert Doppelzählung)
# - Korrekte Molekülbenennung (SOLVIA_1, nicht SIMULATIONS)
# - Symlink-Erstellung für Kompatibilität
# - Position Restraint Generation
# - C-terminale Amidierung (NH2)
```

**Parameter-Erklärung:**
- `-ff martini3001`: MARTINI 3 Kraftfeld
- `-name ${PEPTIDE_ID}`: Molekülname (z.B. SOLVIA_1)
- `-cter NH2-ter`: **C-terminale Amidierung** (NH2 statt COO⁻)
- `-dssp`: Sekundärstruktur aus DSSP
- `-elastic`: Elastic Network für Stabilität
- `-ef 700`: Federkonstante (kJ/mol/nm²)
- `-p backbone`: Position Restraints auf Backbone
- `-pf 1000`: Restraint-Kraft (kJ/mol/nm²)

## Schritt 3: Membranaufbau mit INSANE

### 3.1 Standardisierte RBC-Membran

```bash
# Globale Membran-Template erstellen (einmalig)
python3 scripts/universal/04_build_membrane.py --global-template

# Oder: Membran für spezifischen Run
python3 scripts/universal/04_build_membrane.py ${RUN_DIR}

# Erzeugt:
# data/templates/membrane/  (bei --global-template)
# oder ${RUN_DIR}/membrane/
#   ├── membrane.gro           # Struktur  
#   ├── membrane.top           # Topologie
#   ├── index_leaflets.ndx     # Leaflet-Definitionen
#   └── posre_lipid_z.itp      # Z-Restraints
```

### 3.2 Alternative: Custom-Membran mit INSANE

```bash
# Für spezielle Membranzusammensetzungen
insane \
  -o custom_membrane.gro \
  -p custom_membrane.top \
  -box 10,10,14 \
  -u POPC:0.40 -u PSM:0.15 -u CHOL:0.45 \
  -l POPE:0.45 -l POPS:0.15 -l CHOL:0.40 \
  -sol W \
  -salt 0.15 \
  -seed 42
```

**Leaflet-Komposition:**
- **Outer (RBC-typisch):** 40% POPC, 15% SM, 45% Cholesterol
- **Inner:** 45% POPE, 15% POPS, 40% Cholesterol

## Schritt 4: Peptid-Insertion für PMF

### 4.1 PMF-optimierte Insertion mit Replikaten

```bash
# Single Peptide für PMF (Standard)
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation parallel

# Alternative: Mit verschiedenen Orientierungen für Replikate
# Replikat 1: Parallel
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation parallel

# Replikat 2: 45° Neigung
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation tilt45
```

**Features:**
- Single Peptide für PMF-Studien
- Kontrollierte Orientierung (parallel, tilt45, perpendicular)
- Flexible Platzierung für verschiedene Szenarien

### 4.2 Alternative: Multi-Peptid-System (für Aggregationsstudien)

```bash
# 8 Peptide in 3-2-3 Anordnung
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --membrane membrane_standard/membrane.gro \
  --n-peptides 8 \
  --arrangement "3-2-3" \
  --distance 3.0
```

## Schritt 5: Equilibrierung mit GROMACS

### 5.1 PMF-spezifische Equilibrierung

```bash
# Standard Equilibrierung für PMF
python3 scripts/universal/06_equilibrate.py ${RUN_DIR}

# Mit custom Tag für verschiedene Runs
python3 scripts/universal/06_equilibrate.py ${RUN_DIR} \
  --tag pmf_rep1

# Protokoll (6.5 ns total):
# 1. EM mit Restraints (50k steps)
# 2. NVT (0.5 ns, 100K → 310K)
# 3. NPT-pre (1 ns, mit Restraints)
# 4. NPT-prod (5 ns, nur Lipid-Z-Restraints)
```

### 5.2 Manuelle MDP-Kontrolle

```bash
# Energy Minimization
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr
gmx mdrun -v -deffnm em

# NVT Equilibration
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro
gmx mdrun -v -deffnm nvt

# NPT Equilibration
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro
gmx mdrun -v -deffnm npt
```

## Schritt 6: PMF-Berechnung mit Umbrella Sampling

### 6.1 PMF mit lokaler Patch-Referenz (EMPFOHLEN)

```bash
# Enhanced PMF mit lokaler Patch-Referenz
python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
  --replicate 1 \
  --tag pmf_local

# Für mehrere Replikate:
for REP in 1 2 3; do
  python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
    --replicate ${REP} \
    --tag pmf_rep${REP}
done

# Features:
# - Dynamische Patch-Referenz (2 nm Radius)
# - Adaptive Fenster-Verdichtung
# - SMD-Initialisierung
# - Automatische QC-Gates
```

**Window-Strategie:**
- **Bulk:** +2.8, +2.4 nm
- **Coarse:** 0.2 nm Schritte
- **Dense (Interface):** +0.6 bis -0.6 nm, 0.15 nm Schritte
- **Deep:** bis -2.0 nm

## Schritt 7: PMF-Analyse mit MBAR/WHAM

### 7.1 MBAR-Analyse mit Bootstrap

```bash
# Analyse für jedes Replikat
python3 scripts/analysis/pmf_mbar_analysis.py \
  ${RUN_DIR}/pmf/replicate_1 \
  --method mbar \
  --bootstrap 1000

# Output:
# - pmf_analysis_results.yaml    # Features & Metriken
# - analysis_plots/
#   ├── pmf_profile.png          # PMF mit 95% CI
#   ├── overlap_matrix.png       # Window-Overlap
#   └── convergence.png          # Konvergenz-Check
```

### 7.2 Feature-Extraktion

```yaml
# Extrahierte Features (pmf_analysis_results.yaml):
features:
  delta_g_ads: -12.3        # Adsorptionsenergie (kJ/mol)
  delta_g_insert: -5.7      # Insertionsenergie (kJ/mol)
  delta_g_barrier: 8.2      # Aktivierungsbarriere (kJ/mol)
  z_surf_min: 0.35          # Oberflächenminimum (nm)
  z_head_min: -0.82         # Kopfgruppen-Minimum (nm)
```

## Schritt 8: Qualitätskontrolle

### 8.1 Umfassendes QC-System

```bash
# QC-Report mit Korrekturvorschlägen
python3 scripts/analysis/pmf_qc_system.py \
  ${RUN_DIR}/pmf/replicate_1 \
  --suggest

# Prüft:
# - Window-Overlap (min 10%, target 20%)
# - Konvergenz (Halbzeit-Differenz < 2 kJ/mol)
# - ESS (≥ 200 unabhängige Frames)
# - Replikat-Konsistenz (< 2 kJ/mol Differenz)
```

### 8.2 Manuelle Checks

```bash
# Overlap-Matrix visualisieren
gmx analyze -f ${RUN_DIR}/pmf/replicate_1/overlap.xvg -dist

# Konvergenz prüfen
for window in ${RUN_DIR}/pmf/replicate_1/windows/z_*/; do
  echo "Window: $(basename $window)"
  tail -n 100 $window/pullx.xvg | awk '{sum+=$2; n++} END {print sum/n}'
done
```

## Schritt 9: Toxizitätsklassifikation

### 9.1 Threshold-basierte Klassifikation

```python
# Toxizitätskriterien
toxic = (delta_g_ads < -8.0) and \
        (delta_g_insert <= -3.0 or delta_g_barrier < 12.0)
```

### 9.2 Replicate Averaging

```bash
# Mitteln über Replikate (präzisionsgewichtet)
python3 scripts/analysis/average_replicates.py \
  ${RUN_DIR}/pmf/replicate_*/pmf_analysis_results.yaml \
  --output ${RUN_DIR}/final_results/averaged_features.yaml
```

## Automatisierte Gesamt-Pipeline

### Master-Script für kompletten Workflow

```bash
# Alles in einem Befehl!
./run_pmf_pipeline.sh 1 --replicates 2

# Batch-Processing für mehrere Peptide
for ID in 1 8 12 215; do
  ./run_pmf_pipeline.sh ${ID} --skip-colabfold &
done
wait
```

### Wichtige Dateibenennung

Nach dem Coarse-Graining existieren folgende Dateien:
```
coarse_grain/
├── SOLVIA_1_cg.pdb      # Coarse-grained Struktur  
├── SOLVIA_1.itp         # Molekül-Topologie (nur diese eine ITP!)
├── peptide_cg.pdb       # Symlink → SOLVIA_1_cg.pdb
├── peptide_cg.gro       # Symlink → SOLVIA_1_cg.pdb  
├── peptide_cg.itp       # Symlink → SOLVIA_1.itp
└── peptide_posre.itp    # Position restraints
```

**Wichtig:** 
- Nur EINE ITP-Datei verwenden (sonst Doppelzählung!)
- Martinize2 erstellt initial zwei ITPs (_0.itp und .itp), wir verwenden nur die _0.itp (umbenannt zu .itp)
- Molekülname in der ITP: `SOLVIA_1` (ohne _0)
- Die Symlinks mit "peptide_*" Namen sind für Kompatibilität mit anderen Scripts

## Erwartete Performance-Metriken

### Validierungsziele
- **Spearman ρ ≥ 0.5** für HC50-Korrelation
- **ROC-AUC ≥ 0.80** für binäre Klassifikation
- **95% CI** durch Bootstrap (n=1000)

### Computational Resources
- **Pro Peptid:** 0.64-0.96 μs MD
- **GPU-Zeit:** 24-48h auf NVIDIA L4
- **Storage:** ~50 GB pro Peptid

## Troubleshooting

### Problem: Niedrige Window-Überlappung
```bash
# Automatisch gelöst durch adaptive Verdichtung
# Oder manuell Fenster hinzufügen:
python3 scripts/universal/add_umbrella_window.py \
  --center -0.25 --k 900 --time 20
```

### Problem: Schlechte Konvergenz
```bash
# Simulation verlängern
python3 scripts/universal/extend_window.py \
  --window z_+0.000 --add-ns 10
```

### Problem: Replikat-Divergenz
```bash
# Vergleich der PMF-Profile
python3 scripts/analysis/compare_replicates.py ${RUN_DIR}/pmf
```

## Konfigurationsdateien

### Zentrale PMF-Konfiguration
```yaml
# config/pmf_standard_config.yaml
pmf:
  umbrella:
    ref_mode: "patch"      # Lokale Referenz!
    patch_radius: 2.0      # nm
    force_constant: 900    # kJ/mol/nm²
    production_ns: 20.0
  qc:
    min_neighbor_overlap: 0.10
    target_overlap: 0.20
```

## Wissenschaftliche Referenzen

- MARTINI 3: Souza et al., Nature Methods 2021
- PMF-Konzept: `docs/pmf_concept.md`
- Kritische Verbesserungen: `docs/aenderungen.md`
- Validierungsdaten: `data/pmf_validation/`

## Support

Bei Fragen zur Pipeline:
1. QC-Reports prüfen: `pmf/replicate_*/qc_full_report.yaml`
2. Logs checken: `grep -i error logs/*.log`
3. Dokumentation: `docs/usage_pmf.md` (technische Details)