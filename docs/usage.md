# 🧬 SOLVIA Pipeline - PMF-optimierte Toxizitätsvorhersage

## Übersicht

SOLVIA ist eine **vollständig containerisierte** Pipeline zur Vorhersage der hämolytischen Toxizität antimikrobieller Peptide mittels PMF-Berechnungen (Potential of Mean Force). Die Pipeline nutzt Coarse-Grained Molecular Dynamics mit dem MARTINI 3 Kraftfeld.

**🚀 Wichtig: Alle Tools laufen über Docker Container** - keine lokale Installation erforderlich!
- ColabFold, Martinize, INSANE und GROMACS werden automatisch via Docker aufgerufen
- GPU-Support für GROMACS ist integriert (NVIDIA Runtime)

## Pipeline-Architektur

```
ColabFold → Martinize2 → INSANE → GROMACS → PMF/Umbrella → MBAR (Standard) → Toxizitätsklassifikation
```

## Reproduzierbare Standard-Pipeline (empfohlen)

```bash
# 0) Optional: RBC-Membran-Template einmalig erstellen (validierte Asymmetrie)
python3 scripts/universal/04_build_membrane.py --global-template

# 1) Run-Verzeichnis anlegen (Beispiel: SOLVIA_1)
python3 scripts/universal/01_setup_run.py data/raw/fasta/SOLVIA_1.fasta

# 2) ColabFold (siehe Abschnitt 1.2), danach bestes Modell wählen
python3 scripts/universal/02_run_colabfold.py simulations/solvia_8_run_1

# 3) Coarse-Graining (ITP-Aufbereitung inkl. PosRes)
python3 scripts/universal/03_coarse_grain.py simulations/solvia_8_run_1

# 4) Peptid-Insertion (2 Replikate, evidenzbasiert)
# Rep 1: parallel; Rep 2: kontinuierliche Neigung + Roll
python3 scripts/universal/05_insert_peptides.py simulations/solvia_1_run_1 \
  --n-peptides 1 --orientation parallel
python3 scripts/universal/05_insert_peptides.py simulations/solvia_1_run_1 \
  --n-peptides 1 --orientation continuous

# 5) Equilibrierung (Standard)
python3 scripts/universal/06_equilibrate.py simulations/solvia_1_run_1 --tag pmf_rep1

# 6) PMF/Umbrella pro Replikat (60 ns Fenster, QC/Auto-Densify gemäß Config)
python3 scripts/universal/08_run_pmf.py simulations/solvia_1_run_1 --replicate 1
python3 scripts/universal/08_run_pmf.py simulations/solvia_1_run_1 --replicate 2

# 7) Analyse (MBAR Standard; reproduzierbar ohne Bootstrap; Replikate aggregieren)
# Ein Replikat/Tag analysieren
python3 scripts/analysis/pmf_mbar_analysis.py \
  simulations/solvia_1_run_1/pmf/pmf_midplane \
  --no-bootstrap

# Alternativ: Alle Replikate unter ${RUN_DIR}/pmf analysieren und aggregieren
python3 scripts/analysis/pmf_mbar_analysis.py \
  simulations/solvia_1_run_1/pmf \
  --aggregate \
  --bootstrap 0

# 8) HC50-Migration (2025-Felder in pmf_analysis_results.yaml ergänzen)
# Scannt rekursiv unter ${RUN_DIR}/pmf und schreibt fehlende adsorption.kp/hc50-Felder
# Mit RUN_DIR-Variable (aus obigen Schritten):
python3 scripts/migrate_2025_hc50.py \
  ${RUN_DIR}/pmf \
  --write

Hinweis: Der Aufruf sollte aus dem Repository‑Root erfolgen (damit das Modul
`analysis.hc50` importiert werden kann). Das Migrationsskript setzt den Repo‑Root
zwar automatisch auf den Python‑Pfad, aber am sichersten ist:

# im Repo‑Root ausführen
cd /home/michelhuller/solvia
python3 scripts/migrate_2025_hc50.py ${RUN_DIR}/pmf --write

### Optional: Regionale QC-Schwellen und Auto-Stride

In `config/pmf_standard_config.yaml` können regionenabhängige QC-Grenzen und ein automatischer ESS‑Stride aktiviert werden. Beispiel (Drop‑in, nur relevante Keys):

```
pmf:
  qc:
    # global defaults (werden von region_thresholds überschrieben)
    min_neighbor_overlap: 0.20
    min_ess_frames: 200
    half_energy_tol_kj: 2.0
    half_z_tol_sigma: 2.0
    discard_fraction: 0.1
    ess_stride: "auto"       # automatisch: ~0.2 * tau_int als Abtastabstand
    region_thresholds:
      bulk:
        min_neighbor_overlap: 0.12
        min_ess_frames: 150
        half_energy_tol_kj: 3.0
        half_z_tol_sigma: 2.5
        min_time_ns: 20
      approach:
        min_neighbor_overlap: 0.15
        min_ess_frames: 200
        half_energy_tol_kj: 2.0
        half_z_tol_sigma: 2.0
        min_time_ns: 30
      interface:
        min_neighbor_overlap: 0.20
        min_ess_frames: 300
        half_energy_tol_kj: 1.0
        half_z_tol_sigma: 1.5
        min_time_ns: 40
```
```

Hinweise:
- Membran: validierte RBC-Asymmetrie mit ~42% CHOL ist per Default aktiv; Leaflets werden per Index (`index_leaflets.ndx`) definiert.
- PMF: 60 ns/Fenster; QC/Auto‑Densify/Extend gemäß `config/pmf_standard_config.yaml`.
- Analyse: MBAR ist Standard; deterministisch mit `--no-bootstrap`. Bootstrap optional (z. B. `--bootstrap 100`).
- Reproduzierbarkeit: Seeds für Replikate sind deterministisch; Bootstrap ist per Default deaktiviert.

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
python3 scripts/universal/01_setup_run.py data/raw/fasta/SOLVIA_12.fasta

# Struktur:
# simulations/solvia_1_run_1/
#   ├── input/peptide.fasta
#   └── metadata.yaml
```

### 1.2 ColabFold im Container ausführen

```bash
# Variable für Run-Verzeichnis
RUN_DIR="simulations/solvia_1219_run_1"

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

Evidenzbasiert empfehlen wir eine Kombination aus paralleler und leicht geneigter Startlage mit azimutalem “Rolling”, um Orientierungsbias zu reduzieren und robuste z‑PMFs zu erhalten.

Option A — Ein Replikat (robuster Einzel‑Run)

```bash
# Kontinuierliche Neigung (Tilt ~ N(30°, 20°) in [0°,90°]),
# Roll (Azimut) ~ U(0,360°):
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation continuous
```

Option B — Zwei Replikate (empfohlen)

```bash
# Replikat 1: Parallel (0° Tilt), roll wird in der Dynamik variiert
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation parallel

# Replikat 2: Kontinuierliche Neigung + Roll (siehe oben)
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 1 \
  --orientation continuous
```

Mehrere Peptide (Poisson‑Sampling mit kooperativer Entzerrung)

```bash
# 8 Peptide mit Poisson‑Disk‑Sampling (XY) und MC‑Entzerrung (“cooperative”)
python3 scripts/universal/05_insert_peptides.py ${RUN_DIR} \
  --n-peptides 8 \
  --placement cooperative \
  --orientation continuous
```

Hinweise:

- `--orientation continuous` sampelt Tilt kontinuierlich (um 30°) und Roll uniform (0–360°).
- Für einzelne Peptide reicht `--placement poisson` (Default). Bei mehreren Peptiden verbessert `--placement cooperative` die Startabstände (keine Überschneidungen, realistischere Nachbarschaften).
- Seeds sind deterministisch pro Replikat (über Run‑Infos), Replikate unterscheiden sich in Orientierung/Platzierung.

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

### 6.1 PMF mit lokaler Midplane-Referenz (EMPFOHLEN)

```bash
# Enhanced PMF mit lokaler Midplane-Referenz
python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
  --replicate 1 \
  --tag pmf_midplane \
  --resume

# Für mehrere Replikate:
for REP in 1 2 3; do
  python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
    --replicate ${REP} \
    --tag pmf_rep${REP}
done

# Strenger QC: bricht ab, falls ESS/Overlap-Gates reißen
python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
  --replicate 1 \
  --tag pmf_midplane \
  --strict-qc


> Hinweis: Dauerhaft aktivierbar via `pmf.qc.strict_mode: true` in `config/pmf_standard_config.yaml`.

# Features:
# - Lokale Midplane-Referenz (balancierte Outer/Inner PO4, ~2.5 nm Radius; pbcatom=COM)
# - Adaptive Fenster-Verdichtung
# - SMD-Initialisierung
# - Optionale kurze Pre‑Equilibration am Ziel‑z (pre_equil_ns)
# - Automatische QC-Gates
# - Logging (run.log), Resume & QC-only
```

**Window-Strategie:**
- **Bulk:** +2.8, +2.4 nm
- **Coarse:** 0.2 nm Schritte
- **Dense (Interface):** +0.6 bis -0.6 nm, 0.10 nm Schritte
- **Deep:** bis -2.0 nm

Hinweis zur Referenz:
- Standard ist `ref_mode: midplane_local` (balanciert Outer/Inner PO4 lokal um das Peptid; robust gegen Asymmetrien). Ohne Leaflet-Index (`index_leaflets.ndx`) wird die Trennung fallback-basiert aus der z-Verteilung abgeleitet.
- Produktionsläufe nutzen keine Lipid-Z-Restraints (physikalisch realistischer Reorganisationsraum).

### 6.2 Resume, QC-only und Logging

```bash
# Resume: vorhandene Fenster werden übersprungen, Densify/Extend & QC laufen weiter
python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
  --replicate 1 \
  --tag pmf_midplane \
  --resume

# QC-only: keine weitere MD — erzeugt nur Reports/Metadata aus vorhandenen Fenstern
python3 scripts/universal/08_run_pmf.py ${RUN_DIR} \
  --replicate 1 \
  --tag pmf_midplane \
  --qc-only

# Outputs (pro Tag/Replicate):
# ${RUN_DIR}/pmf/<tag>/run.log         # Laufprotokoll
# ${RUN_DIR}/pmf/<tag>/qc_report.yaml  # QC-Zusammenfassung
# ${RUN_DIR}/pmf/<tag>/pmf_metadata.yaml  # Metadaten für Analyse (inkl. windows[].k)
# ${RUN_DIR}/pmf/<tag>/RUN_INFO.yaml   # Provenienz (Git, Container, Env)
# ${RUN_DIR}/pmf/<tag>/RESIMULATE.flag # Marker bei aktivem strict/qc_audit
```

## Schritt 7: PMF-Analyse mit MBAR

### 7.1 MBAR (Standard)

```bash
# Analyse für jedes Replikat
python3 scripts/analysis/pmf_mbar_analysis.py \
  ${RUN_DIR}/pmf/pmf_midplane \
  --no-bootstrap 

python3 scripts/analysis/pmf_mbar_analysis.py \
  ${RUN_DIR}/pmf/pmf_midplane \
  --bootstrap 200

# Output:
# - pmf_analysis_results.yaml    # Features & Metriken
# - analysis_plots/
#   ├── pmf_profile.png          # PMF mit 95% CI
#   ├── overlap_matrix.png       # Window-Overlap
#   └── convergence.png          # Konvergenz-Check
# Zusatz (vergleichbar zwischen Peptiden):
# - adsorption.peptide_area_nm2 (falls konfiguriert/ableitbar)
# - adsorption.theta_at_1uM, adsorption.theta_at_hc50
# - optional: adsorption.hc50_uM_at_theta_star (bei konfiguriertem theta_star)
```

```bash
# Alternativ: Alle Replikate unter ${RUN_DIR}/pmf analysieren und aggregieren
python3 scripts/analysis/pmf_mbar_analysis.py \
  ${RUN_DIR}/pmf \
  --aggregate \
  --bootstrap 0

# Aggregate-Ausgabe:
# - pmf/pmf_analysis_aggregate.yaml  # gemittelte Features + Replicate-Consistency
# Pro Replikat weiter wie oben: pmf_analysis_results.yaml + analysis_plots/
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
# Automatisierter QC-Audit (gesamte simulations/-Hierarchy)
python3 scripts/analysis/qc_audit.py simulations \
  --output-dir analysis/qc \
  --mark-resim

# Outputs:
# - analysis/qc/qc_report.csv      # tabellarische Übersicht
# - analysis/qc/qc_report.md       # Markdown-Review
# - analysis/qc/qc_resimulate.yaml # Runs + Gründe für Re-Simulation
# - RESIMULATE.flag in pmf/<tag>/  # Marker für erneutes Sampling
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
# config/pmf_standard_config.yaml (Auszug)

Hinweis: HC50 wird ausschließlich aus der Adsorptions‑Thermodynamik (Wasserseite) abgeleitet. ΔG_insert ist optional/diagnostisch und fließt nicht in die Klassifikation ein. Details siehe docs/README_HC50.md.

### CSV‑Export (HC50/Kp Übersicht)

```bash
python scripts/analysis/export_results_csv.py simulations -o out/results.csv
```
Spalten: `peptide,replicate,kp_nm,kp_eff_nm,hc50_low_uM,hc50_high_uM,delta_g_ads,z_ads,label,qc_pass,path`.
pmf:
  umbrella:
    ref_mode: "hybrid"            # LocalMidplane → LocalPatch → UpperPO4 (Fallback)
    consistent_reference: true     # ein gemeinsamer Index für alle Fenster
    patch_radius: 2.5              # nm
    pre_smd_ns: 1.0                # kurzer SMD‑Nudge
    pre_equil_ns: 1.0              # kurzer Vorlauf am Ziel‑z (optional)
    force_constant: 900            # kJ/mol/nm²
    production_ns: 60.0
    lipid_z_posres: false
  qc:
    min_neighbor_overlap: 0.15
    target_overlap: 0.20
    min_ess_frames: 100
    ess_stride: "auto"
    max_extend_ns: 80.0
    region_thresholds:
      interface:
        min_ess_frames: 300
  analysis:
    method: "mbar"                # MBAR ist Standard
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
