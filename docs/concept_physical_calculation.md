# Physikalisches Kalibrierkonzept – Arbeitsstand Oktober 2023

## 1. Ausgangslage

- PMF-Auswertungen liefern aktuell pro Run:
  - ΔG-Profile (MBAR, optional Bootstrap),
  - effektive Partitionierungsintegrale `Kp_eff`,
  - Adsorptionskennwerte (ΔG_ads, θ@1µM, etc.),
  - neue CG-basierte Größen (Fußabdruck, R<sub>g,xy</sub>, Sequenzlänge).
- Die Kalibrierpipeline (`analysis/calibrate_hc50.py`) reichert diese Daten mit DBAASP-HC50 an und erzeugt drei Modelle:
  1. globaler Skalierungsfaktor (1-Parameter),
  2. Gamma-Modell `Γ* = exp(α·ΔG + β)`,
  3. Rang-Modell (Gewichtete Linearkombination, monotone Kalibrierung).
- Aktueller Rang-Modell-Fit (n = 10 Peptide) erreicht Spearman ρ ≈ 0.88 (RMSE_log10 ≈ 0.21), jedoch negative LOO-Spearman → Hinweis auf geringe Datenbasis und empfindliche Ausreißer.

## 2. Ziele

1. **Validität erhöhen**
   - Größeres Datenfundament (weitere Peptide, Replikate bereinigen).
   - Stabilere Modelle (LOOCV Spearman ≥ 0.6).
2. **Physikalische Konsistenz**
   - Fußabdruck-spezifische L/P*-Kalibrierung,
   - strukturierte Einbindung des Bilayer-Faktors,
   - Logging aller QC-Flags für Downstream-Gewichtungen.
3. **Automatisierung**
   - PMF-Analyse + Kalibrierung in einer reproduzierbaren Pipeline (Runner).
   - Möglichkeit, finale Modelle im PMF-Skript zu nutzen, sobald validiert.

## 3. Footprint-spezifisches L/P*

1. `scripts/analysis/pmf_mbar_analysis.py` liefert nun `footprint_cg_nm2` und `sequence_length_cg`.
2. Anpassung der L/P*-Berechnung:
   - aktuelles Modell: `L/P* = 1 / (A_lipid * Γ*)`, mit festem A<sub>lipid</sub> = 0.62 nm².
   - geplantes Konzept: `L/P* = f(footprint) = (Footprint / A_lipid) * packing_factor`.
   - Konfiguration (`config/pmf_standard_config.yaml`) sollte ein Mapping erlauben:
     ```yaml
     pmf:
       adsorption:
         lp_star_mode: auto_from_footprint
         footprint:
           contact_fraction: 0.4
           windows_max: 4
           pack_factor: 1.0
     ```
   - Pipeline speichert sowohl CG- als auch PMF-basierte Footprint-Schätzer; Kalibrierung entscheidet per QC/Regressionsprüfung, welches Maß besser korreliert.
3. In der Kalibrierung `analysis/calibrate_hc50.py`:
   - Neue Features `footprint_area_nm2`, `footprint_cg_nm2` in das Rang-Modell aufnehmen.
   - Ziel: Spearman/MSE für Footprint-gewichtete Varianten auswerten.

## 4. Bilayer-Factor

- Der Bilayer-Faktor (2.0) spiegelt momentan wider, dass die Membran von beiden Seiten zugänglich ist, obwohl pro Run nur einseitig sondiert wird.
- Validierungsstrategie:
  1. Kalibrierung mit Standardwert 2.0 beibehalten, aber in der Auswertung Vergleichsgröße `Kp_eff_single_leaflet = Kp_ads_nm` bereithalten.
  2. Sensitivität testen: Variation des Faktors (1.0 vs. 2.0) auf Rang-Korrelation prüfen.
  3. Option in `config/pmf_standard_config.yaml` dokumentieren (`adsorption.bilayer_factor`).
- Falls sich zeigt, dass die Rang-Korrelation mit Faktor 1.0 konsistenter ist (oder LOOCV stabiler wird), Umstellung diskutieren.

## 5. Weitere physikalische Correctives (Backlog)

1. **Elektrostatische Screening-Effekte**
   - Variation der Salzkonzentration (Parameter `membrane.salt_concentration`) kann in zukünftigen Simulationen für sensible Peptide berücksichtigt werden.
2. **Lipid-Asymmetrie / Cholesterol**
   - Der RBC-Mix ist fix; bei Peptiden, die stärker mit Phosphatidylserin interagieren, könnte eine Variation des Leaflet-Gehalts den Kp-Eff beeinflussen.
   - Idee: Alternativer `membrane` Config-Slot (`pmf.membrane_alternative`) für Sensitivitätsstudien.
3. **PMF-Normalisierung**
   - Bilayer-Factor vs. einzelseitiges Potential.
   - Prüfung, ob Bulk-Minima robust erkannt werden; ggf. `bulk_detection`-Parameter feiner abstimmen.
4. **Bootstrap/Unsicherheiten**
   - Bootstrap zwingend laufen lassen; in YAML/CSV die Konfidenzintervalle für `Kp_eff` und `Γ*` notieren → Rang-Modell kann Unsicherheiten als Feature nutzen.

## 6. Geplanter Ablauf (heute)

1. **PMF-Reanalyse**
   - Für alle Runs `python scripts/analysis/pmf_mbar_analysis.py <pmf_dir>` (mit Bootstrap).
   - Keine manuellen Löschungen nötig – Skript überschreibt YAML & Plots.
   - Sicherstellen: `pymbar` installiert.
2. **CSV-Export & Kalibrierung**
   - `python scripts/analysis/csv_export.py`
   - `python analysis/calibrate_hc50.py --input results_summary.csv --meta data/raw/peptides.csv`
   - Ergebnischeck: `analysis/calibration_summary.txt`, `analysis/results_calibrated.csv`.
3. **Ausreißerdiagnose**
   - Peptide mit großen `HC50_rank_resid_log10` analysieren (PMF-Plot, QC-Werte).
   - Läufe mit schlechten `QC_pass_fraction`, `ESS` etc. priorisiert neu rechnen oder exkludieren.
4. **Footprint/L/P*-Kalibrierung**
   - Zusätzliche Regressionen in `calibrate_hc50.py`: `Γ* = g(ΔG, footprint)`.
   - Vergleich `Lp_auto` vs. fix -> Spearman/LOO bewerten.

## 7. Nächste Iteration (spätestens morgen/kommende Runde)

- Wenn die Rangkalibrierung stabil ist:
  - Modellparameter einfrieren (Gewichte + monotone Kurve in JSON/Config ablegen).
  - PMF-Skript so erweitern, dass es `HC50_rank_uM` auf Basis dieser Parameter direkt berechnet.
- Dokumentation updaten (README_HC50, usage.md).
- Langfristig: zusätzliche physikalische Variablen (z. B. Kontaktflächen, Leitmotiv-Sequence, experimentelle Zelllinie) einführen und testen.

---
Stand: 2025-09-22 21:45 UTC
