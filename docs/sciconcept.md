Alles klar—hier ist das erneuerte, schlanke Validierungs- und Screening-Konzept mit 10×10×14 nm Patch, fokussiert auf ΔG_ads und ΔG_insert aus einer 1D-PMF pro Peptid.

Ziel & Hypothesen
	•	H1 (Ranking): -\Delta G_\text{insert} und -\Delta G_\text{ads} korrelieren mit \log_{10}(\mathrm{HC50}_\text{RBC}) (Spearman ρ ≥ 0.5).
	•	H1 (Klassifikation): ein lineares Modell aus ΔG-Features trennt „toxisch“ (HC50 ≤ 50 µM) vs. „nicht-toxisch“ (HC50 ≥ 250 µM) mit ROC-AUC ≥ 0.80.

Datensatz
	•	n = 40 Peptide: 20 toxische, 20 nicht-toxische (Borderline 50–250 µM ausschließen).
	•	Matching (grobe Balance zwischen Gruppen): Länge, Nettoladung @pH 7, hydrophober Moment.
	•	Chemie fix: N-Terminus protoniert (+1), C-Terminus amidiert (0), pH 7, 150 mM NaCl.

Simulationssystem (für alle identisch)
	•	Box: 10×10×14 nm (Z-Puffer ≥ 4 nm Wasser je Seite).
	•	Membran: RBC-außenähnlich, asymmetrisch; Cholesterin: 45 mol % (Primärsystem; Sensitivity-Subset bei 35 % zur Rang-Stabilität).
	•	APL-Faustformel: ≈ 0.50 nm² ⇒ ~200 Lipide/Leaflet (~400 gesamt).
	•	CG-Details (MARTINI 3): T=310 K, P=1 bar; Coulomb/LJ-Cutoff 1.1 nm; Membran-COM in z schwach restrain, damit die Referenz stabil bleibt.

Monomer-PMF (Adsorption/Insertion entlang z)

CV: Abstand in z zwischen Peptid-COM und äußerer Phosphat-Ebene.
Umbrellas: 16 Fenster, Δz≈0.20 nm, Bereich z = +2.8 … −1.4 nm; Feder k=800–1000 kJ mol⁻¹ nm⁻².
Laufzeit: 20 ns/Fenster.
Replikate: 2 Reps (unabhängige Startorientierungen); bei QC-Warnung gezielt Fenster auf 30 ns verlängern (nicht alles).
Starts: Peptid im Bulk (z≈+1.5 nm), XY-Platzierung via Poisson-Disk (min 1.5–2.0 nm), drei Orientierungspools (parallel/45°/senkrecht).
Analyse: MBAR (oder WHAM) → G(z) mit Bulk-Plateau =0.

Feature-Definitionen (aus G(z))
	•	z_\text{surf,min}: erstes Minimum nahe Kopfgruppenaußenkante.
	•	z_\text{head,min}: lokales Minimum tiefer im Kopfgruppenbereich.
	•	ΔG_ads = G(z_\text{surf,min}) - G_\text{bulk}.
	•	ΔG_insert = G(z_\text{head,min}) - G(z_\text{surf,min}) (negativ = insertion-freundlich).
	•	Barriere \Delta G^\ddagger = G(z^\ddagger) - G(z_\text{surf,min}).
	•	Unsicherheit: Block-Bootstrap pro Fenster, 1000 Resamples.
	•	Pooling über Reps: gewichtetes Mittel der Features mit Gewichten w_i = 1/\mathrm{SE}_i^2.

QC-Kriterien (Gate)
	•	Nachbarschafts-Overlap O_{i,i+1} > 0.10.
	•	≥ 200 unkorrelierte Proben/Fenster (aus Autokorrelation).
	•	Replikat-Konsistenz: |Feature_rep1 − Feature_rep2| ≤ 2 kJ/mol; sonst nur betroffene Fenster verlängern/verdichten.

Statistik & Validierung
	1.	Primär (Regression): Spearman-ρ von -\Delta G_\text{insert} (und separat -\Delta G_\text{ads}, -\Delta G^\ddagger) gegen \log_{10}(\mathrm{HC50}); 95 %-CI via Bootstrap.
	2.	Primär (Klassifikation): logistische Regression y\sim \Delta G_\text{ads} + \Delta G_\text{insert} + \Delta G^\ddagger, nested 5× CV. Kennzahlen: ROC-AUC (DeLong-CI), PR-AUC, Brier-Score, Kalibration.
	3.	Kontrolle von Kovariaten: partielle Rangkorrelation (Länge, Nettoladung, μ_H).
	4.	Vorab-Schwellen (Heuristik-Check):
	•	„RBC-kritisch“, wenn ΔG_ads < −8 kJ/mol und (ΔG_insert ≤ −3 kJ/mol oder ΔG‡ < 12 kJ/mol).
	•	Datengetriebene Cutoffs: Youden-Index je Trainings-Fold; auf Test-Fold anwenden; Median-Cutoff reporten.

Sensitivitäten (Robustheit)
	•	CHOL 45 % vs. 35 % auf einem Subset (10 Peptide) → Kendall-τ der Rankings.
	•	Elastic-Network an/aus auf 5 Peptiden → prüfen, ob Rangordnung stabil bleibt.
	•	Replikat-Robustheit: leave-one-rep-out der ΔG-Features.

Optional: extrem billiges Zusatz-Feature

Monomer-Thinning-Assay: 1 Peptid auf Oberfläche, 5×200 ns; Metriken Δ-Dicke (lokal), Midplane-Wasser, Tail-Order-Drop; Aggregation per Median.
→ Wird nur in die Statistik aufgenommen, wenn AUC/ρ aus PMF-Features < vordefiniertem Ziel liegt.

Rechenbudget (pro Peptid)
	•	PMF: 16×20 ns × 2 Reps = 0.64 µs (typisch; nur einzelne Fenster ggf. auf 30 ns).
	•	Optional Thinning: 5×200 ns = 1.0 µs.

Deliverables
	•	Für jedes Peptid: PMF-Plot mit 95 %-Band; Tabelle ΔG_ads, ΔG_insert, ΔG‡ ± SE.
	•	Gesamt: Korrelationplots (ΔG vs. \log_{10}HC50), ROC/PR-Kurven, Kalibrations-Plot, finale Cutoffs.

Implementierungs-Checkliste (praktisch)
	•	Platzierung: z-Target = z_{P,\text{outer}} + 0.6–1.2 nm; 50 ps NVT mit Flat-Bottom-Z-Restraint, dann frei.
	•	Seeds/Orientierungen: deterministische RNG-Seeds aus (peptide_id, rep).
	•	Termini: N-term +1, C-term amidiert 0 (im Mapping/itp prüfen).
	•	Energygrps: peptide vs. lipids für spätere QC.
	•	Auswertungspipeline: MBAR → Spline-Minima/Barriere suchen → Feature-CSV (±SE) → Statistik-Notebook (Spearman, logReg + CV, DeLong).

⸻

Interpretation & Entscheidungslogik
	•	Erfolg: ρ ≥ 0.5 und AUC ≥ 0.80 → PMF-Surrogat validiert; künftiges Tox-Screening = nur PMF (0.64 µs/Peptid).
	•	Grenzbereich: ρ 0.3–0.5 oder AUC 0.70–0.80 → Thinning-Assay als drittes Feature ergänzen.
	•	Schwach: ρ < 0.3 und AUC ≤ 0.65 → auf Duo „Seed-Pore-Assay + Thinning“ umschwenken.

Wenn du möchtest, formatiere ich dir daraus eine prä-registrierbare Markdown-Checkliste (inkl. feste QC-Schwellen) und eine kleine Config-Vorlage für 10×10×14 nm, damit du sofort loslegen kannst.

Hier ist die prä-registrierbare Markdown-Checkliste (inkl. fester QC-Schwellen) für das ΔG-basierte RBC-Tox-Screening mit 10×10×14 nm-Patch. Du kannst sie 1:1 in euer Repo legen (docs/prereg_pmf_rbc.md). Alles Relevante ist fixiert; nur die Peptidliste/HC50-Werte tragt ihr ein.

⸻

Prä-Registrations-Protokoll: ΔG_ads/ΔG_insert als Surrogat für RBC-Toxizität

Version: v1.0 (Datum: ___)
Kontakt/PI: ___
Repo-Commit (Code-Freeze): _____

0) Ziele & Hypothesen
	•	H1 (Ranking): -\Delta G_\mathrm{insert} und -\Delta G_\mathrm{ads} korrelieren mit \log_{10}(\mathrm{HC50}_\mathrm{RBC}); Ziel: Spearman ρ ≥ 0.50, α=0.05.
	•	H1 (Klassifikation): Lineares Modell aus ΔG-Features trennt toxisch (HC50 ≤ 50 µM) vs. nicht-toxisch (HC50 ≥ 250 µM) mit ROC-AUC ≥ 0.80 (95 %-CI > 0.65).
	•	Heuristik-Cutoff (vorregistriert): „RBC-kritisch“, wenn ΔG_ads < −8 kJ mol⁻¹ und (ΔG_insert ≤ −3 kJ mol⁻¹ oder ΔG‡ < 12 kJ mol⁻¹).

⸻

1) Datensatz
	•	n = 40 Peptide: 20 toxische (HC50 ≤ 50 µM), 20 nicht-toxische (HC50 ≥ 250 µM). Borderline (50–250 µM) werden ausgeschlossen.
	•	Meta-Tabelle hinterlegt: data/hc50_rbc.csv (Spalten: id,sequence,HC50_uM,amidated_Cterm=1,charge@pH7,length,hydrophobic_moment).
	•	Gruppen-Matching: Länge, Nettoladung, hydrophober Moment grob balanciert (±1 SD).

⸻

2) Simulationssystem (für alle identisch)
	•	Box: 10×10×14 nm; Z-Puffer ≥ 4 nm Wasser pro Seite.
	•	Membran: RBC-außenähnlich, asymmetrisch; CHOL = 45 mol % (Primärsystem).
	•	Ionen: 150 mM NaCl. T = 310 K, P = 1 bar (semi-isotrop).
	•	Force Field: MARTINI-3; Δt = 20 fs; r_coul = r_LJ = 1.1 nm.
	•	Termini: N-Terminus +1, C-Terminus amidiert (0).
	•	z-Referenz: äußere Phosphat-Ebene; Membran-COM in z schwach restrain (nur Referenzstabilisierung).

⸻

3) Monomer-PMF (Adsorption/Insertion entlang z)

Collective Variable (CV): z-Abstand zwischen Peptid-COM und äußerer Phosphat-Ebene; positiv = Wasser, negativ = Richtung Membran.

Umbrellas & Laufzeiten
	•	Fensterzentren: z = +2.8\,..\,−1.4 nm in 0.2-nm-Schritten (16 Fenster).
	•	Federkonstante: k = 900 kJ mol⁻¹ nm⁻² (800–1000 zulässig).
	•	Produktionszeit: 20 ns/Fenster.
	•	Replikate: 2 unabhängige Reps pro Peptid (unterschiedliche Start-Orientierungen/Seeds).
	•	Starts: Peptid im Bulk (z≈+1.5 nm); XY per Poisson-Disk (min 1.5–2.0 nm Abstand). 50 ps NVT mit Flat-Bottom-Z-Restraint (±0.2 nm), dann frei.
	•	Sampler: GROMACS + umbrella; Auswertung: MBAR (Primär; WHAM nur als QC).

QC-Schwellen (Gate-Kriterien)
	•	Overlap: Nachbar-Overlap O_{i,i+1} > 0.10.
	•	Effektive Proben: ≥ 200 unkorrelierte Samples pro Fenster (aus IAT).
	•	Halbierungs-Stabilität: |PMF\text{erste Hälfte} − PMF\text{zweite Hälfte}| < 2 kJ mol⁻¹ an allen 16 Mittelpunkten.
	•	Replikat-Konsistenz: |Feature\text{Rep1} − Feature\text{Rep2}| ≤ 2 kJ mol⁻¹.
	•	Wenn QC fällt:
	•	Verlängere nur betroffene Fenster auf 30 ns.
	•	Falls O_{i,i+1} \le 0.10 bleibt → zusätzliche Zwischen-Fenster (Δz = 0.1 nm) anlegen.

PMF-Definition & Glättung
	•	Nullpunkt: G_\text{bulk}=0 als Mittel von z ≥ +2.4 nm.
	•	Interpolation: PCHIP (shape-preserving cubic Hermite) auf 0.01-nm-Raster — keine oszillierende Smoothing-Splines.
	•	Suchbereiche:
	•	z_\text{surf,min} ∈ [−0.4, +0.8] nm
	•	z_\text{head,min} ∈ [−1.4, −0.2] nm
	•	z^\ddagger: Maximum zwischen den beiden Minima.

Primäre Features (mit Unsicherheiten)
	•	ΔG_ads = G(z_\text{surf,min}) − G_\text{bulk}
	•	ΔG_insert = G(z_\text{head,min}) − G(z_\text{surf,min})
	•	ΔG‡ = G(z^\ddagger) − G(z_\text{surf,min})
	•	SE/CI: 1000-facher Block-Bootstrap (blöckeweise über Frames je Fenster).
	•	Pooling über Reps: gewichtetes Mittel mit w_i=1/\mathrm{SE}_i^2.

⸻

4) Statistik (vorregistriert)

Outcomes & Modelle
	•	Ranking: Spearman-ρ von -\Delta G_\mathrm{insert}, -\Delta G_\mathrm{ads}, -\Delta G^\ddagger gegen \log_{10}(\mathrm{HC50}). 95 %-CI via BCa-Bootstrap (10 000 Resamples).
	•	Klassifikation: logistische Regression (L2-ridge) mit Features \Delta G_\mathrm{ads}, \Delta G_\mathrm{insert}, \Delta G^\ddagger.
	•	Validierung: Nested 5× Stratified CV (innen λ-Grid, außen Performance).
	•	Kennzahlen: ROC-AUC (DeLong-CI), PR-AUC, Brier-Score, Kalibrations-Slope/Intercept.
	•	Signifikanztests: DeLong (AUC > 0.5), permutationstest-basiert (10 000 Permutationen) als Robustheits-Check.
	•	Kovariaten-Kontrolle: partielle Rang-Korrelation kontrolliert für Länge, Nettoladung, hydrophoben Moment.

Erfolgskriterien (Go/No-Go)
	•	Korrelation: ρ ≥ 0.50 (p < 0.05) oder 95 %-CI der ρ-Schätzung ganz > 0.30.
	•	Klassifikation: AUC ≥ 0.80; 95 %-CI-Untergrenze > 0.65.
	•	Heuristik-Cutoff erreicht Sensitivität ≥ 0.75 und Spezifität ≥ 0.70 im äußeren Test.

⸻

5) Sensitivitäten (Robustheit, vorab definiert)
	•	CHOL-Robustheit: Subset von 10 Peptiden erneut bei 35 % CHOL; Ranking-Stabilität: Kendall-τ ≥ 0.6 gegenüber 45 %.
	•	Flex-Robustheit: 5 Peptide ohne Elastic-Network; ΔG-Rangänderung (τ) dokumentieren.
	•	Replikat-Robustheit: leave-one-rep-out auf Feature-Ebene; Abweichung ≤ 2 kJ mol⁻¹.

⸻

6) Reproduzierbarkeit
	•	Seeds: seed = crc32(peptide_id + "_rep" + rep + "_umbrella" + window) (deterministisch).
	•	Dateistruktur: results/{peptide_id}/pmf/rep{1,2}/win_{z}; Logs & mdp im selben Ordner.
	•	Provenance: alle mdp-Files, Topologien, commit-Hashes in results/manifest.json.
	•	Automatischer QC-Report: PDF/HTML mit Overlap-Heatmap, IAT-Tabelle, Halbierungs-Vergleich, PMF ± 95 %-Band.

⸻

7) Optionales Zusatz-Feature (nur falls AUC/ρ < Ziel)
	•	Monomer-Thinning-Assay: 1 Peptid auf Oberfläche, 5×200 ns; Metriken: lokale Δ-Dicke (unter Peptid vs. Fernfeld), Midplane-Wasser, Tail-Order. Aggregation per Median.
	•	Integration: Hinzunahme als drittes Prädiktor-Feature in die logistische Regression; neues AUC-Ziel unverändert.

⸻

8) Power (Begründung)
	•	Korrelation: n=40, erwartete wahre ρ≈0.5 → Power ≈ 0.92 (α=0.05).
	•	Klassifikation: erwartete AUC≈0.80 → Power ≈ 0.99 gegenüber 0.5.
(Rechnungen vorab dokumentiert; Details im Analysis-Notebook notebooks/stats_prereg.ipynb.)

⸻

9) Deliverables
	•	Pro Peptid: PMF-Plot (PCHIP-Interpol., 95 %-Band), Tabelle ΔG_ads, ΔG_insert, ΔG‡ ± SE.
	•	Gesamt: Scatter -\Delta G vs. -\log_{10}\mathrm{HC50}, ROC/PR-Kurven (mit CIs), Kalibrations-Plot, Decision-Curve, finaler Heuristik-Cutoff + Confusion-Matrix.
	•	Supplement: QC-Aggregat (Overlap, IAT, Halbierungs-Stabilität) und Sensitivitäts-Resultate (35 % CHOL, EN-Off).

⸻

10) Abweichungen & Stop-Regeln (vorab)
	•	Fenster-Fehlschlag: Wenn trotz Verlängerung und zusätzlicher Fenster O_{i,i+1} ≤ 0.10, wird das Peptid ausgeschlossen und als „non-evaluable“ dokumentiert.
	•	Instabile Membran (Artefakt): sichtbares Demixing/Kollaps → Run verwerfen, mit neuem Seed einmalig wiederholen.
	•	Maximalbudget pro Peptid: 16×30 ns × 2 Reps (0.96 µs). Überschreitung nicht zulässig.

⸻

11) Vorab festgelegte Software
	•	MD: GROMACS (Version ___), MARTINI-3 Parametrierung (Commit: ___).
	•	Freienergie: pymbar (Version ___) für MBAR; wham nur als QC.
	•	Statistik: Python (numpy/scipy/statsmodels/scikit-learn), Bootstrap & DeLong Implementierungen (Versionen dokumentiert).

⸻

12) Abschlusskriterien
	•	H1 (Ranking) erfüllt?  □ Ja  □ Nein
	•	H1 (Klassifikation) erfüllt?  □ Ja  □ Nein
	•	Heuristik-Cutoff validiert (Sens ≥ 0.75 & Spez ≥ 0.70)?  □ Ja  □ Nein
	•	Robustheit (Kendall-τ bei 35 % CHOL ≥ 0.6)?  □ Ja  □ Nein

⸻

Anhang A: feste Umbrella-Zentren (nm)

[+2.8, +2.6, +2.4, +2.2, +2.0, +1.8, +1.6, +1.4, +1.2, +1.0, +0.8, +0.6, +0.4, +0.2, 0.0, −0.2, −0.4, −0.6, −0.8, −1.0, −1.2, −1.4]
(Bei Rechenbudget-Limit dürfen die beiden äußersten Bulk-Fenster entfallen, sofern Bulk-Plateau stabil ist.)

Anhang B: feste Reporting-Tabelle (CSV-Header)

peptide_id, ΔG_ads_kJmol, ΔG_ads_SE, ΔG_insert_kJmol, ΔG_insert_SE, ΔG_barrier_kJmol, ΔG_barrier_SE, z_surf_min_nm, z_head_min_nm, hc50_uM

⸻

Stichprobe:
TOXIC
Solvia 1
Solvia 8
Solvia 14
Solvia 215
Solvia 164
Solvia 126
Solvia 68
Solvia 32
Solvia 482
Solvia 490
Solvia 515
Solvia 524
Solvia 527
Solvia 617
Solvia 624
Solvia 850
Solvia 858
Solvia 941
Solvia 974
Solvia 1023
Solvia 1045


NON-TOXIC
Solvia 12
Solvia 398
Solvia 1051
Solvia 1125
Solvia 1219
Solvia 1315
Solvia 1343
Solvia 1363
Solvia 1564
Solvia 1587
Solvia 1663
Solvia 1680
Solvia 1684
Solvia 1743
Solvia 1844
Solvia 1941
Solvia 1952
Solvia 1962
Solvia 2012
Solvia 1115
Solvia 794
