# HC50 aus Adsorptions‑Thermodynamik (Wasserseite)

Dieses Dokument definiert die Berechnung von Kp_ads und HC50 in SOLVIA.

- Temperatur: 310 K (konfigurierbar unter `pmf.analysis.simulation.temperature`)
- Adsorptionsband: Wasserseite (z ≥ 0) bis `feature_params.ads_max_z`, Energieband `W ≤ W_min + 3 kJ/mol`, zusätzlich nur innerhalb der Membran‑Bounds (`adsorption.membrane.z_lo/z_hi`) und nur Bins mit ausreichender Statistik (`plot.pmf_min_bin_count`).
- Integrand: `exp(-β W(z)) − 1` auf dem Adsorptionsband, Trapezregel → `Kp_ads [nm]`.
- Bilayer‑Normalisierung: `kp_eff_nm = bilayer_factor * kp_nm` (Default `2.0`).
- Konversion zu Konzentration (M ⇄ nm⁻³): `NA = 0.602214`.
- HC50‑Band: `hc50_uM_range = [Γ*_(max)/(NA·Kp_eff), Γ*_(min)/(NA·Kp_eff)] * 1e6` mit `Γ* = 1/(A_lipid · L/P*)` und `L/P* ∈ [154, 515]`.

In der Analyse wird das Adsorptionsband im PMF‑Plot schattiert. ΔG_insert bleibt optional (Diagnose) und fließt nicht in Ranking/Klassifikation ein.

Siehe Implementierung in `analysis/hc50.py` und Anbindung in `scripts/analysis/pmf_mbar_analysis.py`.

