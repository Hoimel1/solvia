# HC50 – FAQ

- Warum nur Wasserseite?
  - Ziel ist die Adsorptions‑Thermodynamik; Einlagerung in die Membran (z < 0) ist peptidspezifisch und nicht notwendig für die Ableitung von HC50.

- Was passiert bei flachen PMFs?
  - Wenn `Kp_eff → 0`, divergiert HC50 → ∞. Das Ergebnis wird als nicht toxisch bzw. „undetermined“ klassifiziert, abhängig von der Policy.

- Einfluss von `area_per_lipid` und `lp_star_range`?
  - Beide skalieren Γ*. Größere `A_lipid` oder größere `L/P*` führen zu niedrigeren Γ* und damit niedrigeren HC50‑Werten. Die Standardwerte sind in der Projekt‑Config festgelegt und können per CLI überschrieben werden.

- Was macht der `bilayer_factor`?
  - Skaliert `Kp` auf die Bilayer‑Geometrie. Standard ist `2.0` (beide Leaflets). Bei explizit bilateraler Auswertung kann `1.0` sinnvoll sein.

