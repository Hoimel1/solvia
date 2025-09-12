Simulationsgestütztes Konzept zur Vorhersage der hämolytischen Aktivität von Peptiden

Ziele und Hypothesen

H1 (Ranking): Peptide mit höherer Membranbindungsaffinität – erkennbar an stärker negativen Freie-Energie-Werten (ΔG) für Adsorption oder Insertionsneigung – sollen stärker hämolytisch sein. Konkret wird erwartet, dass -ΔGinsert und -ΔGads positiv mit der Log10 des HC50-Werts korrelieren (je negativer die Energie, desto niedriger der HC50-Wert). Als Erfolgskriterium gilt Spearman ρ ≥ 0,5 (p < 0,05) für diese Korrelation, mit einem 95%-Konfidenzintervall, dessen untere Grenze > 0,30 liegt.

H1 (Klassifikation): Die ΔG-Parameter sollen es ermöglichen, Peptide als “toxisch” vs. “nicht-toxisch” zu unterscheiden. Ein einfaches lineares Modell (logistische Regression) basierend auf ΔG-Features (siehe unten) soll toxische Peptide (HC50 ≤ 50 µM) von nicht-toxischen (HC50 ≥ 250 µM) mit hoher Güte trennen (Ziel: ROC-AUC ≥ 0,80, unteres 95%-CI > 0,65).

Heuristische Schwellen: Zusätzlich wurde vorab ein praktischer Cutoff definiert, um „RBC-kritische“ Peptide schnell zu identifizieren: ΔGads < –8 kJ/mol und (ΔGinsert ≤ –3 kJ/mol oder ΔG‡ < 12 kJ/mol). Peptide, die diese Kriterien erfüllen, gelten hypothetisch als hämolytisch riskant. Diese Daumenregel wird im Verlauf validiert (erwartet: Sensitivität ≥ 75%, Spezifität ≥ 70%). Sie dient als ergänzende Prüfung, ob eine simple Regel mit festen Schwellen mit dem komplexeren Modell mithalten kann.

Peptid-Datensatz

Wir betrachten 42 Peptide unterschiedlicher Sequenz, für die experimentelle Hämolyse-Daten vorliegen. Sie teilen sich auf in 21 toxische (HC50 ≤ 50 µM) und 21 nicht-toxische (HC50 ≥ 250 µM) Peptide. Fälle im Zwischenbereich (50–250 µM) wurden ausgeschlossen, um eine klare Dichotomie zu erhalten. Wichtig ist, dass die beiden Gruppen in ihren physikochemischen Eigenschaften balanciert sind: Länge, Nettoladung (bei pH 7) und hydrophober Moment unterscheiden sich nicht systematisch zwischen toxisch und nicht-toxisch. Dies verhindert, dass triviale Kovariaten die Ergebnisse verzerren (z.B. „alle langen Peptide sind toxisch“).

Alle Peptide werden in identischer chemischer Form betrachtet: N-Terminus protoniert (+1) und C-Terminus amidiert (neutral). Diese Terminierung ist üblich für synthetische Peptide und stellt sicher, dass die Ladungen konsistent sind. Eine Datendatei (CSV) listet für jedes Peptid ID, Sequenz, HC50-Wert sowie die berechneten Eigenschaften (Amidierung, Nettoladung, Länge, hydrophober Moment). Diese Datenbasis (data/hc50_rbc.csv) dient später auch der Analyse (z.B. Korrelationen mit Länge etc.).

(Beispielhafte Peptid-IDs aus der Auswahl – toxisch: Solvia 1, 8, 14, …; nicht-toxisch: Solvia 12, 398, 1051, … – die vollständige Liste mit Sequenzen ist in der CSV-Datei enthalten.)

Stichprobe (42 Peptide; 21 toxisch / 21 nicht-toxisch):
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

Simulationssystem und Randbedingungen

Um die Peptid-Membran-Interaktion realistisch abzubilden, verwenden wir ein einheitliches Coarse-Grain (CG) Simulationssystem, das die äußere Membranblatt einer roten Blutzelle (RBC) nachstellt. Wesentliche Merkmale dieses Systems:
	•	Membranaufbau: Wir modellieren eine asymmetrische Lipiddoppelschicht, die der Außenleaflet einer RBC-Membran entspricht. Die äußere Schicht besteht überwiegend aus zwitterionischen Phospholipiden (Phosphatidylcholin und Sphingomyelin ähnlich) sowie einem hohen Cholesterin-Anteil (~45 mol%), da Erythrozytenmembranen etwa 40–50 mol% Cholesterin enthalten ￼. (Zum Vergleich: In tierischen Zellmembranen befinden sich PC und SM hauptsächlich im äußeren Leaflet, während PS und PE innen liegen ￼. Cholesterin ist asymmetrisch verteilt, mit leichter Bevorzugung der Außenleaflet in RBCs ￼.) Die innere Leaflet wird in unserem Modell ebenfalls mit neutralen Lipiden bestückt, da wir primär an der Peptid-Wechselwirkung mit der Außenseite interessiert sind – die genaue innere Zusammensetzung ist nachrangig, solange die Gesamtmembran stabil ist.
	•	Simulationsbox: Etwa 10 × 10 × 14 nm mit der Membran horizontal in der Mitte. Ober- und unterhalb der Membran befinden sich jeweils ca. 4 nm Wasserschicht als Puffer. Das System enthält ~400 Lipide (≈200 pro Leaflet) – genug für ein Membranpatch ohne übermäßige Krümmung. Ionen: 150 mM NaCl wird hinzugefügt, um physiologische Ionenstärke nachzuahmen.
	•	Kraftfeld: Wir nutzen das aktuelle MARTINI 3 coarse-grained Kraftfeld für Lipide und Peptide ￼. Dieses ist gegenüber Vorgängerversionen verbessert in Bezug auf Interaktionsbalancen und ermöglicht akkuratere Vorhersagen z.B. von Partitionierung und Protein-Lipid-Interaktionen ￼. Peptide werden auf CG-Ebene gemäß ihrer Sekundärstruktur abgebildet (für Helices idRRA: standardmäßig mit Elastic Network (EN), siehe Robustheitsanalyse).
	•	Bedingungen: Temperatur 310 K (Körpertemperatur) und Druck 1 bar (semi-isotrop, d.h. die x-y-Fläche kann sich unabhängig von der z-Dimension anpassen, um den korrekten Membrandruck zu gewährleisten). Der Zeitschritt beträgt 20 fs. Van-der-Waals- und Coulomb-Wechselwirkungen haben einen Cutoff von ca. 1,1 nm (Coulomb mit Reaction-Field Behandlung jenseits des Cutoffs, gemäß MARTINI-Empfehlungen).
	•	Leaflet-Definition (fix): Die Zuordnung des Außenleaflets ist über eine statische Index- oder Mapping-Datei festgelegt (z.B. `index_leaflets.ndx` mit Gruppen `OuterPO4`/`InnerPO4`). Eine heuristische Online-Erkennung (Median-z) wird nicht verwendet, damit die Definition des Außenleaflets zwischen Systemen/Replikaten absolut stabil bleibt.
	•	Referenzrahmen: Um die Membran als stabilen Bezugspunkt zu haben (wichtig für die nachfolgende z-Koordinate), verwenden wir optional sehr schwache harmonische z-Restraints auf den Lipid-COM (Topologie-Schalter `-DPOSRES_LIPID_Z`, nur z-Komponente, ~10 kJ·mol⁻¹·nm⁻²). Dies verhindert ein langsames Driften des gesamten Bilayers, ohne die lokalen Fluktuationen zu behindern. Die Peptide selbst erhalten keine Positionsrestrains (außer den Umbrella-Biases).

Vor jedem Produktionslauf wird das System energie-minimiert und sorgfältig equilibriert. Jedes Peptid wird einzeln mit einer Membran simuliert, um Monomer-Membran-Wechselwirkungen zu isolieren. Mehrere Peptide gleichzeitig würden Aggregation oder Konkurrenz um Lipidbindestellen verursachen, was wir hier ausschließen (unser Fokus liegt auf der initialen Monomer-Bindung).

Monomer-Umbrella Sampling entlang der Membrannormalen

Ziel: Ermittlung des Freien-Energie-Profils (Potential of Mean Force, PMF) für die Adsorption und ggf. partielle Insertion jedes Peptids in die Membran. Die Reaktionskoordinate ist dabei die z-Position des Peptid-Schwerpunkts relativ zur Membranoberfläche.
	•	Definition der Koordinate (z): Nullpunkt z = 0 ist die **lokale** z-Position der Phosphatgruppen im Außenleaflet unter dem Peptid, nicht der globale Leaflet-COM. Praktisch wird dafür eine dynamische Patch-Referenz aus `OuterPO4`-Atomen innerhalb eines Zylinders (Radius 1.8–2.0 nm) um den projizierten Peptid-COM (x–y) gebildet. Alternative (falls verfügbar): PLUMED-Definition einer zylindrischen Referenz. Positive z: Peptid im Wasser; negative z: Eindringen Richtung Membraninneres.
	•	Umbrella-Ansatz: Variable Fensterdichte von z = +2.8 nm (Bulk) bis z = −2.0 nm (Kopfgruppen). Dichte Zone von +0.6 bis −0.6 nm mit 0.15-nm Schritten, außerhalb 0.2 nm. Ziel-Federkonstante k = 900 kJ·mol⁻¹·nm⁻² (zulässig 800–1100). Eine **adaptive Verdichtung** stellt eine minimale Nachbar-Überlappung O ≥ 0.20 sicher (Script fügt Zwischenfenster ein, wenn Δz > Δz_max(k,T,O)).
	•	Initialisierung/Replikate: Pro Peptid zwei Replikate mit unterschiedlichen Startorientierungen (R1 ~ parallel, R2 geneigt/orthogonal). Seeds (`ld-seed`/`gen-seed`) werden deterministisch aus Peptid-ID und Replikatsindex abgeleitet. Anstelle eines harten „snap to center“ erfolgt eine kurze **SMD-Leiter** vom Bulk zum Ziel-z (5–10 ns total, k=300–500), dann **Pre-Equil pro Fenster** 1 ns bei finalem k.
	•	Simulationslaufzeit pro Fenster: Produktion 20 ns (bei Bedarf bis 30 ns). Pull-Outputs werden mit ≤1 ps aufgezeichnet (für Overlap/ESS/QC). Bereits während der Läufe überwacht die Pipeline die Overlap-Metriken und fügt ggf. Mid-Points ein.

PMF-Berechnung und Merkmalsableitung

Nach Abschluss aller Umbrella-Simulationen für ein Peptid werden die Daten zusammengeführt, um das Freie-Energie-Profil G(z) zu berechnen:
	•	Auswertung mittels WHAM (Standard) und optional MBAR: Standardmäßig verwenden wir WHAM (Weighted Histogram Analysis), um aus den Umbrella‑Histogrammen das kontinuierliche PMF zu bestimmen. MBAR kann optional als ergänzende Methode eingesetzt werden. Die Bulk‑Wasser‑Referenz wird so gewählt, dass G(z) im Bulk auf 0 kJ/mol normiert wird: konkret setzen wir den Mittelwert von G(z) für z ≥ +2,4 nm als Nullniveau (dort ist das Peptid vollständig gelöst, fern der Membran). Negative freie Energien bedeuten also eine energetische Bevorzugung des Peptids in Membrannähe verglichen mit der Bulk‑Lösung.
	•	Rohes Profil und Glättung: Das erhaltene PMF – G als Funktion von z – kann noch statistisches Rauschen aufweisen. Wir glätten es leicht mittels stückweiser kubischer Hermite-Interpolation (PCHIP) auf einem feinen z-Gitter (0,01 nm Schritt). PCHIP ist formtreu, d.h. sie erzeugt keine überschießenden Maxima/Minima und behält die Monotonie zwischen den Stützstellen bei. Wir vermeiden aggressives Glätten, um keine realen Barrieren oder Minima zu verwischen – die Interpolation dient vor allem der genaueren Bestimmung der Extremstellen (Minima/Maxima).
	•	Identifikation charakteristischer Punkte: Aus dem glatten Profil lesen wir folgende Kenngrößen für jedes Peptid ab:
	•	Oberflächenminimum (zsurf,min): Das erste lokale Minimum nahe der Membranoberfläche, typischerweise um z ≈ 0 bis +0,5 nm. Dieser Zustand entspricht einem adsorbierten Peptid, das an der Membran anliegt, aber nicht tief eingedrungen ist.
	•	Kopfgruppen-Minimum (zhead,min): Ein mögliches zweites Minimum tiefer in der Membran (etwa z ≈ –0,5 bis –1,0 nm, also in den Lipidkopfgruppen). Nicht alle Peptide weisen dies auf. Falls vorhanden, deutet es darauf hin, dass das Peptid teilweise in die Membran eindringt und dort einen energetisch günstigen Zustand findet.
	•	Barriere (z‡): Ein lokales Maximum zwischen den obigen Minima (falls zwei Minima existieren). Dies repräsentiert die energetische Übergangsbarriere zwischen oberflächlicher Bindung und tieferem Eindringen. Ist nur ein Minimum vorhanden (kein tieferer Einlagerungszustand), betrachten wir die Barriere zwischen diesem Minimum und dem Bulk.
	•	Berechnung der ΔG-Features: Aus diesen Punkten bestimmen wir die quantitativen Deskriptoren:
	•	ΔGads = G(zsurf,min) – Gbulk: Freie Energie der Adsorption aus dem Bulk ans Membraninterface. (Da Gbulk = 0 gesetzt ist, ist dies einfach G am Oberflächenminimum.) Negative Werte bedeuten spontane, energetisch begünstigte Bindung an die Membran.
	•	ΔGinsert = G(zhead,min) – G(zsurf,min): Zusätzliche Freie-Energie-Änderung, wenn das Peptid von der Oberfläche weiter in die Membran eindringt. Ein negativer Wert bedeutet, dass das Peptid sogar energetisch davon profitiert, tiefer in die Lipidschicht zu gehen – ein starker Hinweis auf Membranaktivität. Positive Werte (oder kein zweites Minimum) bedeuten, dass das Peptid an der Oberfläche “festgehalten” wird und ein tieferes Eindringen ungünstig ist.
	•	ΔG‡ = G(z‡) – G(zsurf,min): Barrierehöhe zwischen dem Oberflächen- und dem tieferen Zustand. Gibt an, welche Energiehürde das Peptid überwinden muss, um von der Oberfläche in die Membran einzudringen. Eine niedrige Barriere (oder gar keine, falls ΔGinsert ≤ 0) lässt vermuten, dass das Peptid leicht in die Membran gelangen kann (häufig korreliert mit pore-forming Tendenzen). Hohe Barrieren deuten auf ein in der Oberfläche verharrendes Peptid hin.
	•	Unsicherheitsabschätzung: Für jedes dieser Features berechnen wir statistische Unsicherheiten. Dazu nutzen wir bootstrapping der Umbrella-Datensätze: Aus jedem Umbrella-Window werden mit Zurücklegen Stichproben (in Blöcken, um Autokorrelation zu berücksichtigen) gezogen und 1000 „synthetische“ PMFs rekonstruiert. Aus der Streuung der resultierenden ΔG-Werte ergibt sich der Standardfehler (SE) und ein 95% Konfidenzintervall. Dies wird pro Replikat gemacht.
	•	Replikate-Mittelung: Liegen für ein Peptid zwei Replikate vor, werden die Feature-Werte anschließend zu einem konsolidierten Wert gemittelt. Wir verwenden ein gewichtetes Mittel nach Präzision: Gewicht w = 1/SE² pro Replikat. So trägt der genauere Wert mehr bei. Voraussetzung ist allerdings, dass die beiden Replikate konsistente Ergebnisse liefern (siehe nächster Abschnitt zur Qualitätskontrolle). Bei Diskrepanzen würde zunächst nachgebessert, statt stumpf zu mitteln. Die gemittelten ΔGads, ΔGinsert, ΔG‡ mit zugehörigem Fehler (kombiniert) bilden schließlich den Deskriptorensatz für jedes Peptid.

Qualitätskontrolle (QC) und Abbruchkriterien

Um sicherzustellen, dass die PMF-Daten vertrauenswürdig und vergleichbar sind, werden mehrere QC-Kriterien angelegt. Alle Umbrella-Datensätze eines Peptids müssen diese erfüllen, andernfalls wird nachgebessert oder das Peptid aus der Auswertung ausgeschlossen. Die wichtigsten Kriterien:
	•	Fenster-Überlappung: Nachbarüberlappung ≥ 10 % (Histogramm-basierte Schnittmenge der projizierten z-Verteilungen). Falls <10 %: automatisches Zwischenfenster (0.1–0.15 nm) und Nachsimulation.
	•	Effektive Stichprobengröße (ESS): Blockweise Autokorrelation; ESS ≥ 200 unkorrelierte Frames pro Fenster (bei 1-ps-Sampling ≈ 200 ps). Falls <200: Fenster um 10 ns verlängern.
	•	Halbzeit-Konvergenz: ΔG am Fensterzentrum zwischen erster/zweiter Hälfte |Δ| ≤ 2 kJ·mol⁻¹. Falls verletzt: +10 ns für betroffene Fenster.
	•	Replikat-Konsistenz: |Δ(ΔG_ads)| ≤ 2 kJ·mol⁻¹ und |Δ(ΔG_insert)| ≤ 2 kJ·mol⁻¹. Bei Verstoß: gezielte Nachverdichtung/Verlängerung in der divergierenden z-Region.

Korrekturmaßnahmen: Sollten eine oder mehrere der obigen Prüfungen fehlschlagen, wird folgendermaßen vorgegangen:
	•	Ungenügende Überlappung: Zusätzliche Umbrella-Fenster einfügen und 20 ns simulieren (wie oben beschrieben). Danach MBAR erneut ausführen.
	•	Schlechte Konvergenz oder wenige unabhängige Samples: Verlängerung der betreffenden Fenster (z.B. auf 30 ns) oder falls generalisiert, aller Fenster eines Replikats auf 30 ns. Erneuter Halbzeit-Check.
	•	Replikat-Divergenz: Identifikation der z-Region, die die Divergenz verursacht (z.B. vergleicht man die PMF-Kurven, oft ist es ersichtlich, wo sie auseinandergehen). Für diesen Bereich zusätzliche Fenster oder Laufzeit hinzufügen.

Wir begrenzen den maximalen Simulationsaufwand pro Peptid jedoch auf ca. 0,96 µs (entspricht 16 Fenster × 30 ns × 2 Replikate). Falls ein Peptid trotz aller Nachbesserungen kein stabiles PMF liefert (z.B. persistente Konvergenzprobleme oder chaotisches Verhalten), wird es als “nicht auswertbar” deklariert und aus der Studie ausgeschlossen. Solche Fälle werden im Bericht transparent gemacht, inklusive Grund (z.B. „Peptid X zeigt in CG-Simulation anhaltende Aggregation an der Membran, keine stabile Monomer-Adsorption erreichbar“ o.ä.).

Glücklicherweise zeigen Vorstudien, dass bei den meisten Peptiden das oben skizzierte Sampling ausreichend ist. Durch die standardisierte Protokollierung (wir erfassen z.B. Overlap-Matrizen, Autokorrelationszeiten und PMF-Hälftenplots für jedes Peptid) stellen wir sicher, dass am Ende nur qualitativ einwandfreie Daten in die Analyse eingehen. Diese QC-Berichte werden auch in den Supplementary Materials bereitgestellt, um die Nachvollziehbarkeit zu gewährleisten.

Statistische Analyse: Korrelation und Klassifikation

Mit dem finalen Datensatz der ΔG-Features und den experimentellen HC50-Werten prüfen wir nun die Hypothesen.

Korrelationsanalyse (Ranking)

Zunächst wird untersucht, ob die Peptide mit stärker negativer freier Energie tatsächlich tendenziell toxischer sind. Dafür berechnen wir Spearman-Rangkorrelationskoeffizienten (ρ) zwischen den ΔG-Parametern und der gemessenen Hämolyse. Im Einzelnen:
	•	ρ(-ΔGinsert, log10(HC50)) – Hauptfokus, da wir erwarten, dass die Insertionsneigung ein besonders guter Prädiktor ist.
	•	Zum Vergleich auch ρ(-ΔGads, log10(HC50)) und ρ(-ΔG‡, log10(HC50)).

Hierbei verwenden wir -ΔG, weil wir annehmen, dass ein Peptid mit z.B. ΔGinsert = -5 kJ/mol (stark begünstigte Einlagerung) einen kleineren HC50 (hoch toxisch) hat als eines mit ΔGinsert = +5 kJ/mol. Eine negative ΔG entspricht also potenziell hoher Toxizität (niedriger HC50), was in einer negativen Pearson-Korrelation resultieren würde – Spearman nimmt uns das Vorzeichenproblem, indem wir -ΔG nutzen.

Wir quantifizieren die Unsicherheit der Korrelation durch Bootstrap (10.000-faches Resampling der Peptid-Paare) und geben ein 95%-Konfidenzintervall für ρ an (bias-corrected, accelerated). Außerdem testen wir die Signifikanz gegen H0: ρ = 0 (kein Zusammenhang) bei α = 0,05. Erfolgreiche Validierung von H1 (Ranking) läge vor, wenn ρ deutlich positiv ist und statistisch signifikant (siehe Kriterien oben).

Um auszuschließen, dass beobachtete Korrelationen bloß indirekt über Peptideigenschaften zustande kommen, führen wir auch partielle Korrelationen durch. Beispielsweise berechnen wir ρ(-ΔGinsert, log HC50 | Kontrolle für Peptidlänge) – also den Zusammenhang, bereinigt um die Effekt von Länge. Ähnliches für Nettoladung oder hydrophoben Moment. Aufgrund der balancierten Datensatzzusammensetzung erwarten wir keine drastische Änderung, aber diese Analyse kann bestätigen, dass ΔG tatsächlich ein eigenständiger Faktor ist.

Wir vergleichen unsere Ergebnisse mit der Literatur: Spaller et al. (2013) fanden experimentell, dass die freie Bindungsenthalpie an POPC-Liposomen gut mit der Hämolyseaktivität korreliert ￼. Zhao et al. (2013) zeigten in einer kombinierten Simulations-/Experimentstudie zu antimikrobiellen Peptiden, dass PMF-Profile (für Bakterien- vs. RBC-Membran) Rückschlüsse auf Wirksamkeit und Toxizität erlauben ￼ ￼. Wir knüpfen hier an und quantifizieren diese Beziehung für unseren Datensatz.

Klassifikationsanalyse (Toxisch vs. Nicht-toxisch)

Neben dem graduellen Zusammenhang interessiert uns, ob wir auch eine klare binäre Vorhersage treffen können. Dazu entwickeln wir ein einfaches logistisches Regressionsmodell mit den ΔG-Features als Prädiktoren.
	•	Features: ΔGads, ΔGinsert, ΔG‡ (alle drei werden angeboten; das Modell kann unwichtige über Koeffizienten ~0 effektiv ignorieren). Gegeben der kleinen Stichprobe (n=42) und zur Vermeidung von Überanpassung nutzen wir eine L2-Regularisierung (Ridge). Die Regularisierungsstärke λ wird via Cross-Validation abgestimmt.
	•	Cross-Validation: Wir führen eine 5-fache Kreuzvalidierung durch, nested: In der inneren Schleife werden verschiedene λ-Werte ausprobiert, in der äußeren Schleife wird die Performance auf einem jeweis gehaltenen Test-Fold gemessen. So erhalten wir robuste Schätzungen der Generalisierungsleistung.
	•	Metriken: Hauptmetrik ist die AUC (Area Under the ROC Curve), die die Trennschärfe des Modells zusammenfasst. Wir streben AUC ≥ 0,80 an. Zusätzlich betrachten wir die PR-Kurve (Precision-Recall), die bei leicht unausgeglichenen Klassen aufschlussreich ist (hier sind Klassen 20/20, also ausgeglichen, aber PR zeigt z.B. den positiven Vorhersagewert bei unterschiedlichen Schwellen). Wir berechnen weiterhin den Brier-Score für die Kalibrierung der prognostizierten Wahrscheinlichkeiten. Ein Kalibrationsplot (vorhergesagte vs. tatsächliche Eintrittswahrscheinlichkeit der Toxizität) gibt Hinweise, ob das Modell z.B. systematisch zu konservativ oder optimistisch schätzt.
	•	Modellgüte und Signifikanz: Für die gemittelte ROC-AUC über die CV-Folds bestimmen wir ein 95%-Konfidenzintervall (mittels DeLongs Verfahren). Auch ein Permutationstest (Zufallsverteilung der AUC bei 10.000 zufälligen Label-Shuffles) kann zeigen, ob die beobachtete AUC weit über Zufall liegt. Die Nullhypothese (AUC = 0,5) sollte klar abgelehnt werden (Ziel p < 0,01).
	•	Schwellenwert-Optimierung: Aus der ROC-Kurve jedes Folds wählen wir den optimalen Trennwert nach dem Youden-Index (maximiert Sensitivität + Spezifität). Daraus leiten wir eine durchschnittliche Schwelle ab, um das Modell in eine Ja/Nein-Entscheidung umzusetzen. Wir berichten mittlere Sensitivität und Spezifität bei diesem Cutoff.
	•	Heuristik-Vergleich: Unabhängig davon prüfen wir die vorgegebene Heuristik (ΔGads < –8 und ΔGinsert ≤ –3 oder ΔG‡ < 12). Diese Regel wenden wir auf die Peptid-Daten an und berechnen, wie viele toxische und nicht-toxische sie jeweils korrekt einordnet. So erhalten wir die Sensitivität/Spezifität dieser festen Regel und können vergleichen, ob sie dem optimierten Modell nahekommt. Falls ja, wäre das praktisch bedeutsam, denn dann könnte man in Zukunft diese einfache Regel zur Vorauswahl verwenden.

Erfolgskriterium H1 (Klassifikation) ist erfüllt, wenn unser Modell eine ROC-AUC ≥ 0,80 erreicht und das untere Ende des 95%-CI über 0,65 liegt. Das würde bedeuten, dass mit hoher Wahrscheinlichkeit eine echte Trennfähigkeit vorliegt, die weit besser als Zufall ist. Zusätzlich sollte die Performance in allen genannten Metriken konsistent gut sein (z.B. PR-AUC hoch, Brier-Score niedrig, etc.).

Robustheits- und Sensitivitätstests

Wir wollen zeigen, dass unsere Schlussfolgerungen nicht empfindlich von spezifischen Modellannahmen oder Parametern abhängen. Daher sind folgende Sensitivitätsanalysen geplant:
	•	Membran-Cholesterin-Gehalt: Da RBC-Membranen etwa 40–50% Cholesterin enthalten, haben wir 45% im Basismodell angenommen. Um zu prüfen, ob die Ergebnisse robust gegen leichte Änderungen sind, simulieren wir einen Subset von 10 Peptiden erneut in einem System mit 35 mol% Cholesterin (ersetzt durch mehr PC/SM). Wir berechnen für diese Peptide erneut ΔGads und ΔGinsert und vergleichen die Rangordnung der Peptide zwischen 45%- und 35%-System. Konkret berechnen wir Kendalls τ zwischen den Rangreihen (nach -ΔGinsert sortiert) in beiden Bedingungen. Erwartung: τ ≥ ~0,6, d.h. die relative Reihenfolge ändert sich kaum. Sollte die Rangfolge vollkommen anders ausfallen (τ nahe 0), würde das bedeuten, dass unser Prädiktor stark vom Cholesterin-Gehalt abhängt. Das wäre problematisch, ist aber unwahrscheinlich – i.d.R. führt weniger Cholesterin zu etwas leichterer Membranpenetration für alle Peptide (abs. ΔGinsert wird etwas negativer), ohne dass Peptide ihre relative Position tauschen.
	•	Elastic Network (EN) im CG-Forcefield: MARTINI 3 nutzt oft ein Elastic Network, um Proteine (hier Peptide) steifer in ihrer experimentellen Struktur zu halten. Das könnte theoretisch die Peptid-Flexibilität einschränken und die Interaktion beeinflussen. Wir testen 5 Peptide exemplarisch ohne EN (d.h. voll flexible CG-Peptide) und führen deren Umbrella-Sampling nochmals durch. Anschließend prüfen wir, ob sich die ΔG-Profile substanziell unterscheiden. Insbesondere vergleichen wir wieder die Rangordnung dieser 5 Peptide untereinander. Erwartung: Die absoluten Werte könnten geringfügig variieren (z.B. ein sehr flexibles Peptid könnte etwas tiefere Insertion finden), aber die Rangfolge bleibt gleich. Falls ein Peptid nur mit EN eine tiefe Insertion zeigt und ohne EN nicht mehr, wäre das ein Hinweis, dass unsere Resultate für dieses Peptid vom Strukturmodell abhängen – das würden wir im Bericht diskutieren. Insgesamt erwarten wir aber Robustheit, da viele antimikrobielle Peptide in Realität auch helikale Strukturen an Membranen annehmen und MARTINI+EN dies gut repräsentiert.
	•	Replikat-Wege: Schließlich überprüfen wir, ob nicht ein einzelnes Replikat ungewollt die Ergebnisse verzerrt. Dazu führen wir eine leave-one-replicate-out Analyse durch: Wir berechnen einmal alle Korrelations- und Klassifikationsergebnisse nur mit den Feature-Werten aus Replikat 1 (d.h. wir ignorieren Replikat 2) und einmal nur mit denen aus Replikat 2. Beide sollten zu ähnlichen Schlüssen kommen innerhalb der statistischen Unsicherheit. Wenn z.B. die Korrelation nur in Replikat 1 signifikant wäre, in Replikat 2 aber nicht, wäre Vorsicht geboten – das würde einen starken Zufallseinfluss nahelegen. Ideal ist, wenn beide getrennt schon eine deutliche Tendenz zeigen (wenn auch mit größerem Fehler, da halbe Daten). Wir erwarten im Erfolgsfall, dass die Replikate konsistent sind, was zusätzlich das Vertrauen stärkt.

Alle diese Sensitivitätschecks werden ebenfalls quantitativ dokumentiert (z.B. in Form von Tabellen oder Scatterplots im Anhang). So können wir zeigen, dass unser Verfahren robust ist und nicht z.B. nur für exakt 45% Cholesterin funktioniert.

Optional: Monomer-“Thinning” Assay als Ergänzung

(Dieser Schritt wird nur verfolgt, falls die Hauptziele mit den bisherigen Features nicht erreicht werden. Er ist hier der Vollständigkeit halber prä-registriert, um im Bedarfsfall darauf zurückgreifen zu können, ohne das Studiendesign nachträglich abzuändern.)

Falls die Korrelationen und Klassifikationsergebnisse wider Erwarten schwach ausfallen (z.B. Spearman ρ < 0,3 oder ROC-AUC < 0,7), erwägen wir, einen zusätzlichen Membranperturbations-Deskriptor einzuführen. Die Hypothese: Einige Peptide verursachen Membranschäden nicht (nur) dadurch, dass sie tief eindringen, sondern indem sie bereits als Monomer auf der Oberfläche die Lipidordnung stören. Um das zu erfassen, führen wir einen unbiased Monomer-Surface-Assay durch:
	•	Jedes Peptid wird einzeln auf ein Membranpatch (wie oben, 45% Chol) gelegt, ohne Zwang (außer ggf. anfänglich positioniert an der Oberfläche), und für 5 × 200 ns MD simuliert (fünf Wiederholungen mit verschiedenen Seeds).
	•	Wir messen daraus strukturelle Veränderungen in der Membran unter dem Peptid: z.B. lokale Bilayer-Dicke (Abstand zwischen äußeren und inneren Leaflet-Phosphaten) unter dem Peptid vs. in einem peptidfreien Referenzbereich; ferner die Wassereindringtiefe (wie viel Wasser dringt bis zur Membranmitte vor unter dem Peptid) und lokale Schwanz-Ordnung (z.B. orientierungsabhängiger Ordnungsparameter der Lipidalkylketten) im Einflussbereich des Peptids.
	•	Aus diesen Metriken ließe sich ein summarischer “Membran-Störungsindex” bilden. Beispielsweise könnten wir für jedes Peptid das mittlere prozentuale Dünnerwerden der Membran unter dem Peptid angeben. Vorarbeiten zeigen, dass manche amphipathische Peptide Membranen merklich ausdünnen oder „weicher“ machen ￼, was ein möglicher Mechanismus der Hämolyse ist.
	•	Diese Kennzahl(en) könnten wir dann als zusätzlichen Prädiktor ins Modell aufnehmen und testen, ob sich die Vorhersage verbessert. Sollte etwa ein Peptid keine tiefe Einlagerung (ungünstiges ΔGinsert) zeigen, aber dennoch stark membran-perturbierend wirken (hoher Störungsindex), könnte dies erklären, warum es hämolytisch ist – das Modell würde davon profitieren, diesen Aspekt zu kennen.

Wichtig: Dieser Schritt wird nur unternommen, wenn die alleinigen ΔG-Features nicht genügen, um die  Hypothesen zu stützen. Wir möchten ein möglichst einfaches Modell behalten, um Rechenaufwand zu sparen und Interpretierbarkeit zu wahren. Die Einführung einer weiteren komplexen Simulation pro Peptid (~1,0 µs zusätzlich) wäre nur gerechtfertigt, wenn klarer Nutzen zu erwarten ist. Falls wir darauf zurückgreifen, wird dies transparent im Bericht gemacht und die Entscheidung begründet.

Abschätzung des Rechenaufwands und Durchführbarkeit

Das Basisszenario umfasst pro Peptid ~0.64 µs MD-Sampling (Umbrellas; 16 Fenster × 20 ns × 2 Replikate). Je nach adaptiver Verdichtung im Bereich −0.6…+0.6 nm können es 18–22 Fenster werden (dann ~0.72–0.88 µs/Peptid). Für 42 Peptide ergibt das ~26.9 µs (Basis, 16 Fenster) bis ~37.0 µs (bei 22 Fenstern). Im Worst-Case (alle Fenster auf 30 ns verlängern) wären es ~0.96 µs pro Peptid, also ~40.3 µs gesamt. Diese Zahlen sind für moderne HPC-Cluster gut handhabbar; die Umbrellas sind „embarrassingly parallel“ (je Peptid/Replikat/Window unabhängige Jobs).

Alle Simulations-Skripte und Parameterfiles sind vorbereit und in einem Versionskontrollsystem (Git) verwaltet. Deterministische Seeds gewährleisten, dass bei Bedarf einzelne Läufe exakt reproduziert werden können. Datenmanagement: Ergebnisse pro Peptid werden in einem strukturierten Ordner abgelegt (results/{peptide_id}/...), inkl. Input-Files, Logs, Trajektorien und Auswerte-Skripten. Ein zentrales Manifest dokumentiert die verwendete Software-Version (z.B. GROMACS 202X, Python-Pakete mit Versionsnummern) und den Git-Commit der Analyseskripte, um die Reproduzierbarkeit sicherzustellen. Somit ist das Vorgehen nicht nur ressourcenmäßig geplant, sondern auch organisatorisch.

Geplante Ergebnisdarstellung

Wir werden die Resultate klar und übersichtlich aufbereiten, um die Hypothesenprüfung zu illustrieren:
	•	PMF-Plots pro Peptid: Für jedes Peptid erstellen wir einen Graphen G(z) vs. z. Darin markieren wir die identifizierten Minima (zsurf, zhead) und ggf. die Barriere. Ein schattierter Bereich zeigt das 95%-Konfidenzband (aus dem Bootstrap) an. Diese Einzelplots wandern in den Anhang oder die Suppl. Info, und exemplarisch werden 1–2 im Haupttext gezeigt, um typische Unterschiede zwischen toxisch vs. nicht-toxisch zu veranschaulichen.
	•	Tabelle der ΔG-Features: Eine kompakte Tabelle (CSV und im Manuskript) listet für jedes Peptid ΔGads, ΔGinsert, ΔG‡ und deren geschätzte Fehler, sowie die Positionen der Minima. Dazu den experimentellen HC50-Wert. Diese Tabelle erlaubt es, die Rohdaten der Analyse nachvollziehen.
	•	Korrelation Scatterplots: Wir plotten -ΔGinsert und -ΔGads jeweils gegen log10(HC50). Auf diesen Streudiagrammen ziehen wir eine Regressionslinie (robuste lineare Regression oder Lowess) ein, um den Trend zu visualisieren. Fehlende oder ausgeschlossene Peptide werden gekennzeichnet. Die Korrelationskoeffizienten und p-Werte werden im Plot oder der Legende angegeben.
	•	ROC- und PR-Kurven: Die Leistung des Klassifikators wird mit einer ROC-Kurve dargestellt (True Positive Rate vs. False Positive Rate). Idealerweise liegt unsere Kurve deutlich über der Diagonale; wir geben AUC und CI an. Eine PR-Kurve (Precision vs. Recall) ergänzt dies, vor allem um zu zeigen, wie präzise ein hoher Score tatsächlich Toxizität vorhersagt.
	•	Kalibrierungsplot: Wir unterteilen die vorhergesagten Wahrscheinlichkeiten in Bins (z.B. 0–0,1, 0,1–0,3, … 0,9–1) und stellen gegenüber, welcher Anteil der Peptide in jedem Bin tatsächlich toxisch war. Ein gut kalibriertes Modell würde nahe der Diagonalen liegen (z.B. von den Peptiden, denen ~80% Toxizitäts-Wahrscheinlichkeit zugeschrieben wurde, sollten ca. 80% tatsächlich toxisch sein). Dieser Plot zeigt, ob unser Modell zur Über-/Untertreibung neigt.
	•	Konfusionsmatrix: Für einen gewählten Schwellenwert (etwa den Youden-Cutoff oder die heuristische Regel) zeigen wir in einer 2×2-Matrix, wie viele Peptide korrekt als toxisch/nicht-toxisch erkannt wurden und wie viele Fehlklassifikationen (False Positives/Negatives) es gab. Dadurch werden Sensitivität und Spezifität anschaulich. Wahrscheinlich erstellen wir zwei solcher Matrizen – eine für den optimierten Schwellenwert, eine für die Heuristik – um deren Performance nebeneinander zu stellen.
	•	Sensitivitäts-Ergebnisse: Z.B. ein Scatterplot, der für die 10 Peptide im Cholesterin-Test -ΔGinsert bei 45% vs. 35% Chol vergleicht (Erwartung: Punkte ~ auf Diagonale). Oder eine Tabelle mit Kendall-τ Werten. Ebenso ein kurzer Bericht ob die Rangfolge mit/ohne EN gleich blieb (kann man textlich beschreiben, ggf. Mini-Plot). Diese Details kommen vermutlich in den Anhang oder methodischen Teil.
	•	QC-Übersicht (Supplement): Um zu demonstrieren, dass alle Daten solide sind, generieren wir einen aggregierten QC-Report. Dieser enthält z.B. eine Heatmap aller Histogram-Overlaps (Fenster i vs. i+1 für jedes Peptid als Matrix), ein Balkendiagramm der effektiven Sample Sizes pro Fenster, sowie exemplarisch je Peptid einen Vergleich der PMF-Hälften und Replikate (z.B. Replikat1/2 Kurven übereinander). Dieser Report zeigt auf einen Blick, dass keine größeren Unstimmigkeiten vorhanden sind und erhöht die Glaubwürdigkeit der Studie.

Insgesamt legen wir Wert auf eine einfache, prägnante Darstellung. Wichtige Ergebnisse (Korrelation, ROC) werden direkt im Haupttext präsentiert, während detailreiche Unterlagen (alle 42 PMFs, QC-Metriken) in den Anhang kommen. Durch diese Visualisierungen wird die Aussagekraft unseres Ansatzes greifbar.

Interpretation und Entscheidungskriterien

Abschließend bewerten wir, ob die Hypothesen gestützt wurden und was dies bedeutet:
	•	Validierung H1 (Ranking): Sollte sich z.B. herausstellen, dass Spearman ρ ~0,6 (p < 0,01) für -ΔGinsert vs. log HC50, dann wäre das ein deutliches Signal, dass stärker membranaffine Peptide tendenziell toxischer sind. Das würde die Grundannahme bestätigen, dass die thermodynamische Bindungsstärke an die Membran ein Schlüsselfaktor für Hämolyse ist – im Einklang mit früheren Befunden ￼. Wir könnten dann argumentieren, dass ein Simulationsscreening auf Membran-PMFs eine praktikable Methode zur Vorauswahl nicht-hämolytischer Kandidaten ist. Im quantitativen Sinne ließe sich sagen: eine Abnahme von ΔGinsert um z.B. 5 kJ/mol entspricht grob einer Verringerung des HC50 um soundso viel (wenn lineare Regressionsmodelle passt, könnte man sogar eine grobe Formel angeben).
	•	Validierung H1 (Klassifikation): Wenn die AUC ~0,85 ist und z.B. bei einem bestimmten Schwellenwert 18/20 toxische korrekt erkannt werden (Sens 90%) und 15/20 nicht-toxische korrekt ausgeschlossen (Spez 75%), dann hätten wir ein nützliches Werkzeug. Das würde bedeuten, wir können mit wenigen berechneten Parametern ziemlich zuverlässig vorhersagen, ob ein neues Peptid problematisch für RBCs wäre. Insbesondere, wenn die Heuristik (fester Cutoff) ähnlich gut funktioniert (angenommen, die o.g. Regel hätte 15/20 toxische und 14/20 nicht-toxische richtig – entspräche 75%/70%), könnte man diese als einfaches Entscheidungskriterium vorschlagen: “Falls ΔGads < -8 kJ/mol und ΔGinsert < -3 kJ/mol, sollte man das Peptid als potentiell hämolytisch betrachten und genauer testen.” Das wäre für Praktiker attraktiv, da es schnell aus einem PMF abgelesen werden kann.
	•	Grenzfälle: Sollte das Ergebnis grenzwertig ausfallen, etwa ρ ~0,4 oder AUC ~0,75, würde das Bild gemischter sein. Wir hätten dann zwar einen Trend, aber keine ausreichende Verlässlichkeit. In diesem Szenario würden wir schauen, ob der Monomer-Thinning-Index (falls berechnet) die Lücken schließt. Möglicherweise entdecken wir, dass einige Ausreißer-Peptide sich durch Membranstörung auszeichnen, obwohl ihr ΔGinsert nicht stark negativ war – dann könnte das erweiterte Modell besser performen. Wir würden berichten, dass mit dem zusätzlichen Feature die AUC z.B. auf 0,85 stieg, was andeutet, dass man für vollständige Vorhersage auch Membranperturbations-Effekte berücksichtigen muss.
	•	Misserfolg: Falls wir weder eine signifikante Korrelation noch eine brauchbare Trennung erreichen (ρ ~0 oder AUC ~0,5–0,6), wäre unsere Hypothese widerlegt: Dann scheint die Einzel-Peptid-Adsorption als Surrogatmaß nicht auszureichen, um Hämolyse abzuleiten. Gründe könnten sein: Hämolyse erfordert möglicherweise Peptid-Oligomerisierung oder Porenbildung, die im Monomer-PMF nicht abgebildet wird; oder andere Faktoren wie Peptid-Sequenzmotive, proteolytische Stabilität etc. spielen mit hinein. In diesem Fall würden wir offen legen, dass das vorgestellte Verfahren in seiner jetzigen Form die Erwartungen nicht erfüllt, und vorschlagen, in Zukunft komplexere Simulationen (z.B. explizite Oligomer-Porenbildungstests) oder Machine-Learning-Modelle mit zusätzlichen Deskriptoren zu nutzen. Selbst ein negatives Ergebnis wäre wissenschaftlich wertvoll, da es die Grenzen der Methode aufzeigt.

Abschließend fassen wir die Studie in einem Fazit zusammen und markieren die vorab definierten Erfolgskriterien als erfüllt oder nicht:
	•	H1 (Ranking) bestätigt? – Ja/Nein (je nachdem, ob ρ und CI den Kriterien genügen).
	•	H1 (Klassifikation) bestätigt? – Ja/Nein (ob AUC-Kriterium erfüllt).
	•	Heuristik ausreichend? – Ja/Nein (ob sie die anvisierten 75/70 erreicht hat).
	•	Robustheit gegeben? – Ja/Nein (z.B. ob Cholesterol-Änderung keinen großen Einfluss hatte, etc.).

Diese Punkte werden im Bericht explizit angesprochen, um den Erfolg des Ansatzes messbar zu machen.

Durch dieses Vorgehen stellen wir sicher, dass das Konzept sowohl wissenschaftlich fundiert als auch möglichst schlank und effizient ist. Jeder Schritt wurde daraufhin optimiert, das nötige Maß an Daten zu liefern, ohne unnötige Komplexität einzuführen. Damit hoffen wir, einen praxistauglichen Workflow zu etablieren, der künftig bei der Entwicklung neuer Peptide helfen kann, frühzeitig toxische Kandidaten auszusortieren – und das mit deutlich geringerem Aufwand als umfangreiche In-vitro-Tests.

⸻


Standardisierung (SOP) – Umsetzung in der Pipeline

• Leaflet-Fixierung: Statische `index_leaflets.ndx` (OuterPO4/InnerPO4); keine Online-Heuristik.
• z-Referenz (lokal): Patch-Referenz (Radius 1.8–2.0 nm) um Peptid-COM; alternativ PLUMED-zylinder.
• Membran-Drift: optional `-DPOSRES_LIPID_Z` (z-only, sehr schwach); 0.5–1 ns Re-Equil bei Aktivierung.
• Replikate: zwei Orientierungen; deterministische Seeds aus Peptid-ID.
• Pfad/Start: SMD-Leiter → 1 ns Pre-Equil/Fenster → 20 ns Produktion (bis 30 ns bei Bedarf).
• Fenster: +2.8…−2.0 nm; 0.15 nm (dicht), 0.2 nm (außen); adaptive Verdichtung auf Overlap O ≥ 0.20.
• Output: Pull-Coordinate/Force ≤1 ps; QC-Gates automatisch enforced (Overlap, ESS, Halbzeit, Replikat).
• Analyse: WHAM primär; MBAR optional/QA; PCHIP (0.01 nm) nur zur Extrempunkt‑Bestimmung; präzisionsgewichtete Replikat‑Mittel.
• Reporting: QC-Report (Overlap-Heatmap, ESS, Halbzeit, R1/R2), PMF mit 95%-Band, Feature-Tabelle.

python3 scripts/universal/02_run_colabfold.py

python3 scripts/universal/02_run_colabfold.py simulations/solvia_8_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_14_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_32_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_68_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_126_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_164_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_215_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_482_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_490_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_515_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_524_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_527_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_617_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_624_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_850_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_858_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_941_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_974_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_1023_run_1
python3 scripts/universal/02_run_colabfold.py simulations/solvia_1045_run_1