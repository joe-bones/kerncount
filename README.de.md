# kerncount: automatisierte Zellkernzählung und Auswertung

*[English version](README.md)*

`kerncount` ersetzt einen manuellen Arbeitsablauf aus Fiji-Makros und Prism für Fluoreszenz-Zellzählexperimente
durch ein einziges Python-Werkzeug. Es leistet Folgendes:

1. Zählt die Zellkerne in jedem Bild (DAPI) mit einem Deep-Learning-Detektor.
2. Schätzt ab, wie verlässlich jede Zählung ist, und erstellt für jedes Bildfeld ein Kontrollbild.
3. Ermöglicht es, Zählungen in einer einzigen CSV-Datei von Hand zu korrigieren. Jede Korrekturrunde wird dabei als eigene „Generation“ gespeichert.
4. Erstellt Tabellen, Prism-fertige CSV-Dateien, Diagramme und Statistiken für jede Behandlung, jedes n und alle n zusammen.
5. Misst die Intensität pro Zelle (in jedem beliebigen Kanal), so wie bisher die Mean-intensity-Makros.
6. Schätzt ab, wie viele Bilder pro Bedingung tatsächlich benötigt werden.

```
kerncount/                 (dieses Repository)
├── README.md              ← diese Anleitung (Englisch); README.de.md = deutsche Fassung
├── docs/technical.md      ← Algorithmen und Implementierungshinweise
├── src/kerncount/         ← der Code
├── tests/                 ← Unit-Tests (pytest)
└── setup.bat / setup.command / kerncount.bat / kerncount-gui.*   ← Hilfsdateien für Installation und Start
```
Ihre Bilder und Ergebnisse liegen in einem separaten **Arbeitsbereich** (Workspace), in dem Eingabe und Ausgabe getrennt sind:
```
my_workspace/                    ← alle Befehle aus diesem Ordner ausführen
├── data/                        ← EINGABE: ein Ordner pro Experiment (Rohbilder; werden nie verändert)
│   ├── my_experiment/
│   └── another_experiment/
└── analysis/                    ← AUSGABE: ein Projekt pro Experiment, gleicher Name wie in data/
    ├── my_experiment/
    └── another_experiment/
```
In `data/` wird nie etwas verändert. Alles, was für ein Experiment erzeugt wird, liegt in `analysis/<experiment>/`.

---

## 1. Eingabe: Wie die Rohdaten organisiert sein müssen

Jedes Experiment ist ein Ordner in `data/`:

```
data/<experiment>/                              z. B. data/my_experiment/
├── <ein Ordner pro n (Versuch)>/               z. B. "05.06.26 DAPI count"
│   ├── <ein Ordner pro Behandlung>/            z. B. "hsamiRNA 21-5p", "LyoVec veh. cntrl", "cntrl"
│   │   ├── 1.tif
│   │   ├── 2.tif
│   │   └── ...
│   └── ...
└── ...
```

- **Ein Ordner pro n**, mit dem Datum im Namen (`05.06.26`, gelesen als Tag.Monat.Jahr). Die n werden
  nach Datum nummeriert: `n1_2026-06-05`, `n2_2026-07-06` usw. Ein einzelnes n funktioniert ebenfalls.
- **Ein Ordner pro Behandlung** innerhalb jedes n. Uneinheitliche Namen werden automatisch bereinigt; die
  Zuordnung können Sie in `layout.csv` überprüfen:
  `hsamiRNA 21-5p` / `hsa miRNA 21-5p` / `21 5p` / `21` → `miR-21-5p`,
  `LyoVec veh. cntrl` / `lyo vec` → `LyoVec`, `cntrl` / `010125 cntrl` / `-PDL cntrl` → `Control`. (`LyoVec`, `CL075` und miRNA-Muster sind eingebaute Regeln; alles
  Unbekannte wird unverändert übernommen und markiert.)
- **Bilder sind nach Nummer benannt** (`1.tif`, `02.tif`; zusätzliche Zeichen wie in `` `10.tif `` sind kein Problem).
  Dateien ohne Nummer (z. B. `5 zoom.tif`) werden übersprungen und in `layout.csv` aufgeführt.
  Die Anzahl der Bilder pro Ordner ist beliebig.
- **Ignoriert werden:** Unterordner namens `output`, Dateien, deren Name `mask` enthält (alte Fiji-Ergebnisse), sowie
  zusätzliche Unterordnerebenen (Bilder werden in jeder Tiefe gefunden).
- **Bildformat:** TIFF, 8 oder 16 Bit, Graustufen, RGB oder mehrkanalig, bei jeder Vergrößerung und von jeder Kamera.
  Die Kerngröße wird pro Ordner gemessen, und der Detektor skaliert entsprechend.
- **Platten:** Standardmäßig gilt Bild 1–6 = Platte 1 und Bild 7–12 = Platte 2 (wird nur für die Sampling-Analyse
  verwendet; änderbar in `settings.toml`).

---

## 2. Installation (einmal pro Computer)

1. **Mac/Linux:** Installieren Sie **Python 3.12** von python.org (ab 3.11 funktioniert es; 3.12 wird für TensorFlow empfohlen).
   **Windows:** Keine Vorinstallation und keine Administratorrechte nötig: `setup.bat` holt Python 3.12 selbst
   (über [uv](https://docs.astral.sh/uv/)) in Ihren Benutzerordner (`%LOCALAPPDATA%\kc`).
2. Laden Sie dieses Repository herunter (grüne Schaltfläche **Code** → *Download ZIP*, oder `git clone`).
3. Doppelklicken Sie auf `setup.bat` (Windows) oder führen Sie `setup.command` aus (Mac/Linux).
   Das dauert einige Minuten und benötigt etwa 1 GB. Das StarDist-Modell (ca. 5 MB) wird beim ersten Lauf
   automatisch heruntergeladen; dafür ist eine Internetverbindung nötig.

Befehle werden aus Ihrem Arbeitsbereich-Ordner ausgeführt. Unter Windows beginnt jeder der folgenden Befehle statt mit
`kerncount` mit dem vollständigen Pfad zu `kerncount.bat` in diesem Ordner (z. B. `C:\tools\kerncount\kerncount.bat run my_experiment`).
Unter Mac/Linux lautet er `<dieser Ordner>/.venv/bin/kerncount`. Wenn Sie diesen Ordner zu Ihrem PATH hinzufügen, genügt `kerncount`.

---

## 3. Bedienung

Führen Sie jeden Befehl aus dem Arbeitsbereich-Ordner aus (dem Ordner, der `data/` enthält). Experimente werden
über ihren Ordnernamen in `data/` angesprochen.

### Ohne Befehle: die Weboberfläche
Doppelklicken Sie auf **`kerncount-gui.bat`** (Windows) bzw. **`kerncount-gui.command`** (Mac). Im Browser öffnet
sich eine Seite mit einer Zeile pro Experiment und Knöpfen für jeden der folgenden Schritte (**1 Set up**,
**2 Run**, **Sampling**) sowie Knöpfen, die Layout, Einstellungen, Kontrollbilder, `counts.csv` und die
Excel-Ergebnisse öffnen. Das Protokoll zeigt den Fortschritt live; **Stop** bricht einen Lauf ab. Lassen Sie das
schwarze Fenster während der Arbeit offen und schließen Sie es zum Beenden. Die Seite funktioniert nur auf diesem
Computer, und es wird nichts nach außen gesendet. Die Einstellungen unter *Options* entsprechen den
Befehlsoptionen am Ende dieses Abschnitts. Dasselbe aus einem Terminal: `kerncount gui`.

### Schritt 1: Projekt anlegen (einmal pro Experiment)
```
kerncount init my_experiment
```
Dies durchsucht `data/my_experiment/` und legt `analysis/my_experiment/` mit zwei Dateien an
(`kerncount init --all` tut dies für jedes Experiment in `data/`):

- **`layout.csv`** enthält eine Zeile pro Rohdatenordner. Prüfen Sie die Spalten `n` und `treatment`. Um einen Namen zu korrigieren,
  bearbeiten Sie ihn; um einen Ordner auszuschließen, setzen Sie `include` auf `no`. Die Spalte `note` zeigt Auffälligkeiten
  (übersprungene Dateien, abgeleitete Namen, „72h“).
- **`settings.toml`** enthält alle Einstellungen (Abschnitt 4). Die Standardwerte sind sinnvoll gewählt.

### Schritt 2: Analyse ausführen
```
kerncount run my_experiment
```
(`kerncount run --all` führt alle initialisierten Experimente aus.)
Der erste Lauf dauert etwa 2–3 s pro Bild (rund 15 min für 370 Bilder). Danach sind die Ergebnisse
zwischengespeichert (Cache), und ein erneuter Lauf dauert unter einer Minute. Ein erneuter Lauf ist **jederzeit gefahrlos möglich**.

### Schritt 3: Zählungen überprüfen
Öffnen Sie `analysis/my_experiment/review/<n>/<treatment>_overview.jpg`. Es gibt eine Übersicht pro Behandlung, die jedes Bildfeld
mit seiner Zählung und einem Rahmen in der Farbe der Qualitätsstufe zeigt: **grün = hoch, orange = mittel, rot = niedrig**.
Für eine genauere Ansicht öffnen Sie das Einzelbild `review/<n>/<treatment>/<NN>.jpg`
(links: Original; rechts: was gezählt wurde). Beginnen Sie mit den roten.

### Schritt 4: Zählungen von Hand korrigieren (optional)
Öffnen Sie **`analysis/my_experiment/counts.csv`** in Excel:

| n | treatment | image | auto_count | qc_ok | qc_level | manual_1 | manual_2 | manual_3 |
|---|---|---|---|---|---|---|---|---|
| n1_2026-06-05 | Control | 3 | 85 | yes | high | | | |
| n1_2026-06-05 | miR-21-5p | 6 | 140 | no | low | **121** | | |
| n1_2026-06-05 | miR-21-5p | 9 | 4 | no | low | **x** | | |

- Tragen Sie in `manual_1` eine korrigierte Zahl ein, oder **`x`**, um das Bild auszuschließen. Lassen Sie alle anderen Zellen leer.
- Speichern Sie und führen Sie `kerncount run my_experiment` erneut aus. Dadurch entsteht Generation **g2** (siehe Abschnitt 5).
- Eine spätere Korrekturrunde kommt in `manual_2`, woraus g3 entsteht, und so weiter.

### Schritt 5: Wie viele Bilder brauchen wir? (optional)
```
kerncount sampling my_experiment
```
Die Ergebnisse werden nach `analysis/my_experiment/sampling/g<N>/` geschrieben (Abschnitt 5.5).

### Einen zweiten Kanal messen (z. B. einen Marker in C2)
Für mehrkanalige TIFFs (C1 = DAPI, C2 = Signal) in `settings.toml` einstellen:
```
[channels]
names = ["DAPI", "marker"]
segment = 1            # Kerne im DAPI-Kanal finden

[intensity]
measure = [2]          # den Marker messen
regions = ["nucleus", "ring", "cell"]
```
Dann wie gewohnt ausführen. Werte pro Zelle landen in `_work/cells.csv` (eine Zeile pro erkanntem Objekt, mit dessen
Kategorie). Mittelwerte pro Bild landen in `qc_details.csv` und durchlaufen dieselben Statistiken pro n und über alle n wie die
Zählungen (Dateien `intensity_marker_cell_*`). Um statt Kernen ganze Zellen aus einem mCherry-Kanal zu segmentieren,
setzen Sie `segment` auf diesen Kanal und `[segmentation] backend = "classical"`.

### Befehlsoptionen
| Befehl | Funktion |
|---|---|
| `kerncount init <experiment>` | Legt `analysis/<experiment>/` aus `data/<experiment>/` an |
| `kerncount init --all` | Dasselbe für jedes Experiment in `data/` (bestehende Projekte bleiben erhalten) |
| `kerncount init <experiment> --force` | Neu durchsuchen und `layout.csv` sowie `settings.toml` überschreiben |
| `kerncount run <experiment>` / `--all` | Zählen, prüfen und auswerten für ein / alle Experimente |
| `kerncount run <experiment> --rebuild g2` | Einen bestehenden Generationsordner überschreiben (geschieht standardmäßig nie). `--rebuild all` für alle. |
| `kerncount run <experiment> --resegment` | Cache ignorieren und jeden Kern neu erkennen |
| `kerncount run <experiment> --no-review` | Kontrollbilder nicht erstellen (schneller) |
| `kerncount sampling <experiment> --generation g2` | Sampling-Analyse für eine bestimmte Generation (Standard: die neueste) |

Für Ordner außerhalb eines Arbeitsbereichs funktionieren auch Pfade: `kerncount init <Rohdatenordner> <Projektordner>`,
`kerncount run <Projektordner>`. Mit der Umgebungsvariable `KERNCOUNT_WORKSPACE` können Sie Experimente
von überall aus über ihren Namen ausführen.

---

## 4. Einstellungen (`settings.toml`)

Mit jedem Texteditor bearbeitbar. Eine gelöschte Einstellung fällt auf ihren Standardwert zurück.

### `[project]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `raw_dir` | (von init gesetzt) | Wo die Rohbilder liegen. Relativ zum Projektordner gespeichert (`../../data/<experiment>`), sodass der gesamte Arbeitsbereich verschoben oder kopiert werden kann. |
| `date_order` | `"DMY"` | Wie Datumsangaben in den Namen der n-Ordner geschrieben sind (`DMY`, `MDY`, `YMD`) |
| `csv_delimiter` | `","` | Verwenden Sie `";"`, wenn Ihr Excel (deutsche/europäische Ländereinstellung) CSV-Dateien in einer einzigen Spalte anzeigt. Beim Einlesen wird beides automatisch erkannt. |

### `[channels]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `names` | `[]` | Kanalnamen, z. B. `["DAPI", "marker"]`. Leer = C1, C2, ... |
| `segment` | `1` | Kanal, in dem die Kerne gesucht werden |

### `[segmentation]`: Kerne finden
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `backend` | `"stardist"` | `stardist` = Deep-Learning-Detektor (empfohlen); `classical` = Schwellenwert + Watershed wie in Fiji |
| `nucleus_diameter_px` | `"auto"` | `auto` misst die Kerngröße pro Ordner (passt sich der Vergrößerung an). Alternativ eine feste Zahl. |
| `target_diameter_px` | `30` | Größe, bei der der Detektor am besten arbeitet; Bilder werden intern darauf skaliert |
| `prob_thresh` | `"default"` | Konfidenzschwelle des Detektors (niedriger = mehr Objekte) |
| `pixel_size_um` | `""` | µm pro Pixel, falls bekannt (aus einem Objektmikrometer). Ergänzt µm² und Zellen/mm². Eine Änderung löst keine Neuerkennung aus. |

### `[counting]`: Was als Kern zählt
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `edge_rule` | `"half"` | `half` = unverzerrter Zählrahmen (vom oberen/linken Rand angeschnittene Kerne zählen, vom unteren/rechten nicht, sodass jeder Randkern im Mittel einmal zählt); `exclude` = alle Randkerne verwerfen; `include` = alle zählen |
| `min_area_frac` | `0.2` | Objekte kleiner als das 0,2-Fache eines typischen Kerns = Debris oder Mikronukleus (wird nicht gezählt) |
| `exclude_out_of_focus` | `true` | Kerne außerhalb der Fokusebene verwerfen (schwimmende oder tote Zellen) |
| `focus_frac` | `0.5` | Kantenschärfe unter dem 0,5-Fachen des typischen Werts ... |
| `texture_frac` | `0.8` | ... UND Chromatintextur unter dem 0,8-Fachen des typischen Werts = unscharf. Beide Werte müssen niedrig sein, damit helle Kerne in Überstrahlung nicht verloren gehen. |
| `merge_area_frac` | `2.0` | Objekte größer als das 2-Fache eines typischen Kerns und nicht konvex = „mögliche Verschmelzung“ (zählt als 1, wird markiert) |

### `[qc]`: Wie streng die Qualitätsstufen sind
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `high_max_uncertainty` | `0.10` | QC hoch erfordert, dass der plausible Bereich ≤ 10 % der Zählung beträgt ... |
| `high_max_disagreement` | `0.15` | ... und die zweite Methode auf 15 % genau übereinstimmt ... |
| `max_saturated_frac` | `0.02` | ... und ≤ 2 % der Kernpixel gesättigt sind |
| `low_min_uncertainty` | `0.25` | QC niedrig, wenn der Bereich ≥ 25 % der Zählung beträgt ... |
| `low_min_disagreement` | `0.35` | ... oder die beiden Methoden um ≥ 35 % abweichen |

Bei sehr kleinen Zählungen wird der Bereich relativ zu mindestens 5 berechnet, sodass ein Bild mit 2 Zellen
nicht automatisch „niedrig“ ist. Ein Bild mit 0 Kernen ist immer „niedrig“ (auf ein leeres oder unscharfes Feld prüfen).

### `[intensity]`: Intensität pro Zelle (ersetzt die Mean-intensity-Makros)
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `enabled` | `true` | Intensitäten messen |
| `measure` | `"all"` | Zu messende Kanäle, z. B. `[2]` nur für Kanal 2 |
| `regions` | `["nucleus", "cell"]` | `nucleus` (Kern); `ring` (Band um den Kern, d. h. Zytoplasma); `cell` (Kern + Ring, entspricht dem „enlarge 6 px“ des Fiji-Makros) |
| `ring_px` | `"auto"` | Breite des Bands; `auto` = 0,25 × Kerndurchmesser |
| `subtract_background` | `true` | Hintergrund jedes Bildes abziehen (Median der leeren Bereiche) |

### `[analysis]`: Statistik
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `reference` | `"LyoVec"` | Hauptvergleichsgruppe |
| `normalise_to` | `["LyoVec", "Control"]` | Gruppen, auf die Ergebnisse in % normiert werden (innerhalb jedes n) und gegen die getestet wird |
| `treatment_order` | `["Control", "LyoVec", "CL075"]` | Reihenfolge in Tabellen und Diagrammen; weitere Behandlungen folgen alphabetisch |
| `exclude_low_qc` | `false` | Bei `true` werden Bilder mit QC niedrig (und ohne manuellen Wert) aus der Statistik ausgeschlossen |

### `[design]`: verwendet von `kerncount sampling`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `plates` | `2` | Aufgenommene Platten (Wells) pro Bedingung und n |
| `fields_per_plate` | `6` | Bildfelder pro Platte |
| `plate_layout` | `"consecutive"` | `consecutive` = Bilder 1–6 Platte 1, 7–12 Platte 2; `alternating` = ungerade/gerade |

### `[review]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `overlay_alpha` | `0.35` | Transparenz der farbigen Überlagerung |
| `jpeg_quality` | `88` | Qualität der Kontrollbilder |

---

## 5. Ausgabe: Was herauskommt und wie man es liest

```
analysis/<experiment>/
├── settings.toml, layout.csv     Ihre Einstellungen und die Ordnerzuordnung
├── counts.csv                    ★ DIE Arbeitsdatei: eine Zeile pro Bild, alle n, manuelle Spalten
├── qc_details.csv                alle Messwerte pro Bild (siehe 5.2)
├── fiji_comparison.csv           alte Fiji-Zählung vs. neue Zählung, pro Bild
├── review/                       Bilder zur Kontrolle der Zählung (5.1)
├── generations/g1, g2, ...       Ergebnisse (5.3, 5.4)
├── sampling/g1/                  wie viele Bilder benötigt werden (5.5)
└── _work/                        Cache, cells.csv (Messwerte pro Zelle), backups/ von counts.csv
                                  (Löschen von _work erzwingt eine Neuerkennung)
```

### 5.1 Kontrollbilder: `review/`
```
review/<n>/<treatment>_overview.jpg     alle Felder einer Behandlung; Rahmen = QC-Stufe
review/<n>/<treatment>/<NN>.jpg         ein Feld: Original | Überlagerung
```
Die Kopfzeile zeigt die Zählung, den plausiblen Bereich, die QC-Stufe und die Gründe dafür. Farben der Überlagerung:

| Farbe | Bedeutung | Gezählt? |
|---|---|---|
| grün | Kern | ja |
| orange | mögliche Verschmelzung von 2+ sich berührenden Kernen | ja, als 1 (erhöht die Obergrenze des Bereichs) |
| magenta | Debris, Mikronukleus oder Fragment | nein |
| blau | unscharf (außerhalb der Fokusebene) | nein |
| gelb | vom unteren/rechten Bildrand angeschnitten | nein (Randregel) |

### 5.2 Werte pro Bild: `counts.csv` und `qc_details.csv`
- **`auto_count`**: die automatische Zählung.
- **`count_low` – `count_high`** (*plausibler Bereich*): was die Zählung ergäbe, wenn jede Grenzfallentscheidung
  (nahe der Größen- oder Fokusschwelle, geringe Detektorkonfidenz, mögliche Verschmelzungen) anders
  ausgefallen wäre. Ein enger Bereich bedeutet ein eindeutiges Bild.
- **`disagreement`**: Abweichung von einer unabhängigen klassischen Zählung (Schwellenwert + Watershed).
  Große Werte bedeuten ein schwieriges Bild (dicht, verklumpt, helles Zytoplasma).
- **`qc_level`** (high / medium / low) fasst das Obige zusammen. **Niedrig heißt „dieses Bild ansehen“, nicht
  „diese Zählung ist falsch“.** Bilder mit niedriger QC häufen sich in spärlich besiedelten, geschädigten Wells (z. B. miRNA-behandelt), wo
  das Zählen tatsächlich schwieriger ist.
- `n_debris`, `n_out_of_focus`, `n_edge`, `n_merged`: wie viele Objekte jeder Art ausgeschlossen oder markiert wurden.
- `diameter_px`: typische Kerngröße in diesem Ordner (nützlich als Kontrolle von Vergrößerung bzw. Biologie).
- `fiji_count`: was das alte Makro ergab, zum Vergleich.
- `intensity_<channel>_<region>`: mittlere, hintergrundkorrigierte Intensität der gezählten Zellen.
- `nuclei_per_mm2`: nur vorhanden, sobald `pixel_size_um` gesetzt ist.
- `stardist_scale`, `saturated_frac`, `source_file`: Skalierungsfaktor des Detektors, Anteil gesättigter
  Kernpixel und die Rohdatei, aus der die Zeile stammt.

*Zu den QC-Stufen:* Sie sind Faustregeln, bis sie mit Handzählungen (Ground Truth) abgeglichen wurden.
Eine Handzählung von etwa 15 Bildern unterschiedlicher Dichte würde sie in echte Fehlerabschätzungen verwandeln.

### 5.3 Generationen: `generations/g1`, `g2`, ...
- **g1** = vollautomatische Zählungen. Immer vorhanden.
- **g2** = g1 + die in `manual_1` eingetragenen Korrekturen. **g3** = g2 + `manual_2`, und so weiter.
- Die Kette endet bei der ersten manuellen Spalte, die vollständig leer ist.
- Ein bestehender Generationsordner wird **nie überschrieben**, außer Sie verlangen es (`--rebuild g2`). Passt counts.csv
  nicht mehr zu einer bestehenden Generation (Sie haben `manual_1` bearbeitet, nachdem g2 erstellt wurde), gibt `run` eine Warnung aus.
- Jede Generation enthält `README.txt` (woraus sie erstellt wurde), `counts_used.csv` (der für jedes Bild verwendete Wert
  und seine Herkunft: `auto`, `manual_1`, `manual_1: excluded`, ...) und `generation.json`.

### 5.4 Innerhalb einer Generation: Behandlung → n → alle n
```
g1/
├── results.xlsx                      ★ alles Folgende in einer übersichtlichen Arbeitsmappe
├── counts_used.csv
├── <n>/                              ein Versuch
│   ├── <treatment>/fields.csv        Werte pro Feld (Zählung + Intensitäten)
│   ├── <treatment>/count_fields.png  Balken pro Feld (x = ausgeschlossen)
│   ├── count_summary.csv             pro Behandlung: Felder, Mittelwert, SD, SEM, Median, % von LyoVec, % von Control
│   ├── count_prism.csv               Prism-Column-Tabelle: Spalten = Behandlungen, Zeilen = Felder
│   ├── count_stats.csv               einfaktorielle ANOVA + Dunnett gegen die Referenz
│   └── count.png / .pdf              Balken = Mittelwert, Fehlerbalken = SD, Punkte = Felder, Sterne = Dunnett (5.4.2)
└── all_n/                            alle Versuche zusammen
    ├── count_summary.csv             pro Behandlung: Wert pro n, N, Mittelwert, SD, SEM, % der Kontrollen
    ├── count_prism.csv               Prism-Column-Tabelle: Spalten = Behandlungen, Zeilen = n
    ├── count_prism_pct_of_LyoVec.csv dasselbe, innerhalb jedes n auf LyoVec normiert (auch _Control)
    ├── count_stats.csv               RM-einfaktorielle ANOVA + gepaarte Ratio-t-Tests gegen jede Kontrolle
    ├── count.png / .pdf              Balken = Mittelwert der n, Fehlerbalken = SEM, Punktfarbe = n (Legende), Sterne (5.4.2)
    ├── count_pct_of_LyoVec.png/.pdf  normierte Version (gestrichelte Linie = 100 %); auch _Control
    └── n_colours.csv
```
Derselbe Satz an Dateien existiert für jede Intensitätsmessgröße (`intensity_<channel>_<region>_*`).

**Prism:** Öffnen Sie eine `*_prism*.csv` und fügen Sie sie in eine Column-Tabelle ein, oder verwenden Sie File → Import.
In den All-n-Dateien steht die n-Bezeichnung als Zeilentitel in der ersten Spalte.

**Tabellenblätter in `results.xlsx`:** *About* (was diese Generation ist) · *Counts used* · *Count - all n* ·
*Count - stats all n* · *Count - per n* · *Count - stats per n* · und, falls die Intensitätsmessung aktiv ist,
*Intensity - all n* · *Intensity - stats all n* · *Intensity - per n*.

#### 5.4.1 Spalten von `*_stats.csv`
| Spalte | Bedeutung |
|---|---|
| `test` | Welcher Test, und gegen welche Kontrolle verglichen wird |
| `comparison` | `<treatment> vs <control>`, oder `all` für eine Gesamt-ANOVA |
| `N` | (alle n) Anzahl der n, die beide Gruppen enthalten |
| `ratio`, `ci_low`, `ci_high` | (alle n) geometrisch gemittelte x-fache Veränderung gegenüber der Kontrolle und ihr 95-%-KI |
| `statistic` | t-, F- oder Dunnett-Statistik |
| `p` / `p_adj` | roher / für multiples Testen adjustierter p-Wert (Dunnett ist per Konstruktion adjustiert; t-Tests nach Holm) |
| `significance` | Sterne aus `p_adj`: ns, \*, \*\*, \*\*\*, \*\*\*\* (< 0,05 / 0,01 / 0,001 / 0,0001) |

#### 5.4.2 Wie Signifikanz in den Diagrammen dargestellt wird
- **Sterne** (\* bis \*\*\*\*) über einem Balken erscheinen nur, wenn der Vergleich dieser Behandlung signifikant ist
  (`p_adj` < 0,05). **Keine Markierung = nicht signifikant.** Die verwendete Kontrolle ist über ihrem Balken mit **„ref“** beschriftet.
- **Diagramme pro n** verwenden Dunnett gegen die Referenz. Fehlt die Referenz in diesem n, wird die nächste Kontrollgruppe
  (aus `normalise_to`) verwendet, und die Fußzeile nennt sie. Diese Sterne beschreiben nur einen einzigen Versuch.
- **All-n-Diagramme** verwenden den gepaarten Ratio-t-Test (Holm-adjustiert) gegen die gezeigte Kontrolle: Das Rohdiagramm verwendet
  `reference`; jedes `pct_of_<X>`-Diagramm verwendet X. Die Fußzeile nennt außerdem den p-Wert der Gesamt-RM-ANOVA.
  Eine signifikante ANOVA ohne Sterne an den Balken bedeutet, dass sich die Behandlungen als Gruppe unterscheiden, aber bei dieser
  Anzahl von n noch kein einzelner Vergleich signifikant ist.
- Die Fußzeile jedes Diagramms nennt den Test, sodass eine Abbildung auch ohne die Tabellen lesbar ist.

**So liest man die Statistik:**
- **Innerhalb eines n** (`<n>/count_stats.csv`): Die Felder sind *technische* Replikate, daher beschreiben diese p-Werte
  nur diesen einen Versuch. Berichten Sie sie nicht als Ergebnis. Dunnett vergleicht gegen `reference`, oder gegen
  die nächste Kontrollgruppe, falls die Referenz in diesem n fehlt.
- **Alle n** (`all_n/count_stats.csv`): Jedes n ist ein biologisches Replikat (der Mittelwert seiner Felder).
  - *RM-einfaktorielle ANOVA* (über die n, die jede Behandlung enthalten): Gibt es überhaupt einen Unterschied?
  - *Gepaarter Ratio-t-Test gegen LyoVec / gegen Control*: Für jede Behandlung das Verhältnis zur Kontrolle innerhalb
    jedes n, getestet über die n hinweg. `ratio` = geometrisch gemittelte x-fache Veränderung (0,30 = 70 % weniger Zellen), mit einem
    95-%-KI. `p_adj` ist innerhalb jeder Kontrollfamilie nach Holm korrigiert. Sterne: ns, \*, \*\*, \*\*\*, \*\*\*\*
    (0,05 / 0,01 / 0,001 / 0,0001).
  - Ein Test verwendet nur n, die *beide* Gruppen enthalten. Daher kann `N` zwischen Vergleichen variieren, wenn eine
    Gruppe in manchen n fehlt.
- **Normierte Werte** (`pct_of_...`) eignen sich im Allgemeinen besser für Vergleiche über die n hinweg, weil sie
  tagesabhängige Unterschiede bei Aussaat, Färbung und Lampenhelligkeit herausrechnen.
- **Intensitätswerte** sind nur vergleichbar, wenn Belichtung und Verstärkung (Gain) innerhalb jedes n fest eingestellt waren (Auto-Belichtung aus).

### 5.5 Wie viele Bilder werden benötigt? `sampling/g1/`
| Datei | Inhalt |
|---|---|
| `REPORT.txt` | Zusammenfassung in Klartext: hier beginnen |
| `images_needed.png/.pdf` | Vorhergesagtes Rauschen, Übereinstimmung beim Subsampling und Fehler pro n in Abhängigkeit von der Bildzahl pro Bedingung |
| `variance_components.csv` | Variation von Feld zu Feld (innerhalb einer Platte) und von Platte zu Platte (innerhalb eines n), als CV² |
| `effects.csv` | Beobachtete x-fache Veränderungen und biologische Varianz (von n zu n) pro Behandlung |
| `design_power.csv` | Für jedes Design (1–2 Platten × 1–6 Felder), pro Behandlung: Rauschen im Vergleich zu jetzt, Power bei Ihrer aktuellen Anzahl von n, und die für 80 % Power nötige Anzahl von n |
| `subsampling.csv` | Echte Bilder, pro Design vielfach neu gezogen: wie oft das Ergebnis mit der Schlussfolgerung aus allen Bildern übereinstimmt, und wie stark sich das Verhältnis jedes n verschiebt |

**So liest man es:** Der Test über die n wird durch zwei Dinge begrenzt. Das erste ist das *Messrauschen*
(durch Felder und Platten), das mit mehr Bildern kleiner wird. Das zweite ist die *biologische Variation von n zu n*,
die nur mit mehr n kleiner wird. Mehr Bilder helfen nur, solange das Messrauschen noch in einer ähnlichen Größenordnung wie die
biologische Variation liegt. Das Rauschen von Platte zu Platte wird nur mit mehr Platten kleiner.

**Typisches Muster:** Wenn Zellen verklumpen, ist die Variation von Feld zu Feld groß. Eine Aufnahme mit geringerer
Vergrößerung (10x deckt etwa die 16-fache Fläche von 40x ab) verringert das Messrauschen dann weit stärker als
zusätzliche Felder bei hoher Vergrößerung. Bei nur wenigen n sind die Ergebnisse als Planungshilfe zu verstehen, nicht als
exakte Werte.

---

## 6. Typische Arbeitsabläufe
**Ein neues n für ein bestehendes Experiment**
1. Legen Sie den neuen n-Ordner in `data/<experiment>/` ab (gleiche Struktur wie die anderen).
2. `kerncount init <experiment>`. Ihre `layout.csv` bleibt erhalten, und ein neuer Scan wird in `layout.new.csv` geschrieben;
   übernehmen Sie die neuen Zeilen in `layout.csv` (Namen prüfen).
3. `kerncount run <experiment>`. Nur die neuen Bilder werden verarbeitet; der Rest kommt aus dem Cache.
   Bestehende Generationen werden nicht überschrieben. Erstellen Sie sie gezielt mit `--rebuild` neu, sobald das neue n geprüft ist.

**Ein neues Experiment**
1. Legen Sie `data/<neues Experiment>/` mit einem Ordner pro n an.
2. `kerncount init <neues Experiment>`, `layout.csv` prüfen, dann `kerncount run <neues Experiment>`.

**Danach, in beiden Fällen:** Prüfen Sie die roten und orangen Kacheln in `review/*_overview.jpg`, korrigieren Sie Zählungen in
`counts.csv` → `manual_1`, und führen Sie den Lauf erneut aus. Berichten Sie aus der neuesten Generation: `generations/g<N>/results.xlsx`,
die `all_n/`-Diagramme und die Prism-CSV-Dateien.

## 7. Fehlerbehebung
| Problem | Lösung |
|---|---|
| „project not found“ / „raw folder not found“ | Führen Sie den Befehl aus dem Arbeitsbereich-Ordner aus (dem Ordner mit `data/`), prüfen Sie den Experimentnamen und führen Sie zuerst `kerncount init <experiment>` aus. |
| Der erste Lauf schlägt beim Herunterladen des Modells fehl | Das StarDist-Modell wird einmalig bei der ersten Verwendung heruntergeladen; mit dem Internet verbinden und erneut ausführen. |
| `counts.csv` öffnet sich in Excel als eine einzige Spalte | `csv_delimiter = ";"` in `settings.toml` setzen (europäisches Excel). Dateien, die Excel mit `;` neu speichert, werden in jedem Fall korrekt gelesen. |
| `run` kann `counts.csv` nicht schreiben (Windows) | Schließen Sie die Datei zuerst in Excel; Excel sperrt geöffnete Dateien. |
| „cannot read manual value“ | Eine manuelle Zelle enthält etwas anderes als eine ganze Zahl oder `x`; die Meldung nennt Zeile und Spalte. |
| „g2 already exists but counts.csv now gives different values“ | Sie haben eine manuelle Spalte geändert, nachdem g2 erstellt wurde. Führen Sie den Lauf mit `--rebuild g2` erneut aus, um sie gezielt zu aktualisieren. |
| Ein Behandlungsname ist falsch oder ein Ordner soll ausgeschlossen werden | `layout.csv` bearbeiten (`treatment`, `include`) und erneut ausführen. |
| „has 1 channel(s) but settings ask for channel 2“ | `[channels] segment` verweist auf einen Kanal, den die Bilder nicht haben. (Nicht vorhandene Kanäle in `[intensity] measure` werden stillschweigend übersprungen; prüfen Sie die Spalten `intensity_*`.) |
| Der Lauf ist langsam | Erste Läufe benötigen ca. 2–3 s pro Bild; spätere Läufe nutzen den Cache. `--no-review` überspringt die Kontrollbilder. |
| Zählungen sehen in einem ganzen Ordner falsch aus | `diameter_px` in `qc_details.csv` prüfen; ist die Größenschätzung daneben (z. B. bei einem fast leeren Ordner), `nucleus_diameter_px` auf eine Zahl setzen. |

## 8. Hinweise und bekannte Einschränkungen
- 8-Bit-Kamerabilder, die als 16-Bit-TIFF gespeichert wurden (jeder Wert × 257), werden korrekt verarbeitet. Sättigung wird in beiden Fällen erkannt.
- Im Vergleich zu einem einfachen Fiji-Schwellenwert-Makro fallen die Zählungen meist niedriger aus. Schwellenwert-Makros zählen
  Hintergrundrauschen in leeren oder unscharfen Feldern, teilen gelappte Kerne auf und zählen jeden Kern, der den Rand berührt.
- Sehr dichte Felder, in denen Zellen helles Zytoplasma haben, sind der schwierigste Fall. Sie werden als QC niedrig/mittel markiert;
  prüfen Sie sie in den Kontrollbildern.
- Ohne bekannte Pixelgröße werden Ergebnisse pro Feld statt pro mm² angegeben.
- Die QC-Stufen sind Heuristiken, bis sie anhand von Handzählungen kalibriert sind.
- Technische Details (Algorithmen, Tests): `docs/technical.md`.

## 9. Entwicklung
```
pip install -e ".[test]"
pytest
```
Hinweise zu Algorithmen und Implementierung: `docs/technical.md`.

## Danksagung
Die Kernerkennung verwendet [StarDist](https://github.com/stardist/stardist) (Schmidt et al., MICCAI 2018;
Weigert et al., WACV 2020) mit dem vortrainierten Modell `2D_versatile_fluo`.
