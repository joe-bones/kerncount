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
7. Misst mCherry-Tumoren in z-Stapeln von Hirnschnitten (Tumorexperimente, Abschnitt 9), mit verblindetem Review und vorab festgelegter Statistik.

```
kerncount/                 (dieses Repository)
├── README.md              ← diese Anleitung (Englisch); README.de.md = deutsche Fassung
├── docs/technical.md      ← Algorithmen und Implementierungshinweise
├── src/kerncount/         ← der Code
├── tests/                 ← Unit-Tests (pytest)
├── launcher/              ← Quellcode der Ein-Klick-Programme kerncount.exe / kerncount.app
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
- **Wells:** Ein n = eine Platte. Jede Bedingung liegt in 2 Wells dieser Platte; standardmäßig gilt Bild 1–6 = Well 1
  und Bild 7–12 = Well 2 (wird nur für die Sampling-Analyse verwendet; änderbar in `settings.toml`).

---

## 2. Installation (einmal pro Computer)

### Am einfachsten: Ein-Klick-Programm (ohne Administratorrechte, ohne Python, ohne Befehle)
Laden Sie aus dem [neuesten Release](https://github.com/joe-bones/kerncount/releases/latest) herunter:
- **Windows:** `kerncount.exe`
- **Mac (Apple Silicon):** `kerncount-mac.zip`. Ein Doppelklick entpackt es zu `kerncount.app`.

Legen Sie es in Ihren Arbeitsbereich-Ordner (neben `data`) und doppelklicken Sie darauf. Der erste Start braucht
Internet, einige Minuten und etwa 1 GB: Python 3.12 und alles, was kerncount benötigt, wird in Ihren Benutzerordner
installiert. Danach öffnet jeder Start kerncount sofort im Browser. Beim ersten Mal warnt das System vor einer
nicht signierten App: unter Windows **Weitere Informationen → Trotzdem ausführen**; am Mac **Fertig**, dann
**Systemeinstellungen → Datenschutz & Sicherheit → Trotzdem öffnen**. Details: [launcher/README.md](launcher/README.md).

### Alternativ: Installation aus diesem Repository (für die Kommandozeile)
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
Doppelklicken Sie auf **`kerncount.exe`** / **`kerncount.app`** (Abschnitt 2) oder, bei einer Repository-Installation,
auf **`kerncount-gui.bat`** (Windows) bzw. **`kerncount-gui.command`** (Mac). Dasselbe aus einem Terminal: `kerncount gui`.
Alles Folgende geht im Browser, ohne Dateien zu bearbeiten:

| Reiter | Was Sie dort tun |
|---|---|
| **Run** | Einrichten, Analyse starten (die Befehlsoptionen unten als Häkchen), Sampling; Live-Protokoll, **Stop** |
| **Layout** | n- / Behandlungsnamen korrigieren, Ordner ein- oder ausschließen |
| **Settings** | Alle Schalter aus Abschnitt 4 als Formular mit Erklärungen; oder `settings.toml` als Text bearbeiten |
| **Review images** | Übersicht pro Behandlung, dann jedes Bild in voller Größe (← / → zum Blättern); Korrektur direkt eintragen |
| **Corrections** | Eine Generation bearbeiten: korrigierte Zählungen oder **x** für eine Runde eintragen (`manual_1` → g2, …), Wert vorher und nachher sehen, dann die Generation **bauen** (etwa eine Minute; die Erkennung wird wiederverwendet) |
| **Results** | Alle Plots, Tabellen und Statistiken pro Generation; `results.xlsx` herunterladen |
| **Sampling** | Die Analyse „Wie viele Bilder?“ starten und lesen |
| **All files** | Alles im Projektordner durchsuchen und ansehen |

Tumorexperimente (Abschnitt 9) nutzen dieselben Seiten: Layout zeigt die Tumoren und was eine Analyse blockiert,
Settings den Analyseplan mit seiner Sperre, Review images und Corrections nur verblindete Codes, Results die Diagramme pro Messgröße.

Lassen Sie das schwarze (Windows) bzw. Terminal-Fenster (Mac) während der Arbeit offen; schließen Sie es oder
drücken Sie **Quit** zum Beenden. Die Seite funktioniert nur auf diesem Computer, und es wird nichts nach außen gesendet.

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
| `kerncount regen <experiment> [--rebuild g2 g3]` | Generationen nur aus `counts.csv` neu bauen (schnell; nach Korrekturen) |
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
| `reference` | `"Control"` | Hauptvergleichsgruppe (Hauptdiagramme) |
| `normalise_to` | `["Control", "LyoVec"]` | Gruppen, auf die Ergebnisse in % normiert werden (innerhalb jedes n) und gegen die getestet wird. Jede aufgeführte Gruppe erhält eigene Tests und Diagramme, Sie bekommen also immer *gegen Control* und *gegen LyoVec* |
| `two_way_anova` | `true` | Zusätzliche zweifaktorielle ANOVA mit Tukey-Vergleichen aller Paare: pro n Behandlung × Well, alle n Behandlung × n (5.4.3) |
| `treatment_order` | `["Control", "LyoVec", "CL075"]` | Reihenfolge in Tabellen und Diagrammen; weitere Behandlungen folgen alphabetisch |
| `exclude_low_qc` | `false` | Bei `true` werden Bilder mit QC niedrig (und ohne manuellen Wert) aus der Statistik ausgeschlossen |

### `[design]`: verwendet von `kerncount sampling`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `wells` | `2` | Aufgenommene Wells pro Bedingung und n (eine Platte pro n) |
| `fields_per_well` | `6` | Bildfelder pro Well |
| `well_layout` | `"consecutive"` | `consecutive` = Bilder 1–6 Well 1, 7–12 Well 2; `alternating` = ungerade/gerade |

Ältere Einstellungsdateien mit `plates` / `fields_per_plate` / `plate_layout` funktionieren weiterhin; der Bericht weist auf die Umbenennung hin.

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
│   ├── count_stats.csv               einfaktorielle ANOVA + Dunnett gegen jede Kontrolle (Control, LyoVec)
│   ├── count.png / .pdf              Balken = Mittelwert, Fehlerbalken = SD, Punkte = Felder, Sterne = Dunnett gegen Control (5.4.2)
│   ├── count_vs_LyoVec.png / .pdf    dasselbe, Sterne = Dunnett gegen LyoVec (nur wenn dieses n LyoVec enthält)
│   └── count_twoway.png / .pdf       zweifaktorielle ANOVA Behandlung × Well + Tukey alle Paare (5.4.3)
│       (+ count_twoway_anova.csv, count_twoway_pairs.csv)
└── all_n/                            alle Versuche zusammen
    ├── count_summary.csv             pro Behandlung: Wert pro n, N, Mittelwert, SD, SEM, % der Kontrollen
    ├── count_prism.csv               Prism-Column-Tabelle: Spalten = Behandlungen, Zeilen = n
    ├── count_prism_pct_of_LyoVec.csv dasselbe, innerhalb jedes n auf LyoVec normiert (auch _Control)
    ├── count_stats.csv               RM-einfaktorielle ANOVA + Dunnett gegen jede Kontrolle (Sterne); gepaarte Ratio-t-Tests (Kontrolle)
    ├── count.png / .pdf              Balken = Mittelwert der n, Fehlerbalken = SEM, Punktfarbe = n (Legende), Sterne (5.4.2)
    ├── count_pct_of_LyoVec.png/.pdf  normierte Version (gestrichelte Linie = 100 %); auch _Control
    ├── count_twoway.png / .pdf       zweifaktorielle ANOVA Behandlung × n + Tukey alle Paare (5.4.3)
    │   (+ count_twoway_anova.csv, count_twoway_pairs.csv)
    └── n_colours.csv
```
Derselbe Satz an Dateien existiert für jede Intensitätsmessgröße (`intensity_<channel>_<region>_*`).

**Prism:** Öffnen Sie eine `*_prism*.csv` und fügen Sie sie in eine Column-Tabelle ein, oder verwenden Sie File → Import.
In den All-n-Dateien steht die n-Bezeichnung als Zeilentitel in der ersten Spalte.

**Tabellenblätter in `results.xlsx`:** *About* (was diese Generation ist) · *Counts used* · *Count - all n* ·
*Count - stats all n* · *Count - per n* · *Count - stats per n* · und, falls die Intensitätsmessung aktiv ist,
*Intensity - all n* · *Intensity - stats all n* · *Intensity - per n* · und die zweifaktoriellen Ergebnisse:
*Two-way anova - all n* · *Two-way pairs - all n* · *Two-way anova - per n* · *Two-way pairs - per n*.

#### 5.4.1 Spalten von `*_stats.csv`
| Spalte | Bedeutung |
|---|---|
| `test` | Welcher Test, und gegen welche Kontrolle verglichen wird |
| `comparison` | `<treatment> vs <control>`, oder `all` für eine Gesamt-ANOVA |
| `N` | (alle n) Anzahl der n, die beide Gruppen enthalten |
| `effect_type`, `effect`, `ci_low`, `ci_high` | (alle n) Dunnett: Mittelwertdifferenz zur Kontrolle mit 95-%-KI (Dunnett-adjustiert); Ratio-t-Test: geometrisch gemittelte x-fache Veränderung mit 95-%-KI |
| `statistic` | t-, F- oder Dunnett-Statistik |
| `p` / `p_adj` | roher / für multiples Testen adjustierter p-Wert (Dunnett ist per Konstruktion adjustiert; t-Tests nach Holm) |
| `significance` | Sterne aus `p_adj`: ns, \*, \*\*, \*\*\*, \*\*\*\* (< 0,05 / 0,01 / 0,001 / 0,0001) |

#### 5.4.2 Wie Signifikanz in den Diagrammen dargestellt wird
- **Sterne** (\* bis \*\*\*\*) über einem Balken erscheinen nur, wenn der Vergleich dieser Behandlung signifikant ist
  (`p_adj` < 0,05). **Keine Markierung = nicht signifikant.** Die verwendete Kontrolle ist über ihrem Balken mit **„ref“** beschriftet.
- **Diagramme pro n** verwenden Dunnett gegen die Kontrolle: `count.png` gegen `reference` (Control), `count_vs_<X>.png`
  gegen jede weitere Kontrolle (LyoVec). Fehlt die Referenz in diesem n, verwendet `count.png` die nächste Kontrollgruppe,
  und die Fußzeile nennt sie. Diese Sterne beschreiben nur einen einzigen Versuch.
- **All-n-Diagramme** verwenden den Prism-Standard: RM-einfaktorielle ANOVA (n = gepaartes Replikat), danach Dunnett
  gegen die gezeigte Kontrolle mit dem gepoolten Fehler aller Behandlungen: Das Rohdiagramm verwendet
  `reference` (Control); jedes `pct_of_<X>`-Diagramm verwendet X, `pct_of_LyoVec` zeigt also die Sterne gegen LyoVec. Die Fußzeile nennt außerdem den p-Wert der Gesamt-RM-ANOVA.
  Eine signifikante ANOVA ohne Sterne an den Balken bedeutet, dass sich die Behandlungen als Gruppe unterscheiden, aber bei dieser
  Anzahl von n noch kein einzelner Vergleich signifikant ist.
- Die Fußzeile jedes Diagramms nennt den Test, sodass eine Abbildung auch ohne die Tabellen lesbar ist.

**So liest man die Statistik:**
- **Innerhalb eines n** (`<n>/count_stats.csv`): Die Felder sind *technische* Replikate, daher beschreiben diese p-Werte
  nur diesen einen Versuch. Berichten Sie sie nicht als Ergebnis. Dunnett wird gegen jede in diesem n
  vorhandene Kontrollgruppe (Control und LyoVec) gerechnet, jeweils als eigene Familie.
- **Alle n** (`all_n/count_stats.csv`): Jedes n ist ein biologisches Replikat (der Mittelwert seiner Felder).
  - *RM-einfaktorielle ANOVA*: Gibt es überhaupt einen Unterschied? Jedes n ist ein gepaartes Replikat. Die Anpassung
    erfolgt nach kleinsten Quadraten, sodass ein n, dem eine Gruppe fehlt (kein LyoVec in n2), für alle anderen Gruppen zählt.
    Ohne fehlende Werte entspricht das genau Prisms RM-einfaktorieller ANOVA; mit fehlenden Werten würde Prism sein
    Mixed-Effects-Modell verwenden, das sehr ähnliche Ergebnisse liefert.
  - *Dunnett gegen Control / gegen LyoVec* (**bestimmt die Sterne**): jede Behandlung gegen die Kontrolle, mit dem gepoolten
    Restfehler der ANOVA, adjustiert für die Anzahl der Vergleiche. `effect` = Mittelwertdifferenz mit 95-%-KI. Das Poolen
    ergibt viel mehr Freiheitsgrade als ein einzelner gepaarter Test (z. B. 20 statt 3 bei 4 n), sodass echte Effekte schon
    mit wenigen n signifikant werden. Annahme: Die Streuung zwischen den n ist für alle Behandlungen ähnlich.
  - *Gepaarter Ratio-t-Test gegen LyoVec / gegen Control* (**Kontrolle**, keine Sterne): Für jede Behandlung das Verhältnis zur Kontrolle innerhalb
    jedes n, getestet über die n hinweg. `ratio` = geometrisch gemittelte x-fache Veränderung (0,30 = 70 % weniger Zellen), mit einem
    95-%-KI. `p_adj` ist innerhalb jeder Kontrollfamilie nach Holm korrigiert. Sterne: ns, \*, \*\*, \*\*\*, \*\*\*\*
    (0,05 / 0,01 / 0,001 / 0,0001).
  - Ein Test verwendet nur n, die *beide* Gruppen enthalten. Daher kann `N` zwischen Vergleichen variieren, wenn eine
    Gruppe in manchen n fehlt.
- **Normierte Werte** (`pct_of_...`) eignen sich im Allgemeinen besser für Vergleiche über die n hinweg, weil sie
  tagesabhängige Unterschiede bei Aussaat, Färbung und Lampenhelligkeit herausrechnen.
- **Intensitätswerte** sind nur vergleichbar, wenn Belichtung und Verstärkung (Gain) innerhalb jedes n fest eingestellt waren (Auto-Belichtung aus).

#### 5.4.3 Zweifaktorielle ANOVA mit Vergleich aller Paare (`*_twoway*`)
Eine zusätzliche Auswertung auf beiden Ebenen, wie die zweifaktorielle ANOVA in Prism mit „compare all pairs“:
- **Pro n:** Behandlung × Well (Wells aus `[design]`); die Felder sind die Replikate.
- **Alle n:** Behandlung × n; die Felder sind die Replikate.

Das Modell ist `Wert ~ Behandlung + zweiter Faktor + Wechselwirkung` mit Typ-II-Quadratsummen, sodass eine Gruppe, die in
einem n fehlt (z. B. kein LyoVec in n2), korrekt behandelt wird. `*_twoway_anova.csv` enthält F, df, p und den Anteil an der
Gesamtvariation für Behandlung, zweiten Faktor und Wechselwirkung. `*_twoway_pairs.csv` vergleicht jedes Paar von Behandlungen
(Tukey–Kramer, `p_adj`) über die Zellmittelwerte innerhalb jedes Wells / n, nur über die Wells / n, die beide Gruppen enthalten
(`matched_levels`), mit Mittelwertdifferenz und 95-%-KI. Das Diagramm zeigt gruppierte Balken (Mittelwert ± SEM der Felder pro
Well / n) neben einer Matrix der adjustierten p-Werte.

**Vorsicht auf der Ebene aller n:** Weil die Felder als Replikate zählen, sind diese p-Werte optimistischer als die gepaarten
RM-ANOVA + Dunnett, die einen Wert pro n verwendet. Berichten Sie RM-ANOVA + Dunnett als Ergebnis; die zweifaktorielle Auswertung zeigt,
ob sich Behandlungen untereinander unterscheiden und ob der Effekt vom n abhängt (Wechselwirkung).

### 5.5 Wie viele Bilder werden benötigt? `sampling/g1/`
| Datei | Inhalt |
|---|---|
| `REPORT.txt` | Zusammenfassung in Klartext: hier beginnen |
| `images_needed.png/.pdf` | Vorhergesagtes Rauschen, Übereinstimmung beim Subsampling und Fehler pro n in Abhängigkeit von der Bildzahl pro Bedingung |
| `batch_by_n.csv` | Pro n: Werte der Referenzgruppen, fehlende Referenzen, systematischer Unterschied zwischen Wells und ob die Behandlungseffekte in diesem n stärker/schwächer sind als in den anderen |
| `batch_by_condition.csv` | Pro n × Behandlung: Mittelwert pro Well, Verhältnis Well zu Well mit Markierung, wenn es über das Feldrauschen hinausgeht, Verhältnis zu jeder Referenz |
| `variance_components.csv` | Variation von Feld zu Feld (innerhalb eines Wells) und von Well zu Well (innerhalb einer Platte), als CV² |
| `effects.csv` | Pro Referenz und Behandlung: beobachtete x-fache Veränderung, p mit allen Bildern, biologische Varianz (von n zu n) |
| `design_power.csv` | Für jedes Design (1–2 Wells × 1–6 Felder), pro Referenz und Behandlung: Rauschen im Vergleich zu jetzt, Power bei Ihrer aktuellen Anzahl von n, und die für 80 % Power nötige Anzahl von n |
| `subsampling.csv` | Echte Bilder, pro Design vielfach neu gezogen: wie oft das Ergebnis mit der Schlussfolgerung aus allen Bildern übereinstimmt, und wie stark sich das Verhältnis jedes n verschiebt |

**So liest man es:** Der Test über die n wird durch zwei Dinge begrenzt. Das erste ist das *Messrauschen*
(durch Felder und Wells), das mit mehr Bildern kleiner wird. Das zweite ist die *biologische Variation von n zu n*,
die nur mit mehr n kleiner wird. Mehr Bilder helfen nur, solange das Messrauschen noch in einer ähnlichen Größenordnung wie die
biologische Variation liegt. Das Rauschen von Well zu Well wird nur mit mehr Wells kleiner.

**Batch-Effekte:** Jedes n ist eine Platte an einem Tag; Platten- und Tageseffekte sind daher Teil der Variation von n zu n
und lassen sich nicht davon trennen. Effekte, die eine ganze Platte skalieren (Aussaat, Färbung, Lampe), heben sich im
Verhältnis zur Referenz auf. Effekte, die die Stärke der Reaktion ändern (z. B. Transfektionseffizienz), heben sich nicht auf;
der Bericht markiert ein n, dessen Behandlungseffekte durchgehend stärker oder schwächer sind. Eine Markierung ist ein Anlass,
die Labornotizen dieses Tages zu prüfen, kein Grund, das n auszuschließen.

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

## 9. Tumorexperimente (mCherry-Tumoren in Hirnschnitten)

Dasselbe Werkzeug wertet auch **Tumorexperimente** aus: mCherry-markierte Tumoren in Hirnschnitten, aufgenommen als
z-Stapel (`[experiment] type = "tumour_volume"`). Ein **Tumor** (Hemisphäre A oder B eines Schnitts) ist ein n; A und
B werden nie gemittelt. Begriffe: `CONTEXT.md` (Tumour experiments); technische Details: `docs/tumour.md`.

### 9.1 Eingabe
Ein Ordner pro Experiment in `data/`, darin **ein Unterordner pro Aufnahmesitzung** (ein Tag am Mikroskop), die
Dateien flach und wie geliefert:
```
data/tumour_TLR7/
├── 5d K07/                                     ← Sitzungsordner: sein Name gibt den Genotyp (K07 = TLR7 KO)
│   ├── 030826 HERV -clodro n3 A.ome.tiff        ← ein Original-z-Stapel = ein Tumor
│   ├── 030826 HERV -clodro n3 Af.tif            ← Fiji-Maske (Altverfahren, optional; nur zum Vergleich)
│   └── Statistics for 030826 HERV -clodro n3 Af.csv   ← Tabelle des Fiji 3D Objects Counter (optional)
└── 5d WT/
    └── 100926 cntrl +clodro n1 B.ome.tif
```
- **Der Genotyp kommt aus dem Namen des Sitzungsordners**, Wort für Wort gelesen (Groß-/Kleinschreibung egal): Er
  muss `WT` oder `K07` / `KO7` / `KO` / `TLR7 KO` / `TLR7-KO` / `TLR7KO` enthalten; andere Wörter (`5d`) werden als
  Sitzungsbezeichnung (`session_label`) festgehalten. Beliebig viele Sitzungsordner sind möglich, auch mehrere mit
  demselben Genotyp (Wiederholungssitzungen). Ein neuer Genotyp ist eine neue Stufe in den Einstellungen plus seine
  Ordner-Schreibweise als Alias. Ist der Genotyp eines Ordners nicht lesbar, stoppt der Lauf mit einer Meldung, wie
  es zu beheben ist. WT und KO wurden an verschiedenen Tagen aufgenommen, d. h. Genotyp und Aufnahmesitzung sind
  vermengt: Verhältnisse innerhalb einer Sitzung werden zwischen den Genotypen verglichen; ein direkter
  WT-gegen-KO-Unterschied wird gezeigt, aber nie getestet. Ob ein *Sitzungsordner* oder ein *Aufnahmetag* (aus den
  Mikroskop-Metadaten) als eine Charge gilt, legt `[tumour] batch` fest (Standard `session_folder`); der Bericht
  nennt, welche Definition verwendet wurde.
- **Umbenennen eines Sitzungsordners** verliert nichts: `kerncount run` findet die Dateien über ihren Namen, behält
  Ihre Änderungen in layout.csv und nutzt alle gespeicherten Messungen weiter; `kerncount init` aktualisiert danach
  layout.csv (die alte Fassung wird gesichert).
- Das 6-stellige Datum am Anfang eines Dateinamens ist optional und **nicht** das Aufnahmedatum; den Aufnahmetag
  liest kerncount aus den Metadaten jedes Stapels (`acquisition_day`).
- **Dateinamen** liefern Datum, Tumorlinie, Clodronat, Schnitt und Hemisphäre. Für **neue Aufnahmen**:
  `DDMMYY <Linie> <±clodro> n<Schnitt> <A|B>.ome.tif`, z. B. `100926 HERV +clodro n3 A.ome.tif`: Datum mit 6 Ziffern,
  Tumorlinie und Clodronat-Stufe genau wie in den Einstellungen (`cntrl`, `HERV`, `N4`; `-clodro`, `+clodro`), die
  Schnittbezeichnung `n<k>` (beginnt pro Gruppe neu; `s3` wird nicht als Schnitt gelesen), Hemisphäre `A` (oben)
  oder `B` (unten). Unterstriche statt Leerzeichen (`100926_HERV_+clodro_n3_A.ome.tif`), `.ome.tiff` und bekannte
  Tippfehler (`-clodron`) werden ebenfalls gelesen und in `layout.csv` vermerkt.
- Die Voxelgröße wird aus den OME-Metadaten jeder Datei gelesen (nie eingetippt); der z-Abstand darf je Datei
  verschieden sein.

### 9.2 Einrichten, Layout prüfen, Analyse starten
In der Weboberfläche das Experiment wählen, unter *Set up* **experiment type: Tumour volume (z-stacks)** einstellen
und **Set up** drücken (Befehlszeile: `kerncount init tumour_TLR7 --type tumour_volume`). Dann:
- **Layout**: eine Zeile pro Tumor mit den Faktorwerten (Genotyp, Tumorlinie, Clodronat), Schnitt, Hemisphäre,
  Tumor-ID (`TLR7 KO / HERV / -clodro / n3A`), include und note. Probleme, die eine Datei zu einem falschen
  Datenpunkt machen würden (unlesbarer Name, doppelte Tumor-ID, Sitzung ohne Genotyp, fehlende Datei), stehen oben
  und **blockieren die Analyse**, bis die Zelle korrigiert oder *include* abgewählt und in *note* begründet ist.
- **Settings** → *Analysis plan*: die Faktoren (Stufen, Referenzstufe, Platz in den Diagrammen), welche Messgrößen
  konfirmatorisch sind, und die Kontrastfamilien (unten). *Tumour measurement*: Glättung, Schwellenregel,
  Dicke des Mittelschnitts, Konturglättung, Schwellen-Sweep, QC- und Review-Einstellungen, jeweils mit Erklärung.
- **Run** (`kerncount run tumour_TLR7`): Optionen *skip review images*, *segment every stack again*, *threshold
  sweep*. Der erste Lauf braucht etwa 20–60 s pro Stapel; Segmentierungen werden zwischengespeichert, spätere Läufe
  sind schnell.

### 9.3 Messgrößen, und warum Fläche statt Volumen
Von jedem ~0,5–1 mm großen Tumor sind nur ~80–150 µm aufgenommen; ein 3D-Volumen misst daher vor allem, wie tief der
Stapel aufgenommen wurde. Die **konfirmatorischen Messgrößen** sind deshalb zweidimensional, auf einem
**Mittelschnitt** fester Dicke (20 µm um die Ebene, in der die Hauptmasse am größten ist – für jeden Tumor dieselbe
physikalische Schicht): die **Querschnittsfläche** (mm²) und der **Umfang** (mm) der Hauptmasse (des größten
zusammenhängenden Objekts). Alles andere ist **explorativ** (Effektgrößen mit KI, gekennzeichnet, ohne Sterne):
Grundfläche (Projektion), Schichtvolumen, Volumen und Oberfläche der Hauptmasse, Tumorlast (Hauptmasse + Satelliten),
Invasionsanteil, Elongation, Solidität, Randrauigkeit, Zirkularität und die Fiji-Altwerte. Eine automatische
Schwellenregel für jeden Tumor, nie pro Tumor von Hand angepasst.

**Sperre des Analyseplans.** Faktoren, Rollen der Messgrößen und Kontrastfamilien werden festgelegt, bevor
Ergebnisse angesehen werden (ADR 0003). Sobald Ergebnisse existieren, ist der Plan in der Weboberfläche **gesperrt**;
eine Änderung braucht **Unlock** mit einer Begründung. Das Entsperren und jede Änderung werden in `lock_log.csv`
geschrieben, und der nächste Lauf druckt den Verlauf (alter → neuer Plan-Hash) in den Bericht und das About-Blatt von
`results.xlsx`. Jede Statistikzeile trägt den Hash, mit dem sie berechnet wurde.

### 9.4 Ausgabe (`analysis/<Experiment>/`)
| Datei / Ordner | Inhalt |
|---|---|
| `results.xlsx` | alles in einer Arbeitsmappe: Tumoren, Statistik, geschätzte Mittelwerte, Familien, QC, Fluss, Sensitivitätsanalysen, Übereinstimmung mit Fiji, Planungshilfe (Download auf der Seite **Results**) |
| `stats_report.md` | der Bericht: Analyseplan und Sperre (mit Verlauf), Ergebnisse pro Messgröße und Familie, QC, Sensitivität, Fiji-Vergleich, Planungshilfe |
| `tumours.csv` | eine Zeile pro ausgewertetem Tumor: Faktorwerte, alle Messgrößen, QC-Markierungen |
| `excluded_tumours.csv`, `flow.csv` | ausgeschlossene Tumoren mit Grund (Layout oder Review); Dateien → Tumoren → ausgewertet |
| `stats_long.csv` | jede Schätzung in einer langen Tabelle (Seite **Statistics**; Sensitivitätszeilen nur auf Wunsch) |
| `plots/` | pro Messgröße: absolut, normiert (% von cntrl bei gleichem Genotyp und Clodronat), Effekte (Forest), Interaktion, A gegen B, Diagnostik, sowie mit Fiji-Daten / Sweep: gegen Fiji, Schwellen-Sweep; `shape_panel` |
| `prism/` | Prism-fertige Tabellen, eine Spalte pro Gruppe (absolut und normiert) |
| `review/` | die **verblindeten** Kontrollbilder und `corrections.csv` (unten) |

Statistik: ein lineares Modell auf log(Messgröße) über alle Faktoren; Effekte als Faktor-Änderungen (fold-change)
mit 95-%-KI, Holm innerhalb jeder Familie. Nur konfirmatorische Familien konfirmatorischer Messgrößen tragen Sterne.
Sensitivitätsanalysen (Welch-artig, Permutation, ohne QC-markierte Tumoren) stehen daneben, nie stattdessen.

### 9.5 Review und Korrekturen (verblindet)
`review/` enthält ein Bild pro Tumor, benannt nach einem zufälligen **Code** (nirgends Gruppe, Sitzung, Datum oder
Dateiname; der Schlüssel `review_key.csv` liegt außerhalb von `review/`): links die Maximumprojektion, rechts die
Ebene des Mittelschnitts; Hauptmasse orange, Satelliten cyan, Objekte nummeriert. Die Seiten **Review images** und
**Corrections** zeigen nur die Codes. In einer Runde (`manual_1` → g2, `manual_2` → g3, …) tragen Sie ein:
- eine **Objektnummer** aus dem Bild, um dieses Objekt zur Hauptmasse zu machen (`1` = zurück zur Automatik);
- **`x: Grund`**, um den Tumor auszuschließen (der Grund kommt in `excluded_tumours.csv`);
- nichts = keine Änderung.

Danach die Analyse erneut starten: Es wird nichts neu segmentiert; die Messgrößen des gewählten Objekts werden neu
berechnet. `review/corrections.csv` lässt sich auch in Excel bearbeiten, und für eine verblindete Kontrolle geben Sie
einer Kollegin nur den Ordner `review/`.

## 10. Entwicklung
```
pip install -e ".[test]"
pytest
```
Hinweise zu Algorithmen und Implementierung: `docs/technical.md`.

## Danksagung
Die Kernerkennung verwendet [StarDist](https://github.com/stardist/stardist) (Schmidt et al., MICCAI 2018;
Weigert et al., WACV 2020) mit dem vortrainierten Modell `2D_versatile_fluo`.
