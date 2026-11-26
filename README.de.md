# kerncount: automatisierte Zellkernzählung und Auswertung

*[English version](README.md)*

`kerncount` macht aus Ordnern mit Fluoreszenz-Mikroskopiebildern geprüfte Zählungen, Tabellen, Diagramme und
Statistiken. Es entstand als Ersatz für einen manuellen Arbeitsablauf in einem Labor (ImageJ/Fiji-Makros zum Zählen
der Kerne und Messen der Intensität, danach Übertragen der Zahlen in GraphPad Prism) und behält dessen Konventionen
bei: ein Ordner pro Versuchstag, Prism-fertige Tabellen und Vergleiche mit den alten Fiji-Ergebnissen, wo es sie
gibt. Fiji oder Prism brauchen Sie für kerncount nicht.

Es leistet Folgendes:

1. Zählt die Zellkerne in jedem Bild (DAPI oder eine andere Kernfärbung) mit einem Deep-Learning-Detektor (StarDist).
2. Schätzt ab, wie verlässlich jede Zählung ist, und erstellt für jedes Bildfeld ein Kontrollbild.
3. Ermöglicht es, Zählungen in einer einzigen CSV-Datei (oder im Browser) von Hand zu korrigieren. Jede
   Korrekturrunde wird dabei als eigene „Generation“ gespeichert.
4. Erstellt Tabellen, Prism-fertige CSV-Dateien, Diagramme und Statistiken für jede Behandlung, jedes n und alle n
   zusammen, mit einem hierarchischen (gemischten) Modell als Haupttest.
5. Misst einen Markerkanal pro Kern (z. B. Ki-67): % positive Kerne, Intensitäten, Flat-Field-Korrektur und einen
   „Ebenenvergleich“, der zeigt, welcher Auswerteschritt die Schlussfolgerung verändert.
6. Ergänzt zusätzliche Messgrößen pro Bild aus den erkannten Kernen (Anteil unscharfer Kerne, Mikronuklei,
   Kerngröße, ...).
7. Plant das nächste Experiment: wie viele n, Wells und Bilder Sie brauchen, und randomisierte Plattenbelegungen.
8. Misst fluoreszierende Tumoren in z-Stapeln von Hirnschnitten (Tumorexperimente, Abschnitt 9), mit verblindetem
   Review, gesperrtem Analyseplan und vorab festgelegter Statistik.

Alles lässt sich über eine Bedienoberfläche im Browser erledigen oder über die Kommandozeile.

## Inhalt

1. [Eingabe: Wie die Rohdaten organisiert sein müssen](#1-eingabe-wie-die-rohdaten-organisiert-sein-müssen)
2. [Installation](#2-installation-einmal-pro-computer)
3. [Bedienung](#3-bedienung)
4. [Einstellungen (`settings.toml`)](#4-einstellungen-settingstoml)
5. [Ausgabe: Was herauskommt und wie man es liest](#5-ausgabe-was-herauskommt-und-wie-man-es-liest)
6. [Typische Arbeitsabläufe](#6-typische-arbeitsabläufe)
7. [Fehlerbehebung](#7-fehlerbehebung)
8. [Hinweise und bekannte Einschränkungen](#8-hinweise-und-bekannte-einschränkungen)
9. [Tumorexperimente](#9-tumorexperimente-fluoreszierende-tumoren-in-hirnschnitten)
10. [Entwicklung](#10-entwicklung)

```
kerncount/                 (dieses Repository)
├── README.md              ← diese Anleitung (Englisch); README.de.md = deutsche Fassung
├── docs/                  ← ausführliche Dokumentation (siehe Abschnitt 10)
├── src/kerncount/         ← der Code
├── tests/                 ← Tests (pytest)
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

Begriffe in dieser Anleitung: Ein **n** ist ein unabhängiger Versuch (eine Platte, an einem Tag aufgenommen); eine
**Behandlung** (Bedingung) liegt in einem oder mehreren **Wells** dieser Platte; jedes Well wird in mehreren
**Bildfeldern** (Bildern) aufgenommen.

---

## 1. Eingabe: Wie die Rohdaten organisiert sein müssen

Dieser Abschnitt gilt für Zellzählexperimente (der Standard). Tumorexperimente sind anders aufgebaut (Abschnitt 9).

Jedes Experiment ist ein Ordner in `data/`:

```
data/<experiment>/                              z. B. data/my_experiment/
├── <ein Ordner pro n (Versuch)>/               z. B. "12.01.26 DAPI count"
│   ├── <ein Ordner pro Behandlung>/            z. B. "Control", "LyoVec veh. cntrl", "miRNA 12-5p"
│   │   ├── 1.tif
│   │   ├── 2.tif
│   │   └── ...
│   └── ...
└── ...
```

- **Ein Ordner pro n**, mit dem Datum im Namen (`12.01.26`, gelesen als Tag.Monat.Jahr; siehe `date_order`). Die n
  werden nach Datum nummeriert: `n1_2026-01-12`, `n2_2026-02-09` usw. Ein einzelnes n funktioniert ebenfalls. Ein
  Name mit nur Tag und Monat erhält das Jahr aus den Dateidaten (mit einem Hinweis zur Prüfung). Ein Zeitpunkt im
  Namen (`72h`, `48 h`, `day 3`) wird in die Spalte `timepoint` von `layout.csv` übernommen.
- **Ein Ordner pro Behandlung** innerhalb jedes n. Uneinheitliche Namen werden automatisch bereinigt; die Zuordnung
  können Sie in `layout.csv` überprüfen. Eingebaute Regeln: `LyoVec` (z. B. `lyo vec`, `LyoVec veh. cntrl`), `CL075`,
  microRNAs (`hsa miRNA 12-5p` / `12 5p` → `miR-12-5p`; eine bloße `12` übernimmt ihren Arm aus anderen n) und
  Kontrollen (`cntrl`, `ctrl`, `control` → `Control`). Alles Unbekannte wird unverändert übernommen und markiert.
- **Bilder sind nach Nummer benannt** (`1.tif`, `02.tif`; zusätzliche Zeichen wie in `` `10.tif `` sind kein Problem).
  Dateien ohne Nummer (z. B. `5 zoom.tif`) werden übersprungen und in `layout.csv` aufgeführt.
  Die Anzahl der Bilder pro Ordner ist beliebig.
- **Ignoriert werden:** Unterordner namens `output`, Dateien, deren Name `mask` enthält (alte Fiji-Ergebnisse), sowie
  zusätzliche Unterordnerebenen (Bilder werden in jeder Tiefe gefunden). Alte Fiji-Zähltabellen in einem
  Behandlungsordner werden zum Vergleich gelesen (`fiji_comparison.csv`).
- **Bildformat:** TIFF, 8 oder 16 Bit, Graustufen, RGB oder mehrkanalig, bei jeder Vergrößerung und von jeder Kamera.
  Die Kerngröße wird pro Ordner gemessen, und der Detektor skaliert entsprechend.
- **Kanäle als getrennte Dateien** (z. B. `1.tif` = Marker, `2.tif` = DAPI, abwechselnd) werden unterstützt:
  `[channels] files = "interleaved"` setzen (Abschnitt 3, *Einen Markerkanal messen*). `kerncount init` erkennt
  Dateien, die abwechselnd dunkel / hell sind, und schlägt dies vor.
- **Wells:** Ein n = eine Platte. Standardmäßig liegt jede Bedingung in 2 Wells, und Bild 1–6 = Well 1, Bild 7–12 =
  Well 2 (abgeleitet aus der höchsten Bildnummer; änderbar in `[design]` oder pro Bild in `wells.csv`).

---

## 2. Installation (einmal pro Computer)

### Am einfachsten: Ein-Klick-Programm (ohne Administratorrechte, ohne Python, ohne Befehle)
Laden Sie aus dem [neuesten Release](https://github.com/joe-bones/kerncount/releases/latest) herunter:
- **Windows:** `kerncount.exe`
- **Mac (Apple Silicon oder Intel):** `kerncount-mac.zip`. Ein Doppelklick entpackt es zu `kerncount.app`. Auf Intel-Macs wird TensorFlow 2.16 installiert, die letzte Version für Intel-Macs.

Legen Sie es in Ihren Arbeitsbereich-Ordner (neben `data`) und doppelklicken Sie darauf. Der erste Start braucht
Internet, einige Minuten und etwa 1 GB: Python 3.12 und alles, was kerncount benötigt, wird in Ihren Benutzerordner
installiert. Danach öffnet jeder Start kerncount sofort im Browser. Beim ersten Mal warnt das System vor einer
nicht signierten App: unter Windows **Weitere Informationen → Trotzdem ausführen**; am Mac **Fertig**, dann
**Systemeinstellungen → Datenschutz & Sicherheit → Trotzdem öffnen**. Aktualisieren = das Programm durch ein
neueres ersetzen. Details (wohin installiert wird, Reparatur mit `--reinstall`): [launcher/README.md](launcher/README.md).

Ein Ein-Klick-Programm enthält den Code des Releases, aus dem es stammt. Beachten Sie die Release-Hinweise; für
Änderungen nach dem neuesten Release installieren Sie aus dem Repository (unten).

### Alternativ: Installation aus diesem Repository (für die Kommandozeile)
1. **Mac/Linux:** Installieren Sie **Python 3.12** von python.org (`setup.command` sucht nach `python3.12`).
   kerncount selbst braucht Python 3.11 oder neuer; 3.12 wird wegen TensorFlow empfohlen.
   **Windows:** Keine Vorinstallation und keine Administratorrechte nötig: `setup.bat` holt Python 3.12 selbst
   (über [uv](https://docs.astral.sh/uv/)) in Ihren Benutzerordner (`%LOCALAPPDATA%\kc`).
2. Laden Sie dieses Repository herunter (grüne Schaltfläche **Code** → *Download ZIP*, oder `git clone`).
3. Doppelklicken Sie auf `setup.bat` (Windows) oder führen Sie `setup.command` aus (Mac/Linux).
   Das dauert einige Minuten und benötigt etwa 1 GB (vor allem TensorFlow). Das StarDist-Modell (ca. 5 MB) wird
   beim ersten Lauf automatisch heruntergeladen; dafür ist eine Internetverbindung nötig.

Oder von Hand, in einer beliebigen Python-Umgebung ab 3.11: `pip install .` (oder `pip install -e .`, um den Code
zu bearbeiten).

Befehle werden aus Ihrem Arbeitsbereich-Ordner ausgeführt. Unter Windows beginnt jeder der folgenden Befehle statt mit
`kerncount` mit dem vollständigen Pfad zu `kerncount.bat` in diesem Ordner (z. B. `C:\tools\kerncount\kerncount.bat run my_experiment`).
Unter Mac/Linux lautet er `<dieser Ordner>/.venv/bin/kerncount`. Wenn Sie diesen Ordner zu Ihrem PATH hinzufügen, genügt `kerncount`.

Wenn TensorFlow nicht lädt (unter Windows meist eine fehlende Microsoft-Visual-C++-Laufzeit; `setup.bat` meldet
das), öffnet sich kerncount trotzdem, und der Detektor `classical` (`[segmentation] backend = "classical"`)
funktioniert.

---

## 3. Bedienung

Führen Sie jeden Befehl aus dem Arbeitsbereich-Ordner aus (dem Ordner, der `data/` enthält). Experimente werden
über ihren Ordnernamen in `data/` angesprochen.

### Ohne Befehle: die Weboberfläche
Doppelklicken Sie auf **`kerncount.exe`** / **`kerncount.app`** (Abschnitt 2) oder, bei einer Repository-Installation,
auf **`kerncount-gui.bat`** (Windows) bzw. **`kerncount-gui.command`** (Mac). Dasselbe aus einem Terminal:
`kerncount gui` (Optionen: ein Arbeitsbereich-Ordner, `--port` (Standard 8765), `--no-browser`). Wird die Seite
nicht innerhalb eines Arbeitsbereichs gestartet, können Sie dort einen auswählen oder anlegen.

Die linke Seitenleiste listet die Experimente; das geöffnete Experiment klappt in diese Seiten auf:

| Seite | Was Sie dort tun |
|---|---|
| **Experiment home** | Das Experiment einrichten (mit seinem **Experimenttyp**), die Analyse starten (die Befehlsoptionen unten als Häkchen), Generationen neu bauen, die Sampling-Analyse starten, eine Wells-Vorlage schreiben; das Feld **Attention** (Duplikate, fehlende Bedingungen, nicht bestandene Prüfungen, seit dem letzten Lauf geänderte Einstellungen, ...); die n, die Design-Übersicht und die Generationen. Live-Protokoll und **Stop** unten. |
| **Results explorer** | Alle Diagramme pro Generation, nach Ebene (alle n / pro n / pro Behandlung), Messgröße und Ansicht (absolut, relativ, Modell, Vergleiche, Robustheit, Diagnostik, Design, Tabellen); `results.xlsx` herunterladen |
| **Statistics table** | Jedes Testergebnis (`stats_long.csv`) als filter- und sortierbare Tabelle; die gefilterten Zeilen exportieren |
| **Layer comparison** | (nur Markerkanäle) Welcher Auswerteschritt die Schlussfolgerung verändert (Abschnitt 5.4) |
| **Intensity & cut-off** | (nur Markerkanäle) Positivitäts-Schwellen, Flat-Field-Karten, Verteilungen, Kanalprüfungen und eine **verblindete Tuning-Ansicht**, um eine Schwelle zu wählen und festzuschreiben |
| **Review images** | Übersicht pro Behandlung, dann jedes Bild in voller Größe (← / → zum Blättern); Korrektur direkt eintragen |
| **Corrections** | Eine Generation bearbeiten: korrigierte Zählungen oder **x** für eine Runde eintragen (`manual_1` → g2, …), Wert vorher und nachher sehen, dann die Generation **bauen** (etwa eine Minute; die Erkennung wird wiederverwendet) |
| **Planner** | Das nächste Experiment interaktiv planen: Power, KI-Breite und Aufwand für beliebige Zahlen von n, Wells und Bildfeldern |
| **Plate maps** | Randomisierte Plattenbelegungen für das nächste Experiment (PNG/PDF/CSV) |
| **Layout** | n- / Behandlungsnamen, Datum und Zeitpunkt korrigieren; Ordner ein- oder ausschließen |
| **Settings** | Alle Einstellungen aus Abschnitt 4 als Formular mit Erklärungen; oder `settings.toml` als Text bearbeiten. Nach dem Speichern zeigt die Seite, was neu laufen muss, mit einer Start-Schaltfläche. |
| **All files** | Alles im Projektordner durchsuchen und ansehen |

Jede Ansicht hat einen eigenen Link (in der Adresszeile), sodass ein Lesezeichen genau diese Ansicht wieder öffnet.
Tumorexperimente (Abschnitt 9) nutzen dieselben Seiten mit eigenem Inhalt; Planner und Plate maps werden für sie
nicht angezeigt.

Lassen Sie das schwarze (Windows) bzw. Terminal-Fenster (Mac) während der Arbeit offen; schließen Sie es oder
drücken Sie **Quit** zum Beenden. Die Seite funktioniert nur auf diesem Computer, und es wird nichts nach außen
gesendet. Alles, was eine möglicherweise von Ihnen bearbeitete Datei überschreibt, fragt vorher nach, und die alte
Datei wird nach `_work/backups/` kopiert.

### Schritt 1: Projekt anlegen (einmal pro Experiment)
```
kerncount init my_experiment
```
Dies durchsucht `data/my_experiment/` und legt `analysis/my_experiment/` mit zwei Dateien an
(`kerncount init --all` tut dies für jedes Experiment in `data/`):

- **`layout.csv`** enthält eine Zeile pro Rohdatenordner. Prüfen Sie die Spalten `n`, `date`, `timepoint` und
  `treatment`. Um einen Namen zu korrigieren, bearbeiten Sie ihn; um einen Ordner auszuschließen, setzen Sie `include`
  auf `no`. Die Spalte `note` zeigt Auffälligkeiten (übersprungene Dateien, abgeleitete Namen, ungepaarte
  Kanaldateien). Sie können Spalten mit den Namen der `[experiment]`-Einstellungen ergänzen (z. B. `cell_line`), um
  diese pro n oder pro Behandlung zu überschreiben.
- **`settings.toml`** enthält alle Einstellungen (Abschnitt 4). Die Standardwerte sind sinnvoll gewählt; füllen Sie
  bei Bedarf die Beschreibung unter `[experiment]` aus (Zelllinie, Vergrößerung, Färbung, ...) und setzen Sie
  `reference` / `normalise_to` auf die Namen Ihrer Kontrollgruppen.

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
(links: Original; rechts: was gezählt wurde). Beginnen Sie mit den roten. Sehen Sie sich auch `design_check.csv`
und `duplicates.csv` (falls vorhanden) oder das Feld **Attention** an: Sie markieren fehlende Bedingungen,
ungewöhnliche Well-Größen und bytegleiche Bilder in verschiedenen Bedingungen.

### Schritt 4: Zählungen von Hand korrigieren (optional)
Öffnen Sie **`analysis/my_experiment/counts.csv`** in Excel (oder nutzen Sie die Seite **Corrections**):

| n | treatment | image | auto_count | qc_ok | qc_level | manual_1 | manual_2 | manual_3 |
|---|---|---|---|---|---|---|---|---|
| n1_2026-01-12 | Control | 3 | 85 | yes | high | | | |
| n1_2026-01-12 | miR-12-5p | 6 | 140 | no | low | **121** | | |
| n1_2026-01-12 | miR-12-5p | 9 | 4 | no | low | **x** | | |

- Tragen Sie in `manual_1` eine korrigierte Zahl ein, oder **`x`**, um das Bild auszuschließen. Lassen Sie alle anderen Zellen leer.
- Speichern Sie und führen Sie `kerncount run my_experiment` erneut aus (oder das schnellere
  `kerncount regen my_experiment`). Dadurch entsteht Generation **g2** (siehe Abschnitt 5.3).
- Eine spätere Korrekturrunde kommt in `manual_2`, woraus g3 entsteht, und so weiter.

### Schritt 5: Wie viele Bilder brauchen wir, und wie geht es weiter? (optional)
```
kerncount sampling my_experiment
```
Die Ergebnisse werden nach `analysis/my_experiment/sampling/g<N>/` geschrieben (Abschnitt 5.5): ein Rückblick
(„Hätten weniger Bilder gereicht?“) und ein Plan für das nächste Experiment (wie viele n, Wells und Bildfelder), den
die Seite **Planner** interaktiv zeigt. Dafür sind mindestens 2 n mit der Referenzgruppe nötig. Danach, für
randomisierte Plattenbelegungen:
```
kerncount platemap my_experiment --plates 3 --format 24
```
Dies schreibt `analysis/my_experiment/plans/platemap_<Datum-Uhrzeit>/` mit einer Belegung pro Platte (PNG, PDF, CSV)
und einer `README.txt` mit dem Seed und dem Befehl, der dieselben Belegungen wieder erzeugt. Details:
[docs/planner.md](docs/planner.md).

### Einen Markerkanal messen (z. B. Ki-67 in C2)
Für mehrkanalige TIFFs (C1 = DAPI, C2 = Marker) in `settings.toml` einstellen:
```
[channels]
names = ["DAPI", "Ki67"]
segment = "DAPI"       # Kerne im DAPI-Kanal finden (Name oder Nummer)
```
Ist jeder Kanal als eigene Datei gespeichert, abwechselnd (1.tif = Ki-67, 2.tif = DAPI, 3.tif = Ki-67, ...):
```
[channels]
files = "interleaved"
file_order = ["Ki67", "DAPI"]
segment = "DAPI"
```
und danach das Layout mit `kerncount init my_experiment --force-layout` neu schreiben. Jedes Dateipaar wird auf
Verschiebungen und Fehlpaarungen geprüft (`intensity/alignment.csv`).

Dann wie gewohnt ausführen. Jeder Kanal außer dem Kernkanal wird **pro Kern** gemessen (`[intensity]
per_nucleus`), auf ungleichmäßige Ausleuchtung (Flat Field) und lokalen Hintergrund korrigiert und als positiv oder
negativ eingestuft (`[positivity]`). Die Hauptmessgröße ist **% positive Kerne** pro Bild; mittlere Intensitäten
sind zweitrangig. Die Schwelle (pro n und Kanal) stammt in dieser Reihenfolge aus Negativkontroll-Wells in diesem n
(`[experiment] negative_control`, z. B. nur Sekundärantikörper), aus einem von Ihnen gesetzten Wert
(`manual_cutoff`) oder einem in der verblindeten Tuning-Ansicht des Browsers festgeschriebenen Wert, oder aus einer
Nachweisgrenzen-Regel. Ohne Negativkontrolle bedeutet „positiv“ „nachweisbar über dem Hintergrund gefärbt“. Die
Ergebnisse landen in `intensity/` und durchlaufen dieselbe Statistik wie die Zählungen. Details:
[docs/intensity.md](docs/intensity.md).

Die älteren regionsbasierten Mittelwerte gibt es weiterhin: `[intensity] regions = ["nucleus", "ring", "cell",
"fiji_cell"]` ergibt Spalten `intensity_<channel>_<region>` (`fiji_cell` = Kern um 6 px vergrößert, wie im alten
Fiji-Makro). Um statt Kernen ganze Zellen aus einem zytoplasmatischen Kanal zu segmentieren, setzen Sie `segment`
auf diesen Kanal und `[segmentation] backend = "classical"`.

### Wells und Feldpositionen (optional)
kerncount muss wissen, aus welchem Well jedes Bild stammt, weil Well-zu-Well- und Feld-zu-Feld-Streuung getrennt
behandelt werden. Die Standardregel (Bild 1–6 = Well 1, 7–12 = Well 2) wird in `[design]` eingestellt. Für alles,
was die Regeln nicht abbilden:
```
kerncount wells-template my_experiment
```
schreibt `wells.csv` (eine Zeile pro Bild: `well`, `position`; oder `wells.new.csv`, falls es schon eine gibt).
Bearbeiten Sie sie und setzen Sie `[design] well_layout = "file"`. Details: [docs/design.md](docs/design.md).

### Befehlsoptionen
| Befehl | Funktion |
|---|---|
| `kerncount init <experiment>` | Legt `analysis/<experiment>/` aus `data/<experiment>/` an. Eine bestehende `layout.csv` bleibt erhalten, ein neuer Scan wird nach `layout.new.csv` geschrieben. |
| `kerncount init --all` | Dasselbe für jedes Experiment in `data/` (bestehende Projekte bleiben erhalten) |
| `kerncount init <experiment> --force` | Neu durchsuchen und `layout.csv` sowie `settings.toml` überschreiben (alte Fassungen kommen nach `_work/backups/`) |
| `kerncount init <experiment> --force-layout` | `settings.toml` behalten, `layout.csv` neu schreiben (z. B. nach Änderung von `[channels] files`) |
| `kerncount init <experiment> --type tumour_volume` | Ein Projekt eines anderen Experimenttyps anlegen (Standard `cell_count`). Mit `--all` werden Experimente, die schon einen anderen Typ haben, gemeldet und übersprungen. |
| `kerncount run <experiment>` / `--all` | Zählen, prüfen und auswerten für ein / alle Experimente |
| `kerncount run <experiment> --rebuild g2` | Einen bestehenden Generationsordner überschreiben (geschieht standardmäßig nie). `--rebuild all` für alle. |
| `kerncount run <experiment> --resegment` | Cache ignorieren und jeden Kern neu erkennen |
| `kerncount run <experiment> --no-review` | Kontrollbilder nicht erstellen (schneller) |
| `kerncount run <experiment> --unlock "Grund"` | Tumorexperimente: einen geänderten Analyseplan mit protokolliertem Grund annehmen (Abschnitt 9.3) |
| `kerncount regen <experiment> [--rebuild g2 g3]` | Generationen nur aus `counts.csv` neu bauen (schnell; nach Korrekturen oder einer geänderten Statistik-Einstellung) |
| `kerncount sampling <experiment> [--generation g2] [--iterations 1000]` | Sampling-Analyse und Planer für eine Generation (Standard: die neueste) |
| `kerncount wells-template <experiment> [--force]` | `wells.csv` aus der aktuellen Well-Zuordnung schreiben (`--force` überschreibt eine bestehende) |
| `kerncount platemap <experiment>` | Randomisierte Plattenbelegungen. Optionen: `--plates` (Standard 3), `--format` 6/12/24/48/96, `--wells` pro Behandlung (Standard 2), `--reference-wells` (eine Zahl oder `auto`), `--avoid-edge`, `--seed`, `--no-balance`, `--negative-control NAME`, `--control-wells` 1/2 |
| `kerncount gui [workspace] [--port N] [--no-browser]` | Die Weboberfläche öffnen |

`regen`, `sampling`, `wells-template` und `platemap` (jeweils auch mit `--all`) sind für Zellzählexperimente.
Für Ordner außerhalb eines Arbeitsbereichs funktionieren auch Pfade: `kerncount init <Rohdatenordner> <Projektordner>`,
`kerncount run <Projektordner>`. Mit der Umgebungsvariable `KERNCOUNT_WORKSPACE` können Sie Experimente
von überall aus über ihren Namen ausführen. `kerncount <Befehl> --help` listet alle Optionen.

---

## 4. Einstellungen (`settings.toml`)

Mit jedem Texteditor oder auf der Seite **Settings** bearbeitbar. Eine gelöschte Einstellung fällt auf ihren
Standardwert zurück. Die Tabellen unten gelten für Zellzählexperimente; Tumorexperimente haben eigene Einstellungen
(Abschnitt 9).

### `[project]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `raw_dir` | (von init gesetzt) | Wo die Rohbilder liegen. Relativ zum Projektordner gespeichert (`../../data/<experiment>`), sodass der gesamte Arbeitsbereich verschoben oder kopiert werden kann. |
| `date_order` | `"DMY"` | Wie Datumsangaben in den Namen der n-Ordner geschrieben sind (`DMY`, `MDY`, `YMD`) |
| `csv_delimiter` | `","` | Verwenden Sie `";"`, wenn Ihr Excel (deutsche/europäische Ländereinstellung) CSV-Dateien in einer einzigen Spalte anzeigt. Beim Einlesen wird beides automatisch erkannt. |

### `[experiment]`: was aufgenommen wurde
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `type` | `"cell_count"` | Experimenttyp: `cell_count` (diese Anleitung) oder `tumour_volume` (Abschnitt 9). Wird von `init --type` gesetzt. |
| `cell_line` | `""` | z. B. `"LN229"`. n verschiedener Zelllinien werden getrennt ausgewertet, nie zusammengelegt. |
| `magnification`, `camera_pixel_um` | `""` | Objektiv (z. B. `20`) und Pixelgröße der Kamera in µm; zusammen ergeben sie die Pixelgröße |
| `seeding_density`, `seeding_unit` | `""`, `"cells/well"` | Wird in den Ausgaben festgehalten |
| `stain`, `notes` | `""` | Freitext, wird in den Ausgaben festgehalten |
| `negative_control` | `""` | Name einer Negativkontroll-Behandlung (z. B. Wells nur mit Sekundärantikörper). Dient nur zum Festlegen der Positivitäts-Schwellen; geht in keine Statistik ein. |
| `negative_control_type` | `"secondary_only"` | `secondary_only`, `isotype` oder `unstained` |

Jeder dieser Werte lässt sich pro n oder pro Behandlung mit einer gleichnamigen Spalte in `layout.csv` überschreiben.

### `[channels]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `names` | `[]` | Kanalnamen, z. B. `["DAPI", "Ki67"]`. Leer = C1, C2, ... |
| `segment` | `1` | Kanal, in dem die Kerne gesucht werden: eine Nummer (1 = erster) oder ein Name |
| `files` | `"single"` | `single` = alle Kanäle in einer Datei; `interleaved` = eine Datei pro Kanal, abwechselnd |
| `file_order` | `[]` | `interleaved`: Kanalnamen in Dateireihenfolge, z. B. `["Ki67", "DAPI"]` |
| `align` | `"check"` | `interleaved`: `check` (Verschiebungen zwischen den Dateien eines Felds melden), `correct` (kleine Verschiebungen auch zurückschieben), `off` |
| `align_warn_px` | `2.0` | Warnen, wenn die Kanäle eines Felds um mehr als diesen Wert verschoben sind |
| `align_max_correct_px` | `10.0` | `correct` korrigiert nur Verschiebungen bis zu diesem Wert; größere werden markiert |

### `[segmentation]`: Kerne finden
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `backend` | `"stardist"` | `stardist` = Deep-Learning-Detektor (empfohlen); `classical` = Schwellenwert + Watershed wie in Fiji |
| `nucleus_diameter_px` | `"auto"` | `auto` misst die Kerngröße pro Ordner (passt sich der Vergrößerung an). Alternativ eine feste Zahl. |
| `target_diameter_px` | `30` | Größe, bei der der Detektor am besten arbeitet; Bilder werden intern darauf skaliert |
| `prob_thresh` | `"default"` | Konfidenzschwelle des Detektors (niedriger = mehr Objekte) |
| `pixel_size_um` | `""` | µm pro Pixel, falls bekannt (aus einem Objektmikrometer). Ergänzt µm² und Zellen/mm². Leer = `camera_pixel_um / magnification`, wenn beide gesetzt sind. Eine Änderung löst keine Neuerkennung aus. |

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

### `[intensity]`: Intensitätsmessung
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `enabled` | `true` | Intensitäten messen |
| `measure` | `"all"` | Kanäle für die Regions-Mittelwerte, z. B. `[2]` nur für Kanal 2 |
| `regions` | `["nucleus", "cell"]` | `nucleus`; `ring` (Band um den Kern, also Zytoplasma); `cell` (Kern + Ring); `fiji_cell` (Kern + `fiji_ring_px`, das alte Fiji-Makro) |
| `ring_px` | `"auto"` | Breite des Bands; `auto` = 0,25 × Kerndurchmesser |
| `subtract_background` | `true` | Den Hintergrund jedes Bildes abziehen (Median der leeren Bereiche) |
| `per_nucleus` | `"auto"` | Pro Kern gemessene Markerkanäle: `auto` = jeder Kanal außer dem Kernkanal; eine Liste; oder `"off"` |
| `erode_px` | `"auto"` | Den Kernumriss vor dem Messen verkleinern (`auto` = 1 px, 2 px bei Kernen ≥ 16 px Durchmesser) |
| `bg_gap_px`, `bg_ring_px`, `bg_min_px` | `3`, `"auto"`, `40` | Lokaler Hintergrundring pro Kern: Abstand zum Kern, Breite (`auto` = 0,5 × Durchmesser), Mindestzahl an Pixeln, bevor auf einen breiteren Ring / den Bildhintergrund ausgewichen wird |
| `fiji_ring_px` | `6` | Vergrößerung für die Region `fiji_cell` |
| `flatfield` | `"auto"` | Korrektur ungleichmäßiger Ausleuchtung: `auto` (pro n und Kanal aus den Bildern dieses n geschätzt), `off` oder ein Ordner mit Referenzbildern |
| `flatfield_dark` | `true` | Das additive Hintergrundmuster (Kamera-Offset, Autofluoreszenz) vor der Schätzung des Flat Fields entfernen |

### `[readouts]`: zusätzliche Messgrößen pro Bild
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `enabled` | `true` | Spalten `readout_*` aus den erkannten Kernen berechnen (keine Neuerkennung) |
| `include` | `"all"` | `all` oder eine Liste wie `["frac_out_of_focus", "coverage"]` |

Die Messgrößen sind `frac_out_of_focus`, `micronuclei_per_nucleus`, `nuclear_area_median` (und `_um2`),
`nuclear_area_cv`, `frac_irregular`, `coverage` und `nn_distance_median`; sie helfen, Zelltod, Schädigung und
vermindertes Wachstum auseinanderzuhalten. Details: [docs/readouts.md](docs/readouts.md).

### `[positivity]`: % positive Kerne pro Markerkanal
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `enabled` | `true` | Kerne als positiv / negativ einstufen |
| `channels` | `"auto"` | `auto` = jeder pro Kern gemessene Kanal, oder eine Liste von Namen |
| `control_percentile` | `99.0` | Mit Negativkontroll-Wells: Schwelle = dieses Perzentil ihrer Kerne |
| `manual_cutoff` | `""` | Eine Zahl (alle Kanäle) oder z. B. `{ Ki67 = 1200 }`. Leer = eine im Browser festgeschriebene Schwelle, sonst automatisch. Vorhandene Negativkontroll-Wells haben Vorrang. |
| `lod_k` | `3.0` | Automatische Regel (Nachweisgrenze): positiv, wenn das Signal über k × dem Hintergrundrauschen liegt |
| `min_nuclei` | `5` | Bilder mit weniger messbaren Kernen erhalten keinen Wert für % positiv |

### `[analysis]`: Statistik
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `reference` | `"Control"` | Hauptvergleichsgruppe (Hauptdiagramme) |
| `normalise_to` | `["Control", "LyoVec"]` | Kontrollgruppen, auf die die Ergebnisse in Prozent bezogen werden (innerhalb jedes n) und gegen die getestet wird. Jede aufgeführte Gruppe erhält eigene Tests und Diagramme; eine aufgeführte Gruppe, die in Ihren Daten fehlt, erhält keine. |
| `treatment_order` | `["Control", "LyoVec", "CL075"]` | Reihenfolge in Tabellen und Diagrammen; andere Behandlungen folgen alphabetisch |
| `exclude_low_qc` | `false` | Bei `true` werden Bilder mit QC niedrig (und ohne manuellen Wert) aus der Statistik ausgeschlossen |
| `two_way_anova` | `true` | Zusätzliche zweifaktorielle ANOVA mit Tukey-Vergleichen aller Paare: pro n Behandlung × Well, alle n Behandlung × n (5.4.3) |
| `primary_test` | `"mixed"` | Was die Sterne in den Diagrammen über alle n setzt: `mixed` (hierarchisches Modell, log-Skala) oder `rm_anova` (Prisms ANOVA mit Messwiederholung + Dunnett) |
| `primary_family` | `""` | Die Kontrolle, deren Vergleiche die Hauptaussage sind; leer = `reference` |
| `vehicle_check` | `"auto"` | Äquivalenztest für eine Vehikelkontrolle, z. B. `"LyoVec vs Control"`; `auto` = jede andere Kontrolle gegen `reference`; `""` = aus |
| `equivalence_margin` | `0.2` | Vehikelprüfung: das Verhältnis muss zwischen 1/1,2 und 1,2 liegen, um als äquivalent zu gelten |
| `equivalence_margin_pp` | `5.0` | Vehikelprüfung für % positiv: Differenz innerhalb von ± 5 Prozentpunkten |
| `positive_control` | `""` | z. B. `"CL075"`: jede Behandlung, die keine Kontrolle ist, wird zusätzlich mit ihr verglichen |
| `model_unit` | `"well"` | Einheit des gemischten Modells: `well` (Well-Mittelwerte; Verhältnisse der mittleren Kerne pro Feld) oder `field` |
| `classes` | `"auto"` | Gruppen von Behandlungen mit eigener n-zu-n-Varianz: `auto` = eine pro Behandlung; `{}` = Kontrollen gegen alle anderen; oder z. B. `{ mirna = ["miR-*"] }` |
| `density_covariate` | `"readout_coverage"` | Für %-Messgrößen zusätzlich eine für die Zelldichte adjustierte Schätzung; `""` = aus |
| `among_classes` | `[]` | Klassen, deren Mitglieder zusätzlich untereinander verglichen werden (alle Paare), z. B. `["mirna"]` |

### `[acceptance]`: Prüfungen pro n
Werden an den Zählungen jedes n geprüft und in jeder Generation in `acceptance.csv` geschrieben. Standardmäßig
werden sie nur berichtet.

| Einstellung | Standard | Bedeutung |
|---|---|---|
| `vehicle_within` | `0.3` | Vehikel / Kontrolle muss in diesem n innerhalb von 1 ± 0,3 liegen (`""` = aus) |
| `positive_control_max_ratio` | `""` | z. B. `0.8`: Positivkontrolle / Referenz muss dieses Verhältnis erreichen (`""` = aus) |
| `min_images_per_condition` | `3` | Jede Bedingung braucht mindestens so viele verwertbare Bilder |
| `max_low_qc_fraction` | `0.5` | Höchstens dieser Anteil der Bilder des n darf QC niedrig sein |
| `no_duplicate_images` | `true` | Ein n fällt durch, wenn eines seiner Rohbilder bytegleich mit einem Bild einer anderen Bedingung ist |
| `enforce` | `false` | `true`: durchgefallene n werden aus der Statistik über alle n ausgeschlossen (und in den Ausgaben genannt) |

### `[design]`: Wells und Bildfelder
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `wells` | `2` | Aufgenommene Wells pro Bedingung und n (eine Platte pro n); `"auto"` = aus `wells.csv`, falls vorhanden, sonst 2 |
| `fields_per_well` | `"auto"` | `auto` = pro Bedingung höchste Bildnummer ÷ Wells (aufgerundet), also 18 Bilder = 9 + 9; oder eine feste Zahl |
| `well_layout` | `"consecutive"` | `consecutive` = Bild 1–6 Well 1, 7–12 Well 2; `alternating` = ungerade/gerade; `file` = `wells.csv` lesen |
| `field_positions` | `"fixed"` | `fixed` = Feld k ist in jedem Well dieselbe Stelle (Positionseffekte werden modelliert); `random` = die Position hat keine Bedeutung |

Ältere Einstellungsdateien mit `plates` / `fields_per_plate` / `plate_layout` funktionieren weiterhin; der Bericht bittet darum, sie umzubenennen.

### `[review]`
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `overlay_alpha` | `0.35` | Transparenz der farbigen Überlagerung |
| `jpeg_quality` | `88` | Qualität der Kontrollbilder |

### `[planner]`: das nächste Experiment planen
| Einstellung | Standard | Bedeutung |
|---|---|---|
| `effects` | `[0.8, 0.7, 0.5]` | Verhältnisse Behandlung / Referenz, für die geplant wird (0,7 = 30 % weniger Kerne) |
| `target_power` | `0.8` | Wahrscheinlichkeit eines signifikanten Ergebnisses, wenn der Effekt real ist |
| `alpha` | `0.05` | Signifikanzniveau |
| `cost_per_n`, `cost_per_well`, `cost_per_field` | `240`, `5`, `1` | Aufwand pro n, pro zusätzlichem Well und pro Bild (beliebige Einheit, z. B. Minuten); damit wird das günstigste Design gesucht |
| `max_n` | `12` | Größte Anzahl n in den Design-Tabellen |
| `plate_format` | `24` | 6, 12, 24, 48 oder 96 Wells |
| `avoid_edge` | `true` | Plattenbelegungen: zuerst innere Wells, Rand-Wells nur bei Bedarf |
| `reference_wells` | `"auto"` | Wells pro Referenzbedingung; `auto` = Dunnett-Regel (Wells × √Vergleiche) |
| `adjust_for_family` | `true` | Power = Wahrscheinlichkeit eines Sterns im Hauptdiagramm (für alle Vergleiche adjustiert); `false` = ein einzelner unadjustierter Test |

Nach einer Änderung der Einstellungen sagt Ihnen die Seite **Settings**, was neu laufen muss: Erkennungs-Einstellungen
erkennen jeden Kern neu (langsam); die meisten anderen brauchen einen Lauf mit `--rebuild all`; `[analysis]` /
`[acceptance]` nur `regen --rebuild all`; `[planner]` nur die Sampling-Analyse.

---

## 5. Ausgabe: Was herauskommt und wie man es liest

```
analysis/<experiment>/
├── settings.toml, layout.csv     Ihre Einstellungen und die Ordnerzuordnung (layout.new.csv = neuer Scan zum Vergleich)
├── counts.csv                    ★ DIE Arbeitsdatei: eine Zeile pro Bild, alle n, Spalten für manuelle Korrekturen
├── qc_details.csv                alle Messwerte pro Bild (siehe 5.2)
├── design_check.csv              eine Zeile pro n × Behandlung: Wells, Bildfelder, fehlende Bilder, Hinweise
├── duplicates.csv                Rohbilder mit identischem Inhalt (falls vorhanden)
├── fiji_comparison.csv           alte Fiji-Zählung vs. neue Zählung, pro Bild (nur wenn alte Fiji-Ergebnisse gefunden wurden)
├── wells.csv                     (optional) welches Well / welche Position jedes Bild ist
├── review/                       Bilder zur Kontrolle der Zählung (5.1)
├── intensity/                    Markerkanäle: Schwellen, Flat Field, Verteilungen, Kanalprüfungen
├── generations/g1, g2, ...       Ergebnisse (5.3, 5.4)
├── sampling/g1/                  wie viele Bilder nötig sind, und der Planer (5.5)
├── plans/                        Plattenbelegungen
└── _work/                        Cache, cells.csv und nuclei_intensity.csv (Werte pro Zelle / pro Kern),
                                  backups/ jeder Datei, die kerncount überschrieben hat (Löschen von _work erzwingt Neuerkennung)
```

### 5.1 Kontrollbilder: `review/`
```
review/<n>/<treatment>_overview.jpg     alle Bildfelder einer Behandlung; Rahmen = QC-Stufe
review/<n>/<treatment>/<NN>.jpg         ein Bildfeld: Original | Überlagerung
```
Die Kopfzeile nennt die Zählung, den plausiblen Bereich, die QC-Stufe und die Gründe. Farben der Überlagerung:

| Farbe | Bedeutung | Gezählt? |
|---|---|---|
| grün | Kern | ja |
| orange | mögliche Verschmelzung von 2+ sich berührenden Kernen | ja, als 1 (erhöht die Obergrenze des Bereichs) |
| magenta | Debris, Mikronukleus oder Fragment | nein |
| blau | unscharf (außerhalb der Fokusebene) | nein |
| gelb | vom unteren/rechten Bildrand angeschnitten | nein (Randregel) |

### 5.2 Werte pro Bild: `counts.csv` und `qc_details.csv`
- **`auto_count`**: die automatische Zählung.
- **`count_low` – `count_high`** (*plausibler Bereich*): was die Zählung wäre, wenn jede Grenzfallentscheidung
  (nahe der Größen- oder Fokusschwelle, geringe Konfidenz des Detektors, mögliche Verschmelzungen) andersherum
  ausgefallen wäre. Ein enger Bereich bedeutet ein eindeutiges Bild.
- **`disagreement`**: Abweichung von einer unabhängigen klassischen Zählung (Schwellenwert + Watershed).
  Große Werte bedeuten ein schwieriges Bild (dicht, verklumpt, helles Zytoplasma).
- **`qc_level`** (high / medium / low) fasst das Obige zusammen. **Niedrig bedeutet „dieses Bild ansehen“, nicht
  „diese Zählung ist falsch“.** Bilder mit niedriger QC häufen sich oft in dünn besiedelten oder geschädigten Wells,
  wo das Zählen tatsächlich schwieriger ist. `legacy` = kein Rohbild, nur eine alte Fiji-Zählung.
- `n_debris`, `n_out_of_focus`, `n_edge`, `n_merged`: wie viele Objekte jedes Typs ausgeschlossen oder markiert wurden.
- `diameter_px`: typische Kerngröße in diesem Ordner (nützlich als Plausibilitätsprüfung für Vergrößerung/Biologie).
- `fiji_count`: was das alte Makro ergab, zum Vergleich (falls vorhanden).
- `intensity_<channel>_<region>`: mittlere hintergrundkorrigierte Intensität der gezählten Zellen.
- `readout_*`: die zusätzlichen Messgrößen (Abschnitt 4, `[readouts]`); `readout_frac_positive_<channel>`: % positive Kerne.
- `channel_check`: Ergebnis der Kanalprüfungen pro Feld (getrennte Kanaldateien).
- `nuclei_per_mm2`: erst vorhanden, wenn die Pixelgröße bekannt ist.
- `stardist_scale`, `saturated_frac`, `source_file`: Skalierungsfaktor des Detektors, Anteil gesättigter
  Kernpixel und die Rohdatei(en), aus der/denen die Zeile stammt.

*Zu den QC-Stufen:* Sie sind Faustregeln, bis sie mit Handzählungen (Ground Truth) abgeglichen wurden.
Etwa 15 von Hand gezählte Bilder über verschiedene Dichten hinweg würden daraus echte Fehlerabschätzungen machen.

### 5.3 Generationen: `generations/g1`, `g2`, ...
- **g1** = vollautomatische Zählungen. Immer vorhanden.
- **g2** = g1 + die in `manual_1` eingetragenen Korrekturen. **g3** = g2 + `manual_2`, und so weiter.
- Die Kette endet bei der ersten manuellen Spalte, die vollständig leer ist.
- Passt counts.csv nicht mehr zu einer bestehenden Generation (Sie haben `manual_1` nach dem Bau von g2 geändert),
  wird der Ordner **nicht überschrieben**: `run` warnt Sie, und `--rebuild g2` aktualisiert ihn bewusst. Haben sich
  nur die Intensitätsergebnisse, das Design oder die Beschreibung des Experiments geändert (Zählungen unverändert),
  wird die Generation automatisch neu gebaut, und das Protokoll nennt den Grund.
- Eine Generation wird in `gN.tmp/` gebaut und erst eingesetzt, wenn sie vollständig ist. Dateien, die Sie selbst
  in einen Generationsordner gelegt haben, werden bei einem Neubau nach `gN/_kept/` verschoben, nie gelöscht.
- Jede Generation enthält `README.txt` (woraus sie gebaut ist), `counts_used.csv` (der für jedes Bild verwendete Wert
  und seine Herkunft: `auto`, `manual_1`, `manual_1: excluded`, ...; dazu Well, Position und die
  Experimentbeschreibung), `generation.json`, `acceptance.csv` (5.4.4), `stats_long.csv` (jedes Testergebnis in einer
  Tabelle) und `manifest.json` (eine Liste aller Diagramme und Tabellen, die der Browser liest).

### 5.4 Innerhalb einer Generation: Behandlung → n → alle n
```
g1/
├── results.xlsx                      ★ alles Folgende in einer lesbaren Arbeitsmappe
├── counts_used.csv, stats_long.csv, acceptance.csv, manifest.json
├── <n>/                              ein Versuch
│   ├── <treatment>/fields.csv        Werte pro Bildfeld (Zählung + Intensitäten + Messgrößen)
│   ├── <treatment>/count_fields.png  Balken pro Bildfeld (x = ausgeschlossen)
│   ├── count_summary.csv             pro Behandlung: Bildfelder, Mittelwert, SD, SEM, Median, % jeder Kontrolle
│   ├── count_prism.csv               Prism-Column-Tabelle: Spalten = Behandlungen, Zeilen = Bildfelder
│   ├── count_stats.csv               einfaktorielle ANOVA + Dunnett gegen jede Kontrolle
│   ├── count.png / .pdf              Balken = Mittelwert, Fehler = SD, Punkte = Bildfelder, Sterne = Dunnett gegen Referenz (5.4.2)
│   ├── count_vs_<control>.png / .pdf dasselbe, gegen jede weitere Kontrolle
│   ├── count_fields_by_well.png, count_position_effect.png   Design-Prüfungen
│   └── count_twoway.png / .pdf       zweifaktorielle ANOVA Behandlung × Well + Tukey alle Paare (5.4.3)
└── all_n/                            alle Versuche zusammen
    ├── count_summary.csv             pro Behandlung: Wert pro n, N, Mittelwert, SD, SEM, % der Kontrollen
    ├── count_prism.csv               Prism-Column-Tabelle: Spalten = Behandlungen, Zeilen = n
    ├── count_prism_pct_of_<control>.csv   dasselbe, innerhalb jedes n auf jede Kontrolle normiert
    ├── count_stats.csv               RM-ANOVA (einfaktoriell) + Dunnett gegen jede Kontrolle; gepaarte Ratio-t-Tests (Kontrolle)
    ├── count.png / .pdf              Balken = Mittelwert der n, Fehler = SEM, Punktfarbe = n (Legende), Sterne (5.4.2)
    ├── count_pct_of_<control>.png    normierte Fassung (gestrichelte Linie = 100 %)
    ├── count_model_absolute.png, count_model_means.csv      gemischtes Modell: geschätzter Mittelwert pro Behandlung
    ├── count_forest_<family>.png, count_contrasts.csv       gemischtes Modell: jeder Vergleich mit seinem KI
    ├── count_variance_components.*, count_diagnostics.png   woher die Streuung kommt; Modellprüfungen
    ├── count_loo*.csv, count_per_n_ratio*                   Robustheit: ein n weglassen; der Effekt in jedem n
    ├── count_twoway.png / .pdf       zweifaktorielle ANOVA Behandlung × n + Tukey alle Paare (5.4.3)
    ├── layers_<channel>_*            Ebenenvergleich (nur Markerkanäle, 5.4.5)
    └── n_colours.csv
```
Derselbe Satz an Dateien existiert für jede weitere Messgröße: jede Intensität (`intensity_<channel>_*`), jede
zusätzliche Messgröße (`readout_*`) und % positiv. Bei mehreren Zelllinien wird die Auswertung über alle n pro
Linie in `all_n/<cell line>/` durchgeführt. Details zu jeder Datei: [docs/statistics.md](docs/statistics.md).

**Prism:** Öffnen Sie eine `*_prism*.csv` und fügen Sie sie in eine Column-Tabelle ein, oder nutzen Sie File → Import.
Die Dateien über alle n enthalten die Bezeichnung des n als Zeilentitel in der ersten Spalte.

**`results.xlsx`** enthält ein Blatt *About* (was diese Generation ist, welcher Test die Sterne gesetzt hat,
Akzeptanz, Duplikate, Hinweise), *Counts used* sowie Blätter für Zählungen, Intensitäten und weitere Messgrößen auf
den Ebenen alle n und pro n, das Modell (Mittelwerte, Verhältnisse, Vergleichsfamilien, Varianzkomponenten),
Robustheit, zweifaktorielle ANOVA und Akzeptanz.

#### 5.4.1 Spalten von `*_stats.csv`
| Spalte | Bedeutung |
|---|---|
| `test` | Welcher Test, und mit welcher Kontrolle verglichen wird |
| `comparison` | `<treatment> vs <control>` oder `all` für eine globale ANOVA |
| `N` | (alle n) Anzahl der n, die beide Gruppen enthalten |
| `effect_type`, `effect`, `ci_low`, `ci_high` | (alle n) Dunnett: mittlere Differenz zur Kontrolle und ihr 95-%-KI (Dunnett-adjustiert); Ratio-t-Test: Fold Change als geometrisches Mittel und sein 95-%-KI |
| `statistic` | t-, F- oder Dunnett-Statistik |
| `p` / `p_adj` | rohes / für multiples Testen adjustiertes p (Dunnett ist per Konstruktion adjustiert; t-Tests nach Holm) |
| `significance` | Sterne aus `p_adj`: ns, \*, \*\*, \*\*\*, \*\*\*\* (< 0,05 / 0,01 / 0,001 / 0,0001) |

Die Ergebnisse des gemischten Modells stehen in `*_contrasts.csv` und `stats_long.csv` (gleiches Prinzip: Schätzwert,
KI, p, adjustiertes p, Sterne, dazu `estimand`, das genau angibt, was die Zahl ist).

#### 5.4.2 Wie Signifikanz in den Diagrammen dargestellt wird
- **Sterne** (\* bis \*\*\*\*) über einem Balken erscheinen nur, wenn der Vergleich dieser Behandlung signifikant ist
  (adjustiertes p < 0,05). **Keine Markierung = nicht signifikant.** Die verwendete Kontrolle ist über ihrem Balken mit **„ref“** beschriftet.
- **Diagramme pro n** verwenden Dunnett gegen die Kontrolle: `count.png` gegen `reference`, `count_vs_<X>.png`
  gegen jede weitere Kontrolle. Fehlt die Referenz in diesem n, verwendet `count.png` die nächste Kontrollgruppe,
  und die Fußzeile nennt sie. Diese Sterne beschreiben nur einen einzigen Versuch.
- **Diagramme über alle n** beziehen ihre Sterne aus `[analysis] primary_test`:
  - `mixed` (Standard): ein hierarchisches Modell auf der log-Skala mit Bildfeldern in Wells in n und einer eigenen
    n-zu-n-Varianz für jede Behandlungsklasse. Vergleiche sind Verhältnisse (0,70 = 30 % weniger Kerne) mit
    95-%-KI, adjustiert für alle Vergleiche gegen dieselbe Kontrolle. Lässt sich das Modell für eine Messgröße nicht
    anpassen, fällt diese Messgröße auf die RM-ANOVA zurück, und die Fußzeile sagt es.
  - `rm_anova`: der Prism-Standard, eine einfaktorielle ANOVA mit Messwiederholung (n = gepaarte Wiederholung), dann
    Dunnett gegen die Kontrolle mit dem gepoolten Fehler aller Behandlungen.
- Die Fußzeile jedes Diagramms nennt den Test, sodass eine Abbildung auch ohne die Tabellen lesbar ist.

**Wie die Statistik zu lesen ist:**
- **Innerhalb eines n** (`<n>/count_stats.csv`): Die Bildfelder sind *technische* Replikate, daher beschreiben diese
  p-Werte nur diesen einen Versuch. Berichten Sie sie nicht als Ergebnis.
- **Alle n** (`all_n/`): Jedes n ist ein biologisches Replikat. Mehrere Tests werden immer nebeneinander angegeben:
  - *Gemischtes Modell* (`*_contrasts.csv`, Standard für die Sterne): nutzt jedes Bildfeld und jedes Well, behandelt
    aber das n als Einheit der Replikation. Seine Freiheitsgrade ergeben sich aus der Anzahl der n, daher bleiben die
    p-Werte auch mit 3–6 n ehrlich. Es nimmt nicht an, dass jede Behandlung gleich stark zwischen den n schwankt.
  - *RM-ANOVA (einfaktoriell) + Dunnett* (`*_stats.csv`): der Prism-Standard. Er poolt einen Fehler für alle
    Behandlungen, was viele Freiheitsgrade ergibt, aber zu optimistisch ist, wenn die Effekte mancher Behandlungen
    zwischen den n viel stärker schwanken als die Kontrollen. Ein n, dem eine Gruppe fehlt, zählt für seine übrigen
    Gruppen trotzdem.
  - *Gepaarter Ratio-t-Test* (Kontrolle, ohne Sterne): für jede Behandlung das Verhältnis zur Kontrolle innerhalb
    jedes n, über die n getestet, Holm-korrigiert. Der vorsichtigste Test; mit 3 n erreicht er selten p < 0,05.
  - *Vehikelprüfung* (Äquivalenztest): „nicht signifikant“ heißt nicht „kein Vehikeleffekt“. Nur **equivalent**
    stützt das; **inconclusive** bedeutet, dass mehr n nötig sind.
  - *Robustheit*: Weglassen je eines n (`stable = yes`, wenn kein einzelnes n die Schlussfolgerung ändert) und der
    Effekt in jedem n.
  - Ein Test verwendet nur die n, die *beide* Gruppen enthalten; daher kann `N` zwischen Vergleichen variieren, wenn
    eine Gruppe in manchen n fehlt. Mit einem einzigen n gibt es Schätzwerte, aber keinen Test über die n
    („needs ≥ 2 n“).
- **Prozentwerte** (% positiv, `readout_frac_*`): Lesen Sie die Differenz in Prozentpunkten und die absoluten
  Prozentwerte; Verhältnisse von Prozentwerten sind zweitrangig.
- **Normierte Werte** (`pct_of_...`) sind für den Vergleich über die n meist besser geeignet, da sie Schwankungen
  von Tag zu Tag bei Aussaat, Färbung und Lampenhelligkeit herausrechnen.
- **Intensitätswerte** sind nur vergleichbar, wenn Belichtung und Verstärkung innerhalb jedes n fest eingestellt waren (Auto-Belichtung aus).

#### 5.4.3 Zweifaktorielle ANOVA mit Vergleich aller Paare (`*_twoway*`)
Eine zusätzliche Auswertung auf beiden Ebenen, wie Prisms zweifaktorielle ANOVA mit „compare all pairs“:
- **Pro n:** Behandlung × Well (Wells aus `[design]`); die Bildfelder sind die Replikate.
- **Alle n:** Behandlung × n; die Bildfelder sind die Replikate.

Das Modell lautet `value ~ treatment + zweiter Faktor + Interaktion` mit Quadratsummen vom Typ II, sodass eine Gruppe,
die in einem n fehlt, korrekt behandelt wird. `*_twoway_anova.csv` liefert F, df, p und den
Anteil an der Gesamtvariation in % für Behandlung, den zweiten Faktor und die Interaktion. `*_twoway_pairs.csv` vergleicht
jedes Behandlungspaar (Tukey–Kramer, `p_adj`) anhand der Zellmittelwerte, gepaart innerhalb jedes Wells / n, und verwendet dabei nur die
Wells / n, die beide Gruppen enthalten (`matched_levels`), mit der mittleren Differenz und ihrem 95-%-KI. Das Diagramm
zeigt gruppierte Balken (Mittelwert ± SEM der Bildfelder pro Well / n) neben einer Matrix der adjustierten p-Werte.

**Vorsicht auf der Ebene aller n:** Werden Bildfelder als Replikate behandelt, fallen diese p-Werte optimistischer aus
als bei den Tests, die das n als Einheit verwenden. Berichten Sie den Haupttest als Ergebnis; nutzen Sie die
zweifaktorielle Auswertung, um zu sehen, ob sich die Behandlungen untereinander unterscheiden und ob der Effekt vom n
abhängt (die Interaktion).

#### 5.4.4 Akzeptanz und Datenprüfungen
`acceptance.csv` führt für jedes n jedes Kriterium aus `[acceptance]` mit Wert, Grenze und bestanden / nicht
bestanden auf. `design_check.csv` (Projektordner) hat eine Zeile pro n × Behandlung mit der Anzahl der Bilder, Wells
und Bildfelder, fehlenden Bildnummern und Hinweisen in Klartext; eine Zeile mit 0 Bildern ist eine Bedingung, die das
n nicht hat (zum Beispiel eine fehlende Referenz). `duplicates.csv` listet Rohdateien mit identischem Inhalt: Die
Ergebnisse der betroffenen Bedingungen sind nicht interpretierbar, bis die echten Bilder gefunden sind. All das
erscheint auch im Feld **Attention**. Ungewöhnliche Designs (ein Well, eine fehlende Referenz, ein leerer Ordner, ein
einzelnes n, ...) sind in [docs/design.md](docs/design.md) beschrieben.

#### 5.4.5 Markerkanäle: `intensity/` und der Ebenenvergleich
`intensity/` (bei jedem Lauf neu geschrieben) enthält die pro n und Kanal verwendete Schwelle (`cutoffs.csv`),
% positiv pro Bedingung mit und ohne Flat-Field-Korrektur, eine Variation der Schwelle, Verteilungen pro Kern, die
Flat-Field-Karten und die Kanalprüfungen. In jeder Generation erhält jeder Markerkanal einen **Ebenenvergleich**
(`layers_<channel>_*`): derselbe Vergleich auf fünf Verarbeitungsstufen (1 ganzes Feld roh, 2 ganzes Feld korrigiert,
3 pro Kern roh, 4 pro Kern korrigiert, 5 % positiv) und eine Tabelle, die den Schritt nennt, an dem sich die
Schlussfolgerung ändert. Details: [docs/intensity.md](docs/intensity.md).

### 5.5 Wie viele Bilder werden benötigt, und Planung: `sampling/g1/`
| Datei | Inhalt |
|---|---|
| `REPORT.txt` | Zusammenfassung in Klartext, beginnend mit **Empfehlungen für das nächste Experiment**: hier anfangen |
| `images_needed.png/.pdf` | Vorhergesagtes Rauschen, Übereinstimmung beim Subsampling und Fehler pro n in Abhängigkeit von der Bildzahl pro Bedingung (gegen die Hauptreferenz; `images_needed_vs_<group>` für die anderen) |
| `batch_by_n.csv` | Pro n: Werte der Referenzgruppen, fehlende Referenzen, systematischer Unterschied zwischen Wells und ob die Behandlungseffekte dieses n stärker/schwächer sind als in den anderen n |
| `batch_by_condition.csv` | Pro n × Behandlung: Mittelwert pro Well, Verhältnis zwischen den Wells mit Markierung, wenn es das Rauschen zwischen den Bildfeldern übersteigt, Verhältnis zu jeder Referenz |
| `variance_components.csv` | Variation von Bildfeld zu Bildfeld (innerhalb eines Wells) und von Well zu Well (innerhalb einer Platte), als CV² |
| `effects.csv` | Pro Referenz und Behandlung: beobachteter Fold Change, p mit allen Bildern, biologische Varianz (von n zu n) |
| `design_power.csv`, `summary_by_design.csv` | Für jedes Design (Wells × Bildfelder), pro Referenz und Behandlung: Rauschen im Vergleich zu jetzt, Power bei der aktuellen Anzahl n und die für 80 % Power benötigte Anzahl n |
| `subsampling.csv` | Echte Bilder, pro Design vielfach neu gezogen: wie oft das Ergebnis mit der Schlussfolgerung aus allen Bildern übereinstimmt und wie stark sich das Verhältnis jedes n verschiebt |
| `planner.json`, `planner_components.csv` | Der vorausschauende Planer (Streuung zwischen Bildfeldern, Positionen, Wells und n mit Intervallen), genutzt von der Seite **Planner** |
| `design_forward.csv`, `mde.csv`, `cost_frontier.csv`, `allocation.csv`, `more_n.csv` | Power und KI-Breite für jedes Design, der kleinste nachweisbare Effekt, die günstigsten Designs, die die Ziel-Power erreichen, die Aufteilung der Wells zwischen Behandlungen und Referenzen und die Wahrscheinlichkeit, dass weitere n einen Vergleich klären (+ Diagramme) |

Jede Auswertung wird einmal pro Referenzgruppe durchgeführt (`reference` + `normalise_to`), wie die Statistik in
Abschnitt 5.4: Ein n ohne eine der Kontrollen wird bei den Vergleichen gegen diese Kontrolle ausgelassen, zählt aber
gegen die anderen weiterhin. Vorhersagen werden pessimistisch / erwartet / optimistisch angegeben, weil mit 3–4 n jede
Schätzung der Streuung unsicher ist. Zelllinien werden getrennt geplant (`sampling/g1/<cell line>/`). Details:
[docs/planner.md](docs/planner.md).

**Wie es zu lesen ist:** Der Test über alle n wird durch zwei Dinge begrenzt. Das erste ist das *Messrauschen*
(aus Bildfeldern und Wells), das mit mehr Bildern abnimmt. Das zweite ist die *Variation von n zu n*,
die nur mit mehr n abnimmt. Mehr Bilder helfen nur, solange das Messrauschen noch vergleichbar groß ist wie die
Variation von n zu n. Das Rauschen zwischen den Wells nimmt nur mit mehr Wells ab.

**Batch-Effekte:** Jedes n ist eine Platte an einem Tag, daher sind Platten- und Tageseffekte Teil der Variation
von n zu n und lassen sich nicht davon trennen. Effekte, die eine ganze Platte skalieren (Aussaat, Färbung, Lampe), heben sich
im Verhältnis zur Referenz auf. Effekte, die verändern, wie stark die Zellen reagieren (z. B. Transfektionseffizienz),
heben sich nicht auf; der Bericht markiert ein n, dessen Behandlungseffekte durchgehend stärker oder schwächer sind
als in den anderen. Eine Markierung ist ein Anlass, die Laborunterlagen dieses Tages zu prüfen; sie ist kein Grund, das n auszuschließen.

**Typisches Muster:** Wenn Zellen verklumpen, ist die Variation zwischen den Bildfeldern groß. Eine Aufnahme mit geringerer
Vergrößerung (10x deckt etwa die 16-fache Fläche von 40x ab) verringert das Messrauschen dann weit stärker als
zusätzliche Bildfelder bei hoher Vergrößerung. Mit nur wenigen n sind die Ergebnisse als Orientierung für die Planung
zu verstehen, nicht als exakte Werte.

---

## 6. Typische Arbeitsabläufe
**Ein neues n für ein bestehendes Experiment**
1. Legen Sie den neuen n-Ordner in `data/<experiment>/` (gleiche Struktur wie die anderen).
2. `kerncount init <experiment>`. Ihre `layout.csv` bleibt erhalten, und ein neuer Scan wird in `layout.new.csv`
   geschrieben; übernehmen Sie die neuen Zeilen in `layout.csv` (Namen prüfen), oder verwenden Sie den neuen Scan auf
   der Seite **Layout**.
3. `kerncount run <experiment>`. Nur die neuen Bilder werden verarbeitet; der Rest kommt aus dem Cache.
   Generationen mit manuellen Korrekturen werden nicht überschrieben. Bauen Sie sie bewusst mit `--rebuild` neu, sobald
   das neue n geprüft ist.

**Ein neues Experiment**
1. Legen Sie `data/<neues Experiment>/` mit einem Ordner pro n an.
2. `kerncount init <neues Experiment>`, `layout.csv` und die Einstellungen unter `[experiment]` / `[analysis]`
   prüfen, dann `kerncount run <neues Experiment>`.

**Danach, in beiden Fällen:** Prüfen Sie das Feld **Attention** und die roten und orangen Kacheln in
`review/*_overview.jpg`, korrigieren Sie Zählungen in `counts.csv` → `manual_1` und führen Sie die Analyse erneut aus.
Berichten Sie aus der neuesten Generation: `generations/g<N>/results.xlsx`, die Diagramme in `all_n/` und die Prism-CSV-Dateien.

**Das nächste Experiment planen:** `kerncount sampling <experiment>`, dann die Seite **Planner** (oder
`sampling/g<N>/REPORT.txt`), dann `kerncount platemap <experiment>` für die Plattenbelegungen.

## 7. Fehlerbehebung
| Problem | Lösung |
|---|---|
| „project not found“ / „raw folder not found“ | Führen Sie den Befehl aus dem Arbeitsbereich-Ordner aus (dem Ordner mit `data/`), prüfen Sie den Experimentnamen und führen Sie zuerst `kerncount init <experiment>` aus. |
| Der erste Lauf scheitert beim Herunterladen des Modells | Das StarDist-Modell wird einmalig bei der ersten Verwendung heruntergeladen; stellen Sie eine Internetverbindung her und starten Sie erneut. |
| TensorFlow lädt nicht (Windows) | Installieren Sie das Microsoft Visual C++ Redistributable, oder verwenden Sie `[segmentation] backend = "classical"`. |
| `counts.csv` öffnet sich in Excel als eine einzige Spalte | Setzen Sie `csv_delimiter = ";"` in `settings.toml` (europäisches Excel). Dateien, die Excel mit `;` neu speichert, werden in jedem Fall korrekt gelesen. |
| `run` kann `counts.csv` nicht schreiben, oder „cannot replace gN“ (Windows) | Schließen Sie die Datei zuerst in Excel; Excel sperrt geöffnete Dateien. Die neuen Ergebnisse warten in `gN.tmp`, bis Sie erneut starten. |
| „cannot read manual value“ | Eine manuelle Zelle enthält etwas anderes als eine ganze Zahl oder `x`; die Meldung nennt Zeile und Spalte. |
| „g2 already exists but counts.csv now gives different values“ | Sie haben eine manuelle Spalte geändert, nachdem g2 gebaut wurde. Führen Sie die Analyse mit `--rebuild g2` erneut aus, um sie bewusst zu aktualisieren. |
| Ein Behandlungsname ist falsch oder ein Ordner soll ausgeschlossen werden | Bearbeiten Sie `layout.csv` (`treatment`, `include`) und führen Sie die Analyse erneut aus. |
| „has 1 channel(s) but settings ask for channel 2“ | `[channels] segment` verweist auf einen Kanal, den die Bilder nicht haben. (Kanäle in `[intensity] measure`, die nicht existieren, werden stillschweigend übersprungen; prüfen Sie die Spalten `intensity_*`.) |
| Marker-Messgrößen für einen ganzen Ordner leer | Die Kanalprüfung ist gescheitert (Dateien nicht so gepaart, wie `file_order` sagt, oder eine fehlende Datei verschiebt jedes folgende Paar). Siehe `intensity/alignment.csv` und das Feld Attention. |
| Sterne fehlen, „needs ≥ 2 n“ | Ein Test über die n braucht mindestens 2 n, die beide Gruppen enthalten. |
| Der Lauf ist langsam | Erste Läufe brauchen ca. 2–3 s pro Bild; spätere Läufe nutzen den Cache. `--no-review` überspringt die Kontrollbilder. |
| Zählungen wirken in einem ganzen Ordner falsch | Prüfen Sie `diameter_px` in `qc_details.csv`; liegt die Größenschätzung daneben (z. B. bei einem fast leeren Ordner), setzen Sie `nucleus_diameter_px` auf eine Zahl. |
| Tumor-Lauf „blockiert“ | Siehe die Liste oben auf der Seite **Layout** (oder die Ausgabe des Laufs): die Zelle korrigieren, oder *include* abwählen und den Grund in *note* angeben (Abschnitt 9.2). |
| Tumor-Lauf verweigert einen geänderten Analyseplan | Der Plan ist gesperrt, sobald Ergebnisse existieren; mit Begründung entsperren (Abschnitt 9.3). |

## 8. Hinweise und bekannte Einschränkungen
- 8-Bit-Kamerabilder, die als 16-Bit-TIFF gespeichert wurden (jeder Wert × 257), werden korrekt verarbeitet. Sättigung wird in beiden Fällen erkannt.
- Im Vergleich zu einem einfachen Fiji-Schwellenwert-Makro fallen die Zählungen meist niedriger aus. Schwellenwert-Makros zählen
  Hintergrundrauschen in leeren oder unscharfen Bildfeldern mit, zerteilen gelappte Kerne und zählen jeden Kern, der den Rand berührt.
- Sehr dichte Bildfelder, in denen die Zellen helles Zytoplasma haben, sind der schwierigste Fall. Sie werden mit QC niedrig/mittel markiert;
  prüfen Sie sie in den Kontrollbildern.
- Ohne bekannte Pixelgröße werden die Ergebnisse pro Bildfeld statt pro mm² angegeben.
- Die QC-Stufen sind Heuristiken, bis sie mit Handzählungen kalibriert wurden.
- Mit 3–4 n ist jede Schätzung zwischen den n (Varianzen, Power) grob; lesen Sie die Konfidenzintervalle, nicht nur
  die Sterne.
- Die eingebauten Regeln für Behandlungsnamen und die Standardnamen der Kontrollen (`Control`, `LyoVec`, `CL075`)
  stammen aus dem Labor, für das das Werkzeug geschrieben wurde; andere Namen funktionieren, setzen Sie aber
  `reference`, `normalise_to` und `treatment_order` auf Ihre eigenen.
- Algorithmen und Implementierungshinweise: [docs/technical.md](docs/technical.md).

## 9. Tumorexperimente (fluoreszierende Tumoren in Hirnschnitten)

Dasselbe Werkzeug wertet auch **Tumorexperimente** aus: fluoreszenzmarkierte (z. B. mCherry) Tumoren in Hirnschnitten,
aufgenommen als konfokale z-Stapel (`[experiment] type = "tumour_volume"`). Ein **Tumor** (Hemisphäre A oder B eines
Schnitts) ist ein n; A und B werden nie gemittelt. Technische Details: [docs/tumour.md](docs/tumour.md).

### 9.1 Eingabe
Ein Ordner pro Experiment in `data/`, mit **einem Unterordner pro Aufnahme-Sitzungsordner wie geliefert** (z. B.
`5d KO`, `5d WT`), die Dateien wie geliefert:
```
data/my_tumours/
├── 5d KO/                                      ← Sitzungsordner: sein Name ergibt den Genotyp
│   ├── 010126 lineA +drug n3 A.ome.tiff        ← ein originaler z-Stapel = ein Tumor
│   ├── 010126 lineA +drug n3 Af.tif            ← Fiji-Maske (Altbestand, optional; nur zum Vergleich gelesen)
│   └── Statistics for 010126 lineA +drug n3 Af.csv   ← Tabelle des Fiji 3D Objects Counter (Altbestand, optional)
└── 5d WT/
    └── 150126 control -drug n1 B.ome.tif
```
- **Faktoren.** Welche Gruppen es gibt, wird in `settings.toml` als `[[tumour.factors]]` festgelegt: Jeder Faktor hat
  seine Stufen, eine Referenzstufe, andere Schreibweisen (`aliases`), seinen Platz in den Diagrammen und seine
  Quelle: den **Namen des Sitzungsordners** oder ein Wort im **Dateinamen**. Die Einstellungsvorlage enthält ein
  Beispieldesign mit drei Faktoren: einen Genotyp (aus dem Namen des Sitzungsordners) und zwei Faktoren aus dem
  Dateinamen (eine Tumorlinie und eine ±-Behandlung). Passen Sie die Stufen an Ihr Experiment an; eine neue Stufe oder
  ein neuer Faktor ist ein Eintrag in den Einstellungen, keine Codeänderung. Das Beispiel oben nimmt Stufen namens
  `lineA` / `control` und `+drug` / `-drug` an.
- **Sitzungsordner** werden Wort für Wort gelesen (Groß-/Kleinschreibung egal); Wörter, die keine Stufe sind (`5d`),
  bleiben als Sitzungsbezeichnung erhalten. Beliebig viele Sitzungsordner sind möglich, auch mehrere eines Genotyps
  (Wiederholungssitzungen). Ein Ordner, dessen Genotyp nicht lesbar ist, stoppt den Lauf mit einer Meldung, die sagt,
  wie es sich beheben lässt.
- **Unterordner spielen keine Rolle**: Unterhalb eines Sitzungsordners dürfen die Dateien in beliebiger
  Ordnerstruktur liegen; alles außer den Faktoren aus dem Sitzungsordner wird aus dem Dateinamen gelesen. Endungen
  `_s1`, `1n` für `n1` und ein fehlendes Datum werden verarbeitet; eine `.czi` neben der `.ome.tiff` desselben Tumors
  wird vermerkt, ein Tumor, der nur als `.czi` vorliegt, wird gemeldet (als OME-TIFF exportieren). Ein Sitzungsordner
  darf ein Link auf ein externes Laufwerk sein; ist das Laufwerk nicht angeschlossen, meldet der Lauf das, und
  vorhandene Ergebnisse bleiben einsehbar.
- **Umbenennen eines Sitzungsordners** verliert nichts: `kerncount run` findet die Dateien über ihren Namen, behält
  Ihre Änderungen in layout.csv und verwendet jede zwischengespeicherte Messung wieder; `kerncount init` aktualisiert
  danach layout.csv (die alte Fassung wird gesichert).
- **Dateinamen** liefern Datum, die Faktoren aus dem Dateinamen, Schnitt und Hemisphäre. Für **neue Aufnahmen**
  verwenden Sie `DDMMYY <Faktorstufen> n<Schnitt> <A|B>.ome.tif`: Datum als 6 Ziffern, jede Faktorstufe genau wie in
  den Einstellungen, die Schnittbezeichnung `n<k>` (beginnt pro Gruppe neu; `s3` wird nicht als Schnitt gelesen),
  Hemisphäre `A` (oben) oder `B` (unten). Unterstriche statt Leerzeichen, `.ome.tiff` und Aliase (z. B. bekannte
  Tippfehler) werden ebenfalls gelesen und in `layout.csv` vermerkt.
- **Batches.** Was als Batch gilt, legt die Einstellung `[tumour] batch` fest: standardmäßig der **Kultur-Batch** (das
  6-stellige Datum am Anfang der Dateinamen = die Kulturpräparation), alternativ der *Sitzungsordner* oder der
  *Aufnahmetag* (aus den Metadaten jedes Stapels). Mehrere Batches, die eine Gruppe gemeinsam haben, gehen als fester
  Blockeffekt ins Modell ein. Ein Vergleich, der sich nicht vom Batch trennen lässt (z. B. zwei Genotypen aus
  verschiedenen Batches), wird gezeigt, als „confounded“ beschriftet, aber nie getestet; Verhältnisse innerhalb eines
  Genotyps lassen sich trotzdem zwischen Genotypen vergleichen. Der Bericht nennt die verwendete Definition.
- Voxelgröße und Aufnahmeeinstellungen werden aus den OME-Metadaten jeder Datei gelesen (nie von Hand eingegeben);
  die z-Schritte dürfen sich zwischen Dateien unterscheiden.

### 9.2 Einrichten, Layout prüfen, Analyse starten
Wählen Sie in der Weboberfläche das Experiment, unter *Set up* den **Experimenttyp: Tumour volume (z-stacks)** und
drücken Sie **Set up** (Kommandozeile: `kerncount init my_tumours --type tumour_volume`). Dann:
- **Layout**: eine Zeile pro Tumor mit seinen Faktorwerten, Schnitt, Hemisphäre, Tumor-ID (Faktorstufen + Schnitt und
  Hemisphäre, z. B. `KO / lineA / +drug / n3A`), include und note. Probleme, die eine Datei zu einem falschen
  Datenpunkt machen würden (ein unlesbarer Name, eine doppelte Tumor-ID, eine Sitzung ohne Genotyp, eine fehlende
  Datei, ein Stapel im Rohdatenordner, den layout.csv noch nicht aufführt, ein Ausschluss ohne Grund), stehen oben und
  **blockieren den Lauf**, bis Sie die Zelle korrigieren oder *include* abwählen und den Grund in *note* angeben.
  Später hinzugekommene Stapel: **Set up** erneut ausführen (`kerncount init my_tumours`); das ergänzt sie und behält
  Ihre Änderungen.
- **Settings** → *Analysis plan*: die Faktoren (Stufen, Referenzstufe, Platz in den Diagrammen), welche Messgrößen
  konfirmatorisch sind, und die Kontrastfamilien (unten). *Tumour measurement* (`[tumour]`): Glättung,
  Schwellenwertregel, minimale Objektgröße, Dicke des Mittelschnitts, Glättung des Umrisses, Schwellen-Sweep,
  Batch-Definition, QC-, Review- und Planungshilfe-Einstellungen, jeweils mit Erklärung.
- **Run** (`kerncount run my_tumours`): Optionen *skip review images* (`--no-review`), *segment every stack again*
  (`--resegment`), *threshold sweep*. Der erste Lauf dauert etwa 20–60 s pro Stapel, plus etwa noch einmal so viel
  für den **Schwellen-Sweep** (standardmäßig an: jeder konfirmatorische Effekt wird bei Schwelle × 0,8–1,2 neu
  berechnet; abwählen oder `threshold_sweep = false` setzen, um diese Zeit zu sparen). Segmentierungen werden
  zwischengespeichert, spätere Läufe sind daher schnell.

### 9.3 Messgrößen, und warum Fläche statt Volumen
Oft wird nur eine Scheibe (z. B. ca. 80–150 µm) eines viel größeren Tumors aufgenommen; ein 3D-Volumen misst dann vor
allem, wie tief der Stapel aufgenommen wurde. Die **konfirmatorischen Messgrößen** sind daher 2D, auf einem
**Mittelschnitt** fester Dicke (`mid_section_um`, Standard 20 µm, um die Ebene, in der die Hauptmasse des Tumors am
größten ist; für jeden Tumor exakt dieselbe physikalische Dicke, unabhängig vom z-Schritt; die verwendete Dicke steht
in `tumours.csv`): die **Querschnittsfläche** (mm²) und der **Umfang** (mm) der Hauptmasse (des größten
zusammenhängenden Objekts). Alles andere ist **explorativ** (Effektgrößen mit KI, beschriftet, ohne Sterne):
Grundfläche, Scheibenvolumen, Volumen und Oberfläche der Hauptmasse, Tumorlast (Hauptmasse + Satelliten),
Elongation, Solidität, Randrauigkeit, Zirkularität und die Fiji-Altwerte. Eine einzige automatische
Schwellenwertregel für jeden Tumor, nie von Hand pro Tumor angepasst.

**Sperre des Analyseplans.** Faktoren, Rollen der Messgrößen, Kontrastfamilien, die Batch-Regel und die
Messeinstellungen hinter den konfirmatorischen Messgrößen (Glättung, Schwellenwertregel, minimale Objektgröße, Dicke
des Mittelschnitts, Glättung des Umrisses) werden festgelegt, bevor Ergebnisse angesehen werden. Sobald Ergebnisse
existieren, ist der Plan **gesperrt**; eine Änderung erfordert **Unlock** mit einer Begründung — in der
Weboberfläche oder auf der Kommandozeile `kerncount run my_tumours --unlock "Grund"` (`kerncount run` verweigert
einen geänderten Plan ohne diese Option). Das Entsperren und jede Änderung werden in `lock_log.csv` geschrieben, und
der nächste Lauf gibt den Verlauf (alte → neue Plan-Sperre) im Bericht und im About-Blatt von `results.xlsx` aus.
Jede Statistikzeile trägt die Sperre, mit der sie erstellt wurde. QC-, Sweep-, Sensitivitäts-, Review- und
Planungshilfe-Einstellungen gehören nicht zum Plan und können jederzeit geändert werden.

### 9.4 Ausgabe (`analysis/<experiment>/`)
| Datei / Ordner | Inhalt |
|---|---|
| `results.xlsx` | alles in einer Arbeitsmappe: Tumoren, Statistik, geschätzte Mittelwerte, Familien, QC, Fluss, Sensitivitätsanalysen, Schwellen-Sweep, Übereinstimmung mit Fiji, Planungshilfe (Download auf der Seite **Results**) |
| `stats_report.md` | der Bericht: Analyseplan und seine Sperre (mit Verlauf), Ergebnisse pro Messgröße und Familie, QC, Aufnahmeeinstellungen, Sensitivität, Übereinstimmung mit Fiji, Planungshilfe |
| `tumours.csv` | eine Zeile pro ausgewertetem Tumor: Faktorwerte, Batch, jede Messgröße, QC-Markierungen |
| `excluded_tumours.csv`, `flow.csv` | ausgeschlossene Tumoren mit Grund (Layout oder Review); Dateien → Tumoren → ausgewertet, pro Gruppe |
| `stats_long.csv` | jeder Schätzwert in einer langen Tabelle (die Seite **Statistics**; Sensitivitätszeilen nur auf Wunsch sichtbar) |
| `acquisition.csv` | pro Tumor die Mikroskopeinstellungen aus den Dateimetadaten (Detektorverstärkung, Laser, Pinhole, Voxelgröße, Aufnahmetag, ...) |
| `concordance.csv`, `threshold_sweep.csv` | neue Methode vs. Fiji pro Kontrast; jeder konfirmatorische Effekt bei jedem Schwellenfaktor |
| `lock_log.csv`, `corrections_used.csv` | Entsperrungen und Änderungen des Plans; die angewendeten Review-Korrekturen |
| `plots/` | pro Messgröße: absolut, normiert (% der Referenzstufe innerhalb der übrigen Faktoren), Effekte (Forest), Interaktion, A vs. B, Residuen, und mit Fiji-Daten / Sweep: vs. Fiji, Schwellen-Sweep; `shape_panel` |
| `prism/` | Prism-fertige Tabellen, eine Spalte pro Gruppe (absolut und normiert) |
| `review/` | die **verblindeten** Review-Bilder und `corrections.csv` (unten) |
| `review_key.csv` | Code → Tumor (außerhalb von `review/`, damit der Review-Ordner verblindet bleibt) |

**Statistik:** ein lineares Modell auf log(Messgröße) über alle Faktoren; Effekte sind Fold Changes mit 95-%-KI,
Holm innerhalb jeder Familie. Nur konfirmatorische Familien konfirmatorischer Messgrößen erhalten Sterne.
**Sensitivitätsanalysen** (Welch-artig, Permutation, ohne QC-markierte Tumoren) werden daneben berichtet, nie
stattdessen. Die **Planungshilfe** nennt den kleinsten nachweisbaren Fold Change bei den aktuellen Gruppengrößen und
die Tumoren pro Gruppe, die für `planner_fold_change` (Standard 1,5×) nötig sind; sie fließt nie in die Ergebnisse
ein. Die **Aufnahmeprüfung** sagt, ob jede helligkeitsrelevante Mikroskopeinstellung identisch war und welche Gruppen
an nur einem Tag aufgenommen wurden; mit `shared_threshold_check = true` (und identischen Einstellungen) wird die
Auswertung bei einer gemeinsamen Schwelle wiederholt. Der Bericht weist außerdem darauf hin, dass das Fiji-Altvolumen
vor allem die aufgenommene Stapeltiefe widerspiegelt, und gibt die Schnittdicke sowie an, wie der z-Bereich am
Mikroskop gewählt wurde (Einstellungen `section_thickness_um`, `z_range_rule`; „not recorded“, solange sie leer
sind).

**QC:** Ein Tumor wird **markiert**, wenn seine Hauptmasse den Bildrand berührt oder sein Mittelschnitt nicht in den
Stapel passt; eine **Sättigungswarnung** gibt es oberhalb von `qc_saturation_warn`. Markierte Tumoren bleiben in der
Hauptauswertung (hohl gezeichnet) und werden in einer Sensitivitätsanalyse weggelassen.

### 9.5 Review und Korrekturen (verblindet)
`review/` enthält ein Bild pro Tumor, benannt nach einem zufälligen **Code** (nirgends Gruppe, Sitzung, Datum oder
Dateiname): links die Maximumprojektion, rechts die Ebene des Mittelschnitts; die Hauptmasse orange, Satelliten cyan,
Objekte nummeriert. Die Seiten **Review images** und **Corrections** zeigen nur die Codes. In einer Runde
(`manual_1` → g2, `manual_2` → g3, …) tragen Sie ein:
- eine **Objektnummer** aus dem Bild, um dieses Objekt zur Hauptmasse zu machen (`1` = zurück zur Automatik);
- **`x: Grund`**, um den Tumor auszuschließen (der Grund kommt in `excluded_tumours.csv`);
- nichts = keine Änderung.

Danach die Analyse erneut ausführen: Nichts wird neu segmentiert; die Messgrößen des gewählten Objekts werden neu
berechnet (auch der Schwellen-Sweep folgt dem gewählten Objekt). Sie können `review/corrections.csv` auch in Excel
bearbeiten und nur den Ordner `review/` an eine Kollegin oder einen Kollegen zur verblindeten Prüfung weitergeben.

## 10. Entwicklung
```
pip install -e ".[test]"   # Python 3.11 oder neuer
pytest                     # schnelle Testsuite
pytest -m slow             # aufwendige Simulationen und Durchläufe von Anfang bis Ende (standardmäßig übersprungen)
```
`pytest -m "slow or not slow"` führt alles aus. Die Weboberfläche ist reines HTML/JavaScript in
`src/kerncount/web/` (ohne Build-Schritt); die Ein-Klick-Programme werden aus `launcher/` gebaut (Go; siehe
[launcher/README.md](launcher/README.md)).

Ausführliche Dokumentation (auf Englisch):

| Dokument | Inhalt |
|---|---|
| [docs/technical.md](docs/technical.md) | Wie das Zählen funktioniert, Algorithmen, Experimenttypen |
| [docs/statistics.md](docs/statistics.md) | Jede statistische Ausgabe, das gemischte Modell, `stats_long.csv`, `manifest.json`, wann man was vertrauen kann |
| [docs/intensity.md](docs/intensity.md) | Markerkanäle: getrennte Kanaldateien, Messung pro Kern, Flat Field, Positivität, Ebenenvergleich |
| [docs/readouts.md](docs/readouts.md) | Die zusätzlichen Messgrößen pro Bild |
| [docs/design.md](docs/design.md) | Wells, Feldpositionen, `wells.csv`, Zeitpunkte, ungewöhnliche Designs |
| [docs/planner.md](docs/planner.md) | Der Planer und die Plattenbelegungen |
| [docs/ui.md](docs/ui.md) | Die Weboberfläche |
| [docs/tumour.md](docs/tumour.md) | Tumorexperimente im Detail |

Einige dieser Dokumente verweisen auf ein Glossar (`CONTEXT.md`) und Designentscheidungen (ADRs), die außerhalb
dieses Repositorys liegen; die Dokumente sind auch ohne sie vollständig.

## Danksagung
Die Kernerkennung verwendet [StarDist](https://github.com/stardist/stardist) (Schmidt et al., MICCAI 2018;
Weigert et al., WACV 2020) mit dessen vortrainiertem Modell `2D_versatile_fluo`.

## Lizenz
MIT: siehe [LICENSE](LICENSE).
