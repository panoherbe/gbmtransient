# gbmtransient

**Ein Programm findet kurze Gammablitze in den offenen Messdaten des NASA-Satelliten Fermi, wählt seine
Einstellungen per Optimierung und belegt seine Trefferquote am offiziellen Katalog: 49 von
54 Gammablitzen gefunden.** Es läuft heute: ausgewertet am 03.11.2026 auf 61 Tagen Messdaten
von zwölf Detektoren, 24 Tests bestanden unter Linux und Windows.

*Prototype: finds short gamma-ray bursts in open NASA Fermi GBM data with an optimised, statistically calibrated
trigger and measures its hit rate against the official catalogue. Documentation in German. Not affiliated with NASA.*

## Rundgang (69 Sekunden, ohne Ton, mit Untertiteln)

![Rundgang: Tests, Auswertung von 61 Tagen, ein echter Gammablitz, gemessene Empfindlichkeit, Abgleich mit dem Katalog, Grenzen](medien/rundgang.gif)

[Video herunterladen (MP4, 1,0 MB)](medien/rundgang.mp4) · [Standbild](medien/poster.png)

<details><summary>Untertitel als Text</summary>

- **0:00** Prototyp · läuft heute. gbmtransient findet kurze Gammablitze in den offenen Messdaten des NASA-Satelliten Fermi und misst, wie zuverlässig das gelingt: mit optimierten Einstellungen und gemessener Fehlalarmrate.
- **0:05** Echte, offene Messdaten: 61 Tage (01.09. bis 31.10.2025) von zwölf Detektoren, 793 Dateien mit 9,3 GB, jede unverändert abgelegt; 182 Ereignisse im offiziellen Katalog, davon 54 Gammablitze, als Maßstab; 21.805 künstliche Ausbrüche mit bekannter Stärke.
- **0:11** Erst die Tests: 24 von 24 bestanden, ebenso unter Windows in der automatischen Prüfung.
- **0:14** Dann ein Befehl für alle 61 Tage, nur aus der unveränderten Rohablage.
- **0:19** Zur Kontrolle ein zweiter Lauf mit anderer Aufteilung auf die Prozessoren: alle Ergebnisse identisch.
- **0:23** Der Bericht: Jede Zahl stammt aus der Ergebnisdatei des Laufs; Grenzen und Fehlschläge stehen offen darin.
- **0:28** Ein echter Gammablitz vom 22.10.2025: In mehreren Detektoren zugleich steigen die Zählraten über den Untergrund; die Prüfgröße überschreitet die Schwelle deutlich.
- **0:34** Künstliche Ausbrüche bekannter Stärke wiedergesucht: Die optimierte Einstellung findet die Hälfte schon ab 87 statt 108 Zählungen pro Sekunde (kurze Pulse), ab 50 statt 67 (lange).
- **0:41** 162 Einstellungen verglichen, gewählt auf 31 Tagen, geprüft auf 30 anderen: Die Rangfolge hält auf ungesehenen Tagen.
- **0:47** Das Modell der Detektorausrichtung passt zu den gemessenen Mustern von 54 Gammablitzen (Korrelation im Median 0,92).
- **0:52** Was nachgewiesen ist: 49 von 54 Gammablitzen gefunden; 61,2 % der künstlichen Ausbrüche erkannt statt 55,6 %; 24 Tests unter Linux und Windows; zweiter Lauf identisch.
- **0:58** Ehrlich gesagt: 0,33 statt 0,1 zufällige Fehlalarme je Tag an den Prüftagen; 5 Gammablitze verpasst, meist sehr kurze; 80 Meldungen ohne Katalogeintrag, nicht geprüft; Detektormodell vereinfacht.
- **1:04** Selbst prüfen: Code, Bericht, Testprotokoll und Abrufliste liegen offen. Apache 2.0.

</details>

## Bilder

| | |
|---|---|
| ![Ein echter Gammablitz in mehreren Detektoren zugleich](medien/02-gammablitz.png) | ![Nachweiswahrscheinlichkeit nach Stärke, gewählte Einstellung gegen Grundeinstellung](medien/03-empfindlichkeit.png) |
| **Ein echter Gammablitz** (22.10.2025): Die Zählraten steigen in mehreren Detektoren zugleich; die Prüfgröße liegt weit über der Schwelle. | **Gemessen statt behauptet:** Künstliche Ausbrüche bekannter Stärke wiedergefunden; die optimierte Einstellung (blau) sieht schwächere als die Grundeinstellung (orange). |
| ![Gemessenes Muster über die Detektoren gegen das Modell](medien/05-antwortmodell.png) | ![Terminal: Tests, Auswertung, Kontrolllauf](medien/01-lauf.png) |
| **Physik geprüft:** Welcher Detektor wie viel sieht, folgt bei 54 echten Gammablitzen dem Modell der Detektorausrichtung (Korrelation im Median 0,92). | **Nachvollziehbar:** 24 Tests bestanden, ein Befehl für alle 61 Tage, ein zweiter Lauf mit identischen Ergebnissen. |

## Was Sie brauchen → läuft hier schon

Gemessen am 03.11.2026 auf den Daten vom 01.09. bis 31.10.2025; Quelle jeder Zahl:
[`nachweis/2026-11-03/ergebnisse.json`](nachweis/2026-11-03/ergebnisse.json).

| Was Sie brauchen | Läuft hier schon |
|---|---|
| Auffällige Ereignisse in einem Verbund vieler Detektoren erkennen | Koinzidenz über zwölf Detektoren: 49 von 54 Gammablitzen des offiziellen Katalogs gefunden (91 %; eine an den Bord-Trigger angelehnte Einstellung: 41) |
| Fehlalarme messen | Rate zufälliger Fehlalarme über Zeitverschiebungen gemessen; Grenzen je Tag als Nebenbedingung der Optimierung (an den Prüftagen überschritten, siehe „Ehrlich gesagt“) |
| Einstellungen begründet wählen statt raten | Optimierung über 162 Einstellungen, gewählt auf 31 Tagen, geprüft auf 30 anderen |
| Güte belegen statt behaupten | 21.805 künstliche Ausbrüche bekannter Stärke eingebaut und wiedergesucht: 61,2 % erkannt statt 55,6 % (Prüftage) |
| Zahlen, die jeder nachrechnen kann | 793 Rohdateien (9,3 GB) mit SHA-256 im Manifest; ein zweiter Lauf ergab dieselben Ergebnisse, Datei für Datei |
| Python unter Linux und Windows | 24 automatisierte Tests, bestanden in Docker unter Linux und in GitHub Actions unter Ubuntu und Windows (Python 3.11 bis 3.13) |

## In zwei Minuten selbst prüfen

1. **Bericht lesen:** [`nachweis/2026-11-03/BERICHT.pdf`](nachweis/2026-11-03/BERICHT.pdf), 8 Seiten, Testprotokoll in Abschnitt 8.
2. **Testlauf ansehen:** [GitHub Actions, Lauf 41718520863](https://github.com/panoherbe/gbmtransient/actions/runs/41718520863),
   sechs Jobs (Ubuntu und Windows, je Python 3.11, 3.12, 3.13), in jedem 24 von 24 Tests bestanden.
3. **Tests selbst ausführen** (Python ab 3.11, kein Netz nötig, rund 10 bis 30 Sekunden):

   ```bash
   git clone https://github.com/panoherbe/gbmtransient && cd gbmtransient
   pip install ".[plots,test]"
   python -m pytest tests
   ```

4. **Jede Zahl nachschlagen:** [`ergebnisse.json`](nachweis/2026-11-03/ergebnisse.json),
   [`katalog_abgleich.csv`](nachweis/2026-11-03/katalog_abgleich.csv) (jeder der 182 Katalogeinträge mit Ergebnis).

## Ihr Vorteil: Sie bewerten ein laufendes System, kein Versprechen

Das Verfahren ist gebaut, auf echten Messdaten ausgewertet und gegen eine unabhängige Referenz geprüft. Was es kann
und was nicht, steht in Zahlen im Bericht, einschließlich der Fehlschläge.

## Ehrlich gesagt

- **Zufällige Fehlalarme:** An den Prüftagen 0,33 je Tag statt der auf den Trainingstagen
  eingehaltenen 0,1. Eine Schwelle aus wenigen Wochen trägt die Schwankung von Tag zu Tag nicht; ein Betrieb
  bräuchte eine laufende Nachkalibrierung.
- **Verpasst:** 5 von 54 Gammablitzen, meist sehr kurze; für sie sind
  Daten mit 0,256 s Auflösung zu grob.
- **Nicht geprüft:** 80 Meldungen ohne Katalogeintrag in 61 Tagen, nur nach Merkmalen
  eingeordnet (Rand der Strahlungszone über dem Südatlantik, vermutlich schwache Sonneneruptionen, ohne Zuordnung).
- **Vereinfacht:** Das Modell der Detektoren (Kosinus, ohne Streuung); die Nachweisgrenzen gelten für dieses Modell.
- **Was ein Auftrag ergänzen würde:** Anpassung an die Daten und Meldewege des Auftraggebers, laufende
  Nachkalibrierung der Schwellen, Schutzabstand vor Datenlücken, Betrieb in Echtzeit.

## Kontakt

Panagiotis Herbe · [panoherbe.com](https://panoherbe.com)

---

## Für die IT: Technik

### Was er kann

- **Rohdaten holen und unverändert ablegen:** CTIME-Tagesdateien der zwölf NaI-Detektoren (acht Energiekanäle,
  0,256 s), die Positionsdatei je Tag und den Trigger- und Burstkatalog des HEASARC (TAP-Abfrage). Jede Datei mit
  URL, SHA-256 und Abrufzeit in `manifest.jsonl`; ausgewertet wird nur aus der Ablage, jede Datei wird vor dem
  Lesen gegen ihre Prüfsumme geprüft.
- **Datenmodell:** gemeinsames Zeitraster (0,256 s) mit Zählungen und Lebenszeit je Zelle, Detektor und
  Energieband; feinere Intervalle nach Bord-Triggern, Übergänge und Lücken gehen mit ihrer echten Messzeit ein.
- **Signifikanz je Detektor:** Likelihood-Quotienten-Test für Poisson-Zählungen (Li & Ma 1983) gegen einen
  zweiseitigen oder kausalen Untergrund aus Nebenfenstern, für Zeitskalen von 0,256 bis 8,192 s.
- **Koinzidenz über Detektoren:** k-größte Signifikanz der zwölf Detektoren (Ordnungsstatistik).
- **Fehlalarmrate durch Zufall:** Zeitverschiebungen der Detektoren gegeneinander (jedes Paar mindestens 43 min
  versetzt).
- **Gütemessung:** künstliche Ausbrüche (Pulsform nach Norris et al. 2005, Detektorantwort ∝ cos θ) werden in die
  echten Zählungen eingebaut und wiedergesucht.
- **Optimierung:** höchste Nachweiswahrscheinlichkeit unter zwei Nebenbedingungen (zufällige Fehlalarme, Zahl
  unbestätigter Kandidaten je Tag), vollständige Suche über 162 Konfigurationen, Auswahl auf Trainingstagen,
  Bewertung auf Testtagen, Bootstrap-Stabilität.
- **Abgleich:** Wiederfindung der Bord-Trigger je Einordnung im Katalog, Prüfung des Antwortmodells an echten
  Gammablitzen, Einordnung unbestätigter Kandidaten nach Merkmalen.
- **Bericht:** `ergebnisse.json` (jede Zahl des Berichts), `BERICHT.md`, Abbildungen, CSV-Listen.

### Start

```bash
pip install ".[plots]"
gbmtransient holen --daten daten --von 2025-09-01 --bis 2025-10-31     # 793 Dateien, rund 10 GB
gbmtransient auswerten --raw daten/raw/<Abrufdatum> --aus nachweis/<Datum> --prozesse 10
```

Mit Docker: `docker compose run --rm gbmtransient holen --daten /daten --von … --bis …`, dann
`docker compose run --rm gbmtransient auswerten --raw /daten/raw/<Abrufdatum> --aus /nachweis/<Datum>`;
Tests: `docker compose run --rm test` oder `python -m pytest tests`.

### Nachweis

`nachweis/<Datum>/` enthält je Auswertung Bericht (Markdown und PDF), Kennzahlen, Abbildungen, Listen und
Testprotokoll. Die Tests laufen bei jedem Push in GitHub Actions unter Ubuntu und Windows (Python 3.11 bis 3.13).

- `nachweis/2026-11-03/`: Auswertung von 61 Tagen (01.09.–31.10.2025): `BERICHT.pdf` (8 Seiten, Testprotokoll in
  Abschnitt 8), `ergebnisse.json` (jede Zahl des Berichts), Abbildungen, `ereignisse.csv`, `katalog_abgleich.csv`,
  `konfigurationen.csv`, `tage.csv`, das Manifest der Rohablage und die Testberichte der CI (`ci/`).

### Grenzen

Vereinfachte Detektorantwort (Kosinus, ohne Streuung), zwei Pulsformen, ein Zeitraum von 61 Tagen; unbestätigte
Kandidaten sind nach Merkmalen eingeordnet, nicht geprüft. Der Prototyp zeigt Verfahren und Messmethodik; er ist
kein Ersatz für die Auswertung des GBM-Teams.

### Daten und Lizenzen

- Fermi GBM, HEASARC (NASA/GSFC): https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/ und die Tabellen
  `fermigtrig`, `fermigbrst`. GBM-Daten waren zu keiner Zeit proprietär (Fermi Science Support Center, Fermi Data
  Policy); „HEASARC materials are all available freely for your use“ (HEASARC, Use of HEASARC Data).
- This research has made use of data and/or software provided by the High Energy Astrophysics Science Archive
  Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC.

Die Rohdaten sind nicht Teil des Repositorys. Code: Apache 2.0 (`LICENSE`, `NOTICE`).
