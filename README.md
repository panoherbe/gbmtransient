# gbmtransient (Prototyp)

Suche nach **kurzen Gamma- und Röntgentransienten** (Gammablitze, Sonneneruptionen, Magnetar-Ausbrüche) in den
offenen Messdaten des **Fermi Gamma-ray Burst Monitor (GBM)**, mit einem statistisch kontrollierten Trigger, dessen
Parameter durch ein Optimierungsverfahren gewählt und gegen den offiziellen GBM-Triggerkatalog geprüft werden.
Reines Python (NumPy, Astropy für FITS und VOTable, Matplotlib für Abbildungen), getestet unter Linux und Windows.
Lizenz Apache 2.0.

*Prototype: search for short gamma-ray transients in open Fermi GBM data with an optimised, statistically
controlled trigger, validated against the GBM trigger catalogue. Documentation in German. Not affiliated with NASA.*

## Was er kann

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

## Start

```bash
pip install ".[plots]"
gbmtransient holen --daten daten --von 2025-09-01 --bis 2025-10-31     # 793 Dateien, rund 10 GB
gbmtransient auswerten --raw daten/raw/<Abrufdatum> --aus nachweis/<Datum> --prozesse 10
```

Mit Docker: `docker compose run --rm gbmtransient holen --daten /daten --von … --bis …`, dann
`docker compose run --rm gbmtransient auswerten --raw /daten/raw/<Abrufdatum> --aus /nachweis/<Datum>`;
Tests: `docker compose run --rm test` oder `python -m pytest tests`.

## Nachweis

`nachweis/<Datum>/` enthält je Auswertung Bericht (Markdown und PDF), Kennzahlen, Abbildungen, Listen und
Testprotokoll. Die Tests laufen bei jedem Push in GitHub Actions unter Ubuntu und Windows (Python 3.11 bis 3.13).

- `nachweis/2026-11-03/`: Auswertung von 61 Tagen (01.09.–31.10.2025): `BERICHT.pdf` (8 Seiten, Testprotokoll in
  Abschnitt 8), `ergebnisse.json` (jede Zahl des Berichts), Abbildungen, `ereignisse.csv`, `katalog_abgleich.csv`,
  `konfigurationen.csv`, `tage.csv`, das Manifest der Rohablage und die Testberichte der CI (`ci/`).

## Grenzen

Vereinfachte Detektorantwort (Kosinus, ohne Streuung), zwei Pulsformen, ein Zeitraum von 61 Tagen; unbestätigte
Kandidaten sind nach Merkmalen eingeordnet, nicht geprüft. Der Prototyp zeigt Verfahren und Messmethodik; er ist
kein Ersatz für die Auswertung des GBM-Teams.

## Daten und Lizenzen

- Fermi GBM, HEASARC (NASA/GSFC): https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/ und die Tabellen
  `fermigtrig`, `fermigbrst`. GBM-Daten waren zu keiner Zeit proprietär (Fermi Science Support Center, Fermi Data
  Policy); „HEASARC materials are all available freely for your use“ (HEASARC, Use of HEASARC Data).
- This research has made use of data and/or software provided by the High Energy Astrophysics Science Archive
  Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC.

Die Rohdaten sind nicht Teil des Repositorys. Code: Apache 2.0 (`LICENSE`, `NOTICE`).
