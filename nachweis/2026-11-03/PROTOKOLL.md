# Testprotokoll gbmtransient 0.1.0, 03.11.2026

Verfasser: Panagiotis Herbe. Alle Angaben sind gemessen; die Belege liegen im Repository im Ordner
`nachweis/2026-11-03/` (`pytest-linux.txt`, `ci/`, `ergebnisse.json`, `manifest.jsonl.gz`) und im Lauf der GitHub
Actions.

## Umgebung

| | |
|---|---|
| Linux | eigener Server mit Ubuntu 24.04.5 LTS, 16 Kerne, 62 GB Arbeitsspeicher; Docker-Image auf `python:3.12-slim` mit Python 3.12.13, NumPy 2.5.4, Astropy 8.0.2, Matplotlib 3.11.2 (Versionen im `Dockerfile` festgelegt) |
| Windows und Linux in der CI | GitHub Actions, Workflow `.github/workflows/tests.yml`, Lauf 41718520863 zum Commit `1223818` (03.11.2026, 21:05–21:07 UTC): `windows-latest` (Image windows-2025-vs2026) mit Python 3.11.9, 3.12.10, 3.13.17 und `ubuntu-latest` (Image ubuntu-24.04) mit Python 3.11.17, 3.12.16, 3.13.17; Astropy 8.0.2, NumPy 2.5.4 (unter Python 3.11: 2.4.6). Windows ist ein Server der CI, kein eigener Windows-Rechner |
| Code | Rechnung mit Commit `da8a612`; Bericht, Listen und Abbildungen aus dem gespeicherten Zwischenstand neu geschrieben mit Commit `1223818` (geändert haben sich nur Berichtstext und Darstellung) |
| Daten | Rohablage `raw/2026-11-03/` auf dem Linux-Server; ihr Manifest (URL, SHA-256, Größe, Abrufzeit je Datei) liegt im Ordner als `manifest.jsonl.gz` |

## Befehle

```
(1) gbmtransient holen --daten daten --von 2025-09-01 --bis 2025-10-31 \
        --abrufdatum 2026-11-03
(2) python -m pytest -v tests
(3) gbmtransient auswerten --raw /daten/raw/2026-11-03 --aus /nachweis/2026-11-03 \
        --prozesse 10 --zwischen /zwischen/lauf-a.pkl.gz
(4) gbmtransient auswerten --raw /daten/raw/2026-11-03 --aus /kontrolle/2026-11-03-lauf2 \
        --prozesse 6 --zwischen /zwischen/lauf-b.pkl.gz
    python werkzeuge/vergleich.py <Ergebnis von (3)> kontrolle/2026-11-03-lauf2
(5) gbmtransient bericht --zwischen /zwischen/lauf-a.pkl.gz --aus /nachweis/2026-11-03
```

(1) bis (5) liefen in Docker auf dem Linux-Server (der Vergleich in (4) mit dem Python des Servers), (2) außerdem in
GitHub Actions.

## Ergebnisse

| Schritt | Ergebnis |
|---|---|
| (1) Abruf, 03.11.2026, 18:02–18:38 UTC | 856 Einträge im Manifest: 732 CTIME-Dateien (12 je Tag), 61 Positionsdateien, 61 Verzeichnisseiten, 2 Katalogabfragen; 9,26 GB; keine Fehlerdatei; 2.137 s mit 12 parallelen Abrufen |
| (2) Tests unter Linux (Docker) | 24 von 24 Tests bestanden in 10,4 s (`pytest-linux.txt`) |
| (2) Tests in GitHub Actions, sechs Jobs | in jedem Job 24 von 24 Tests bestanden; Dauer der Tests unter Windows 28,5 / 21,3 / 28,4 s, unter Ubuntu 16,1 / 17,4 / 10,2 s (Python 3.11 / 3.12 / 3.13); JUnit-Berichte in `ci/`; danach `gbmtransient --version` |
| (3) Auswertung, 19:27–20:53 UTC | Durchgang 1 (alle Konfigurationen, Zeitverschiebungen, künstliche Ausbrüche) 4.856 s, Durchgang 2 (Ereignisse, Messungen) 286 s, gleichzeitig mit (4). Ergebnis: `ergebnisse.json`, `BERICHT.md`, fünf Abbildungen mit 300 dpi, `ereignisse.csv` (390 Ereignisse), `katalog_abgleich.csv` (182 Trigger), `konfigurationen.csv` (162), `tage.csv` (61) |
| (4) Kontrolllauf mit 6 statt 10 Prozessen, 19:27–21:05 UTC | `ergebnisse.json` gleich (ohne Zeitstempel, Laufzeiten, Rechnerangaben); alle 10 Dateien (CSV, PNG, MD) Byte für Byte gleich |
| (5) Bericht aus dem Zwischenstand mit Commit `1223818` | `ergebnisse.json` gleich bis auf ein neues Feld (`bewertung.test.zwilling.bindend`); CSV-Dateien Byte für Byte gleich; geändert nur `BERICHT.md` und die Legende in Abb. 4 |

Der Bericht `BERICHT.pdf` ist aus `BERICHT.md` und diesem Protokoll gerendert (`werkzeuge/pdf.py`: Markdown, HTML,
PDF; Abbildungen eingebettet).

## Was die 24 Tests prüfen

| Bereich | Prüfung |
|---|---|
| Zeitrechnung | Missionszeit und UTC am GBM-Trigger zu GW170817 (MET 524666471) und am Tagesbeginn des 01.06.2025 (fünf Schaltsekunden) |
| Statistik | Formel nach Li & Ma gleich dem direkt berechneten Likelihood-Quotienten (5 Fälle); Nullverteilung auf Poisson-Zufallszahlen; Signifikanz mit Datenlücke; linearer Untergrund: zweiseitig unverzerrt, kausal verzerrt; Ordnungsstatistik; Ereignisbildung; Wilson- und Poisson-Intervalle gegen Tabellenwerte; logistische Anpassung findet vorgegebene Parameter |
| Datenmodell und Abruf | Zählungen und Messzeit bleiben auf dem Raster erhalten (feine Intervalle, Lücke, QUALITY); falsche Prüfsumme wird erkannt; höchste Dateiversion im Verzeichnis; Abruf ohne Netz mit Wiederaufnahme |
| Modell und Optimierung | Detektorachsen, Pulsform, Einbau mit dem erwarteten Mittelwert; Zeitverschiebungen trennen jedes Detektorpaar; 162 Konfigurationen; Schwelle aus Histogrammen |
| Ende zu Ende | zwei synthetische Tage im Format der GBM-Dateien und der Kataloge mit einem eingebauten Ausbruch: erkannt, Bericht, Listen und Abbildungen geschrieben; Bericht aus dem Zwischenstand mit denselben Kennzahlen |

## Was dieser Lauf nicht zeigt

- Die Auswertung der 61 Tage lief unter Linux. Unter Windows liefen die Tests einschließlich des Ende-zu-Ende-Laufs auf
  synthetischen Dateien, nicht die volle Auswertung.
- Die Laufzeiten gelten für zwei gleichzeitige Läufe auf einem Rechner.
