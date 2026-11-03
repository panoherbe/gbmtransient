"""Zwei Auswertungen derselben Rohablage vergleichen.

ergebnisse.json ohne Erstellungszeit, Laufzeiten und Rechnerangaben; CSV-Dateien und Abbildungen Byte für Byte.

    python werkzeuge/vergleich.py nachweis/<lauf-1> nachweis/<lauf-2>
"""

import hashlib
import json
import sys
from pathlib import Path

IGNORIEREN = {"erstellt", "rechenzeit_s", "rechner"}


def laden(p: Path) -> dict:
    r = json.loads((p / "ergebnisse.json").read_text(encoding="utf-8"))
    return {k: v for k, v in r.items() if k not in IGNORIEREN}


def unterschiede(a, b, pfad=""):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            yield from unterschiede(a.get(k), b.get(k), f"{pfad}.{k}")
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            yield from unterschiede(x, y, f"{pfad}[{i}]")
    elif a != b:
        yield pfad


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    a, b = map(Path, sys.argv[1:3])
    diff = list(unterschiede(laden(a), laden(b)))
    print(f"ergebnisse.json: {'gleich' if not diff else f'{len(diff)} Unterschiede'} (ohne Zeitstempel, Laufzeiten)")
    for d in diff[:20]:
        print("  ", d)
    dateien = sorted(p.relative_to(a) for p in a.rglob("*") if p.suffix in (".csv", ".png", ".md"))
    abweichend = [str(f) for f in dateien if not (b / f).exists() or sha(a / f) != sha(b / f)]
    print(f"{len(dateien)} Dateien (CSV, PNG, MD) verglichen, {len(abweichend)} abweichend: {abweichend}")
    return 1 if diff else 0


if __name__ == "__main__":
    sys.exit(main())
