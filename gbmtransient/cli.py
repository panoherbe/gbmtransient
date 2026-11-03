"""Kommandozeile: gbmtransient holen | auswerten."""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

from . import __version__


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="gbmtransient",
                                 description="Suche nach kurzen Gamma- und Röntgentransienten in Fermi-GBM-Daten")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="befehl", required=True)

    h = sub.add_parser("holen", help="Tagesdaten und Kataloge vom HEASARC unverändert ablegen")
    h.add_argument("--daten", type=Path, required=True, help="Datenordner (darin raw/<Abrufdatum>/)")
    h.add_argument("--von", type=date.fromisoformat, required=True, help="erster Tag (JJJJ-MM-TT)")
    h.add_argument("--bis", type=date.fromisoformat, required=True, help="letzter Tag (JJJJ-MM-TT)")
    h.add_argument("--abrufdatum", help="Name des Rohordners (Standard: heute, UTC)")
    h.add_argument("--arbeiter", type=int, default=12, help="gleichzeitige Abrufe (Standard 12)")

    a = sub.add_parser("auswerten", help="Rohablage auswerten: Suche, Optimierung, Abgleich, Bericht")
    a.add_argument("--raw", type=Path, required=True, help="Rohordner raw/<Abrufdatum>")
    a.add_argument("--aus", type=Path, required=True, help="Zielordner")
    a.add_argument("--prozesse", type=int, default=8, help="parallele Prozesse (Standard 8)")
    a.add_argument("--max-tage", type=int, help="nur die ersten N Tage (zum Ausprobieren)")
    a.add_argument("--ohne-abbildungen", action="store_true")

    args = ap.parse_args(argv)
    if args.befehl == "holen":
        from .quellen import holen
        r = holen(args.daten, args.von, args.bis, args.abrufdatum, args.arbeiter)
        return 1 if r["fehler"] else 0
    from .bericht import auswerten
    r = auswerten(args.raw, args.aus, args.prozesse, args.max_tage, not args.ohne_abbildungen)
    print(f"{r['daten']['tage']} Tage ausgewertet, Bericht: {args.aus / 'BERICHT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
