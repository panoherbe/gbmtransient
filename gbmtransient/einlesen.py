"""Rohablage lesen: CTIME-Dateien der zwölf NaI-Detektoren auf ein gemeinsames Zeitraster, Positionsdatei,
Trigger- und Burstkatalog.

Datenmodell eines Tages: Zeitraster mit 0,256 s Zellen ab ``t0`` (MET); je Detektor und Energieband die Zahl der
Ereignisse je Zelle und die Messzeit (Lebenszeit, EXPOSURE) je Zelle. CTIME liefert normal 0,256 s, nach einem
Bord-Trigger einige Minuten 0,064 s und an den Übergängen einzelne kürzere Intervalle; jedes Intervall wird der
Zelle seines Mittelpunkts zugeordnet, Zählungen und Messzeit werden addiert. Datenlücken (Südatlantische Anomalie,
Telemetrie) sind Zellen mit Messzeit 0. Intervalle mit QUALITY != 0 oder ohne Dauer werden verworfen.

Vor dem Lesen wird jede Datei gegen die SHA-256-Summe im Manifest geprüft.
"""

from __future__ import annotations

import hashlib
import io
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .zeit import met_aus_mjd

DT = 0.256  # s, Rasterweite = Standardauflösung von CTIME
NAI = ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb"]
# Energiebänder als CTIME-Kanäle (Kanalgrenzen stehen je Datei in EBOUNDS und werden mitgeführt)
BAENDER = {"12-27": (1, 1), "27-50": (2, 2), "50-300": (3, 4)}


@dataclass
class Tag:
    datum: str
    t0: float  # MET des Rasterbeginns
    zaehl: dict[str, np.ndarray]  # Band -> (12, n) int32
    dauer: np.ndarray  # (12, n) float32, Lebenszeit je Zelle in s
    ebounds: np.ndarray  # (12, 8, 2) keV
    dateien: list[dict] = field(default_factory=list)
    pos: dict[str, np.ndarray] | None = None

    @property
    def n(self) -> int:
        return self.dauer.shape[1]

    def zeiten(self) -> np.ndarray:
        """MET der Zellgrenzen 0..n (Index j = Ende eines Fensters, das bei Zelle j-1 endet)."""
        return self.t0 + DT * np.arange(self.n + 1)


class PruefsummenFehler(RuntimeError):
    pass


def manifest(raw: Path) -> dict[str, dict]:
    out = {}
    for zeile in (raw / "manifest.jsonl").read_text(encoding="utf-8").splitlines():
        if zeile.strip():
            e = json.loads(zeile)
            out[e["pfad"]] = e
    return out


def lies_geprueft(raw: Path, eintrag: dict) -> bytes:
    b = (raw / eintrag["pfad"]).read_bytes()
    if hashlib.sha256(b).hexdigest() != eintrag["sha256"]:
        raise PruefsummenFehler(f"SHA-256 stimmt nicht: {eintrag['pfad']}")
    return b


def tagesdateien(man: dict[str, dict]) -> dict[str, dict[str, dict]]:
    """Datum (JJJJ-MM-TT) -> {Detektor oder 'poshist' -> Manifesteintrag}."""
    tage: dict[str, dict[str, dict]] = {}
    for rel, e in man.items():
        m = re.search(r"daily/(\d{4})/(\d{2})/(\d{2})/current/glg_(?:ctime_(n[0-9ab])|poshist_all)_", rel)
        if m:
            tage.setdefault(f"{m.group(1)}-{m.group(2)}-{m.group(3)}", {})[m.group(4) or "poshist"] = e
    return dict(sorted(tage.items()))


def ctime_lesen(b: bytes) -> dict[str, np.ndarray]:
    from astropy.io import fits

    with fits.open(io.BytesIO(b), memmap=False) as h:
        d = h["SPECTRUM"].data
        eb = h["EBOUNDS"].data
        return {"zaehl": np.asarray(d["COUNTS"], dtype=np.int32), "exposure": np.asarray(d["EXPOSURE"], float),
                "quality": np.asarray(d["QUALITY"], int), "start": np.asarray(d["TIME"], float),
                "ende": np.asarray(d["ENDTIME"], float),
                "ebounds": np.stack([np.asarray(eb["E_MIN"], float), np.asarray(eb["E_MAX"], float)], axis=1)}


def auf_raster(start, ende, exposure, quality, zaehl, t0: float, n: int):
    """Intervalle der Zelle ihres Mittelpunkts zuordnen; liefert (Zählungen je Kanal (n, k), Lebenszeit (n,))."""
    ok = (quality == 0) & (ende > start) & (exposure > 0)
    mitte = 0.5 * (start[ok] + ende[ok])
    j = np.floor((mitte - t0) / DT).astype(np.int64)
    drin = (j >= 0) & (j < n)
    j = j[drin]
    auswahl = zaehl[ok][drin]
    z = np.zeros((n, zaehl.shape[1]), dtype=np.int64)
    for k in range(zaehl.shape[1]):
        z[:, k] = np.rint(np.bincount(j, weights=auswahl[:, k], minlength=n)[:n]).astype(np.int64)
    e = np.bincount(j, weights=exposure[ok][drin], minlength=n)[:n]
    return z, e


def positionen_lesen(b: bytes) -> dict[str, np.ndarray]:
    from astropy.io import fits

    with fits.open(io.BytesIO(b), memmap=False) as h:
        d = h[1].data
        return {"met": np.asarray(d["SCLK_UTC"], float), "flags": np.asarray(d["FLAGS"], int),
                "lat": np.asarray(d["SC_LAT"], float), "lon": np.asarray(d["SC_LON"], float),
                "q": np.stack([np.asarray(d[f"QSJ_{i}"], float) for i in range(1, 5)], axis=1)}


def tag_lesen(raw: Path, datum: str, eintraege: dict[str, dict], mit_position: bool = True) -> Tag:
    roh = {det: ctime_lesen(lies_geprueft(raw, eintraege[det])) for det in NAI if det in eintraege}
    if len(roh) < len(NAI):
        fehlt = [d for d in NAI if d not in roh]
        raise FileNotFoundError(f"{datum}: CTIME fehlt für {','.join(fehlt)}")
    t0 = min(float(r["start"][0]) for r in roh.values())
    t1 = max(float(r["ende"][-1]) for r in roh.values())
    n = int(np.ceil((t1 - t0) / DT))
    zaehl = {b: np.zeros((len(NAI), n), dtype=np.int32) for b in BAENDER}
    dauer = np.zeros((len(NAI), n), dtype=np.float32)
    eb = np.zeros((len(NAI), 8, 2))
    for i, det in enumerate(NAI):
        r = roh[det]
        z, e = auf_raster(r["start"], r["ende"], r["exposure"], r["quality"], r["zaehl"], t0, n)
        for b, (k0, k1) in BAENDER.items():
            zaehl[b][i] = z[:, k0:k1 + 1].sum(axis=1)
        dauer[i] = e
        eb[i] = r["ebounds"]
    pos = None
    if mit_position and "poshist" in eintraege:
        pos = positionen_lesen(lies_geprueft(raw, eintraege["poshist"]))
    dateien = [{"pfad": eintraege[k]["pfad"], "sha256": eintraege[k]["sha256"]} for k in sorted(eintraege)]
    return Tag(datum, t0, zaehl, dauer, eb, dateien, pos)


def katalog_lesen(raw: Path, man: dict[str, dict], tabelle: str) -> list[dict]:
    """VOTable aus der Ablage -> Liste von Zeilen (dict) mit zusätzlich 'met' (Triggerzeit in MET)."""
    from astropy.io.votable import parse_single_table
    import warnings

    eintrag = next(e for rel, e in man.items() if re.search(rf"tap/{tabelle}_.*\.vot$", rel))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = parse_single_table(io.BytesIO(lies_geprueft(raw, eintrag))).to_table()
    zeilen = []
    for r in t:
        z = {}
        for c in t.colnames:
            v = r[c]
            if np.ma.is_masked(v):
                v = None
            elif isinstance(v, (bytes, np.bytes_)):
                v = v.decode()
            elif isinstance(v, np.generic):
                v = v.item()
            z[c] = v
        z["met"] = float(met_aus_mjd(z["trigger_time"]))
        zeilen.append(z)
    return sorted(zeilen, key=lambda z: z["met"])
