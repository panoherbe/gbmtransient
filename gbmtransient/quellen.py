"""Rohdaten holen: Fermi-GBM-Tagesdaten (CTIME der zwölf NaI-Detektoren, Positionsdatei) und die Trigger- und
Burstkataloge aus dem HEASARC der NASA (Goddard Space Flight Center).

Jede Antwort wird zuerst unverändert auf die Platte geschrieben, mit URL, SHA-256, Größe und Abrufzeit im
Manifest. Ausgewertet wird danach nur aus diesen Dateien; ein zweiter Lauf liest dieselben Bytes und kommt zu
denselben Zahlen. Bereits vorhandene, im Manifest verzeichnete Dateien werden nicht noch einmal geholt.
"""

from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

from . import __version__
from .zeit import mjd_aus_datum

HEASARC = "https://heasarc.gsfc.nasa.gov"
TAGESDATEN = HEASARC + "/FTP/fermi/data/gbm/daily"
TAP = HEASARC + "/xamin/vo/tap/sync"
NAI = ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb"]
KATALOGE = ("fermigtrig", "fermigbrst")  # Triggerkatalog (alle Trigger mit Einordnung), Burstkatalog (GRB)
USER_AGENT = f"gbmtransient/{__version__} (+https://github.com/panoherbe/gbmtransient)"


def jetzt() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def tag_url(d: date) -> str:
    return f"{TAGESDATEN}/{d:%Y/%m/%d}/current/"


def tap_url(abfrage: str) -> str:
    return TAP + "?" + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable",
                                               "QUERY": abfrage})


def katalog_abfrage(tabelle: str, von: date, bis: date) -> str:
    # ein Tag Rand auf jeder Seite; trigger_time ist MJD (UTC)
    return (f"SELECT * FROM {tabelle} WHERE trigger_time >= {mjd_aus_datum(von) - 1} "
            f"AND trigger_time < {mjd_aus_datum(bis) + 2} ORDER BY trigger_time")


def dateien_im_verzeichnis(html: str, d: date) -> dict[str, str]:
    """Dateinamen der CTIME-Dateien (je NaI-Detektor) und der Positionsdatei, jeweils die höchste Version."""
    tag = d.strftime("%y%m%d")
    beste: dict[str, tuple[int, str]] = {}
    for name in set(re.findall(r"glg_(?:ctime_n[0-9ab]|poshist_all)_\d{6}_v\d{2}\.(?:pha|fit)", html)):
        m = re.match(r"glg_(ctime_(n[0-9ab])|poshist_all)_(\d{6})_v(\d{2})\.", name)
        if not m or m.group(3) != tag:
            continue
        schluessel = m.group(2) or "poshist"
        v = int(m.group(4))
        if schluessel not in beste or v > beste[schluessel][0]:
            beste[schluessel] = (v, name)
    return {k: v[1] for k, v in sorted(beste.items())}


def hole(url: str, timeout: float = 300, versuche: int = 4) -> bytes:
    letzter: Exception | None = None
    for i in range(versuche):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            letzter = e
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            letzter = e
        time.sleep(5 * (i + 1))
    raise RuntimeError(f"{url}: {letzter}")


class Archiv:
    """Rohablage unter <wurzel>/raw/<abrufdatum>/ mit manifest.jsonl; der Pfad einer Datei folgt ihrer URL."""

    def __init__(self, wurzel: Path, abrufdatum: str):
        self.dir = wurzel / "raw" / abrufdatum
        self.dir.mkdir(parents=True, exist_ok=True)
        self.manifest = self.dir / "manifest.jsonl"
        self._lock = threading.Lock()
        self.bekannt: dict[str, dict] = {}
        if self.manifest.exists():
            for zeile in self.manifest.read_text(encoding="utf-8").splitlines():
                if zeile.strip():
                    e = json.loads(zeile)
                    self.bekannt[e["pfad"]] = e

    @staticmethod
    def rel_aus_url(url: str) -> str:
        teile = urllib.parse.urlsplit(url)
        pfad = teile.path.lstrip("/")
        if pfad.endswith("/"):
            pfad += "index.html"
        return "heasarc/" + pfad

    def vorhanden(self, rel: str) -> bool:
        e = self.bekannt.get(rel)
        p = self.dir / rel
        return e is not None and p.exists() and p.stat().st_size == e["bytes"]

    def ablegen(self, rel: str, url: str, daten: bytes, dauer: float, abfrage: str | None = None) -> dict:
        ziel = self.dir / rel
        ziel.parent.mkdir(parents=True, exist_ok=True)
        tmp = ziel.with_name(ziel.name + ".part")
        tmp.write_bytes(daten)
        tmp.replace(ziel)
        e = {"pfad": rel, "url": url, "bytes": len(daten), "sha256": hashlib.sha256(daten).hexdigest(),
             "abgerufen": jetzt(), "dauer_s": round(dauer, 2)}
        if abfrage:
            e["abfrage"] = abfrage
        with self._lock:
            with self.manifest.open("a", encoding="utf-8") as f:
                f.write(json.dumps(e, ensure_ascii=False) + "\n")
            self.bekannt[rel] = e
        return e

    def holen(self, url: str, abfrage: str | None = None, neu: bool = False) -> dict:
        rel = self.rel_aus_url(url) if abfrage is None else abfrage
        if not neu and self.vorhanden(rel):
            return self.bekannt[rel]
        t = time.time()
        daten = hole(url)
        return self.ablegen(rel, url, daten, time.time() - t, abfrage=None if abfrage is None else url)


def tage(von: date, bis: date) -> list[date]:
    return [von + timedelta(days=i) for i in range((bis - von).days + 1)]


def holen(daten: Path, von: date, bis: date, abrufdatum: str | None = None, arbeiter: int = 12) -> dict:
    abrufdatum = abrufdatum or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    a = Archiv(daten, abrufdatum)
    t0 = time.time()
    fehler: list[str] = []
    # 1. Kataloge (klein): Triggerkatalog und Burstkatalog für den Zeitraum
    for tab in KATALOGE:
        q = katalog_abfrage(tab, von, bis)
        rel = f"heasarc/tap/{tab}_{von:%Y-%m-%d}_{bis:%Y-%m-%d}.vot"
        if not a.vorhanden(rel):
            t = time.time()
            a.ablegen(rel, tap_url(q), hole(tap_url(q)), time.time() - t, abfrage=q)
    # 2. Verzeichnisse der Tage, dann die Dateien
    auftraege: list[str] = []
    for d in tage(von, bis):
        try:
            e = a.holen(tag_url(d), neu=False)
            html = (a.dir / e["pfad"]).read_text(encoding="utf-8", errors="replace")
        except Exception as ex:  # noqa: BLE001 - jeder Fehler wird protokolliert, der Lauf geht weiter
            fehler.append(f"{d}: Verzeichnis: {ex}")
            continue
        namen = dateien_im_verzeichnis(html, d)
        fehlt = [k for k in NAI + ["poshist"] if k not in namen]
        if fehlt:
            fehler.append(f"{d}: fehlt im Verzeichnis: {','.join(fehlt)}")
        auftraege += [tag_url(d) + n for n in namen.values()]

    def eins(url: str) -> str | None:
        try:
            a.holen(url)
            return None
        except Exception as ex:  # noqa: BLE001
            return f"{url}: {ex}"

    with cf.ThreadPoolExecutor(arbeiter) as ex:
        for i, r in enumerate(ex.map(eins, auftraege), 1):
            if r:
                fehler.append(r)
            if i % 120 == 0:
                print(f"  {i}/{len(auftraege)} Dateien, {time.time() - t0:.0f} s", flush=True)
    if fehler:
        (a.dir / "fehler.txt").write_text("\n".join(fehler) + "\n", encoding="utf-8")
    bytes_ = sum(e["bytes"] for e in a.bekannt.values())
    print(f"{len(a.bekannt)} Dateien im Manifest, {bytes_ / 1e9:.2f} GB, {len(fehler)} Fehler, "
          f"{time.time() - t0:.0f} s; Ablage {a.dir}")
    return {"dateien": len(a.bekannt), "bytes": bytes_, "fehler": fehler}
