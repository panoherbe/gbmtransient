"""Zweiter Durchgang je Tag, nur für die gewählten Konfigurationen: Ereignisliste mit Merkmalen, Messungen an den
Katalog-Triggern (Spitzenrate, Muster über die Detektoren, Härte) und die Lichtkurve eines Beispiels."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import einlesen, simulation
from .einlesen import DT
from .statistik import ereignisse, kumulieren, ordnungsstatistik, ueberschuss
from .suche import (ABGLEICH, BAND_SAETZE, G, RAND, TRIGGERFENSTER, VERBINDEN, ZEITSKALEN, ZS_SAETZE, Konfig,
                    abgedeckt, alle_signifikanzen)

MESS_W = 128  # Untergrund für Messungen: zweiseitig, 32,8 s
MESSFENSTER = (-5.0, 60.0)  # s um die Triggerzeit


def glieder(cfg: Konfig) -> list[tuple]:
    return [(b, m, cfg.modus, cfg.w) for b in BAND_SAETZE[cfg.band] for m in ZEITSKALEN if m <= ZS_SAETZE[cfg.zs]]


def _luecken_abstand(dauer: np.ndarray, j: int, min_det: int = 6) -> float:
    """Abstand in s von Zelle j zur nächsten Zelle mit weniger als min_det Detektoren (Datenlücke, SAA)."""
    luecke = np.flatnonzero((dauer > 0).sum(axis=0) < min_det)
    if luecke.size == 0:
        return float("inf")
    return float(np.min(np.abs(luecke - j)) * DT)


def _position(pos: dict | None, met: float) -> dict:
    if pos is None or pos["met"].size == 0:
        return {}
    i = int(np.clip(np.searchsorted(pos["met"], met), 0, pos["met"].size - 1))
    f = int(pos["flags"][i])
    return {"sonne_sichtbar": bool(f & 1), "saa": bool(f & 2), "breite": round(float(pos["lat"][i]), 2),
            "laenge": round(float(pos["lon"][i]), 2)}


def tag_details(aufgabe: tuple) -> dict:
    raw, datum, eintraege, tag_trigger, alle_trigger, gewaehlt, beispiel = aufgabe
    tag = einlesen.tag_lesen(Path(raw), datum, eintraege, mit_position=True)
    n = tag.n
    zeit = tag.t0 + DT * np.arange(n + 1)
    voll = float(np.median(tag.dauer[tag.dauer > 0]))
    kombis = set()
    for _, cfg, _ in gewaehlt:
        kombis |= set(glieder(cfg))
    S = alle_signifikanzen(tag.zaehl, tag.dauer, kombis)
    ce = kumulieren(tag.dauer, RAND)
    mess = {}
    for b in ("27-50", "50-300"):
        cz = kumulieren(tag.zaehl[b], RAND)
        for m in (1, 4, 16):
            mess[(b, m)] = ueberschuss(cz, ce, n, RAND, m, MESS_W, G, voll)
        del cz
    alle_t = np.array(sorted(t["met"] for t in alle_trigger), float)
    alle_name = [t["trigger_name"] for t in sorted(alle_trigger, key=lambda t: t["met"])]
    alle_typ = [t["trigger_type"] for t in sorted(alle_trigger, key=lambda t: t["met"])]

    def haerte(j: int, dets) -> float | None:
        h = np.nansum(mess[("50-300", 16)][0][list(dets), j])
        w = np.nansum(mess[("27-50", 16)][0][list(dets), j])
        return float(w / h) if h > 0 else None

    out: dict = {"datum": datum, "ereignisse": [], "trigger": [], "beispiel": None}
    for label, cfg, s in gewaehlt:
        gl = glieder(cfg)
        stapel = np.stack([ordnungsstatistik(S[g], (cfg.k,))[cfg.k] for g in gl])
        c = stapel.max(axis=0)
        welches = stapel.argmax(axis=0)
        a, e, p, spitze = ereignisse(c, s, VERBINDEN)
        for ai, ei, pi, sp in zip(a, e, p, spitze):
            g = gl[int(welches[pi])]
            sig = S[g][:, pi]
            ueber = np.flatnonzero(np.nan_to_num(sig, nan=-np.inf) >= s)
            hell = np.argsort(np.nan_to_num(mess[("50-300", 16)][0][:, pi], nan=-np.inf))[::-1][:3]
            lo = np.searchsorted(alle_t, zeit[ai] - ABGLEICH[0])
            hi = np.searchsorted(alle_t, zeit[ei] + ABGLEICH[1], side="right")
            ev = {"konfig": label, "met_anfang": float(zeit[ai]), "met_ende": float(zeit[ei]),
                  "met_spitze": float(zeit[pi]), "spitze": round(float(sp), 2), "band": g[0],
                  "zeitskala_s": round(g[1] * DT, 3), "detektoren_ueber": [einlesen.NAI[d] for d in ueber],
                  "dauer_s": round(float(zeit[ei] - zeit[ai]) + g[1] * DT, 2),
                  "katalog": [{"name": alle_name[i], "typ": alle_typ[i]} for i in range(lo, hi)],
                  "luecke_abstand_s": _luecken_abstand(tag.dauer, int(pi)), "haerte": haerte(int(pi), hell)}
            ev.update(_position(tag.pos, float(zeit[pi])))
            out["ereignisse"].append(ev)

    for t in tag_trigger:
        ab = abgedeckt(tag.dauer, tag.t0, t["met"])
        rec = {"name": t["trigger_name"], "typ": t["trigger_type"], "met": t["met"], "abgedeckt": ab}
        lo = np.searchsorted(zeit, t["met"] + MESSFENSTER[0])
        hi = np.searchsorted(zeit, t["met"] + MESSFENSTER[1], side="right")
        if ab and hi > lo:
            x4, e4 = mess[("50-300", 4)]
            rate = x4[:, lo:hi] / e4[:, lo:hi]
            spitzen = np.nanmax(np.where(np.isfinite(rate), rate, -np.inf), axis=1)
            x16 = mess[("50-300", 16)][0][:, lo:hi]
            summe = np.nansum(np.where(np.isfinite(x16), x16, 0.0), axis=0)
            jm = lo + int(np.argmax(summe))
            muster = np.nan_to_num(mess[("50-300", 16)][0][:, jm], nan=0.0)
            hell = np.argsort(muster)[::-1][:3]
            rec.update({"spitzenrate_cts_s": float(np.max(spitzen)), "hellster": einlesen.NAI[int(np.argmax(spitzen))],
                        "muster": [round(float(v), 1) for v in muster], "haerte": haerte(jm, hell),
                        "met_muster": float(zeit[jm])})
            th, ph = t.get("theta"), t.get("phi")
            if th is not None and ph is not None:
                modell = simulation.antwort(float(th), float(ph))
                rec["modell"] = [round(float(v), 3) for v in modell]
                if np.std(muster) > 0 and np.std(modell) > 0:
                    rec["korrelation"] = float(np.corrcoef(muster, modell)[0, 1])
                rec["hellster_in_modell_top2"] = bool(int(np.argmax(muster)) in np.argsort(modell)[::-1][:2])
        out["trigger"].append(rec)
        if beispiel and t["trigger_name"] == beispiel:
            lo = np.searchsorted(zeit, t["met"] - 40.0)
            hi = np.searchsorted(zeit, t["met"] + 80.0)
            x1 = mess[("50-300", 1)][0]
            dets = np.argsort(np.nan_to_num(rec.get("muster", np.zeros(12)), nan=0.0))[::-1][:3]
            zell = tag.zaehl["50-300"][:, lo:hi]
            label, cfg, s = gewaehlt[0]
            stapel = np.stack([ordnungsstatistik(S[g], (cfg.k,))[cfg.k] for g in glieder(cfg)])
            out["beispiel"] = {"name": beispiel, "met_trigger": t["met"], "zeit_rel": (zeit[lo:hi] - t["met"]).tolist(),
                               "detektoren": [einlesen.NAI[d] for d in dets],
                               "zaehl": zell[dets].tolist(), "untergrund": (zell[dets] - x1[dets, lo + 1:hi + 1]).tolist(),
                               "koinzidenz": np.where(np.isfinite(stapel.max(0)[lo:hi]), stapel.max(0)[lo:hi], np.nan).tolist(),
                               "schwelle": s, "konfig": label}
    return out
