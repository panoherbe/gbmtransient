"""Suche je Tag. Für alle 162 Konfigurationen eines Triggers (k, Untergrund, Zeitskalen, Band) in einem Durchgang:

1. Nulllauf auf den echten Daten: Ereignisse über der Bodenschwelle, getrennt nach "im Triggerkatalog" und
   "unbestätigt", und die höchste Koinzidenz im Fenster jedes Katalog-Triggers;
2. Zeitverschiebungen: jeder Detektor wird um ein anderes Vielfaches von 258 s verschoben; echte gleichzeitige
   Signale fallen auseinander, übrig bleiben zufällige Koinzidenzen (Fehlalarmrate durch Zufall);
3. künstliche Ausbrüche: Pulse mit bekannter Stärke und Richtung werden in die echten Zählungen eingebaut und
   wiedergesucht (Nachweiswahrscheinlichkeit).

Jeder Tag wird für sich gerechnet und liefert nur Histogramme und Maxima je Konfiguration; die Auswahl der
Schwelle und der Konfiguration (optimierung.py) liest nur diese Zusammenfassungen.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import einlesen, simulation
from .einlesen import DT
from .statistik import ereignisse, kumulieren, ordnungsstatistik, signifikanz

TRIGGERBAENDER = ("27-50", "50-300")
ZEITSKALEN = (1, 2, 4, 8, 16, 32)  # Zellen à 0,256 s: 0,256 ... 8,192 s
HINTERGRUENDE = (("zweiseitig", 32), ("zweiseitig", 64), ("zweiseitig", 128),
                 ("kausal", 32), ("kausal", 64), ("kausal", 128))  # Fensterbreite w in Zellen (8,2 / 16,4 / 32,8 s)
G = 8  # Abstand zwischen Signal- und Untergrundfenster, Zellen (2,048 s)
KS = (1, 2, 3)
ZS_SAETZE = {"0,256-1,024 s": 4, "0,256-4,096 s": 16, "0,256-8,192 s": 32}  # alle Zeitskalen bis einschließlich
BAND_SAETZE = {"50-300": ("50-300",), "27-50": ("27-50",), "27-300": ("27-50", "50-300")}
RAND = max(ZEITSKALEN) + G + max(w for _, w in HINTERGRUENDE) + 1
BODEN = 3.5
KANTEN = np.round(np.arange(BODEN, 25.0 + 1e-9, 0.02), 2)  # Histogrammkanten der Ereignis-Spitzen
VERBINDEN = 40  # Zellen (10,24 s): näher beieinander liegende Überschreitungen sind ein Ereignis
N_VERSCHIEBUNGEN = 20
MIN_ABSTAND = 0.03  # kleinster gegenseitiger Versatz zweier Detektoren, Anteil des Tages (0,03 = 43 min)


def verschiebungsanteile(anzahl: int = N_VERSCHIEBUNGEN, eps: float = MIN_ABSTAND) -> np.ndarray:
    """Anteile phi_v des Tages: Detektor d wird um d*phi_v verschoben. Gewählt so, dass q*phi_v für q = 1..11
    (Abstand zweier Detektornummern) mindestens eps von einer ganzen Zahl entfernt ist: Jedes Detektorpaar ist in
    jedem Durchlauf um mindestens eps Tage gegeneinander versetzt, auch lange Ereignisse (Sonneneruptionen) fallen
    auseinander."""
    kand = np.arange(1, 20000) / 20000
    q = np.arange(1, 12)[:, None]
    abst = np.abs(q * kand - np.round(q * kand)).min(axis=0)
    gut = kand[abst >= eps]
    return gut[np.linspace(0, gut.size - 1, anzahl).round().astype(int)]


PHI = verschiebungsanteile()
ABGLEICH = (20.0, 5.0)  # Ereignis [a, e] passt zu einem Katalog-Trigger t, wenn a - 20 s <= t <= e + 5 s
TRIGGERFENSTER = (-5.0, 20.0)  # s um die Triggerzeit, in denen die höchste Koinzidenz gesucht wird
INJ_ABSTAND = 200.0  # s zwischen zwei künstlichen Ausbrüchen
INJ_FENSTER = (-1.0, 30.0)  # s ab Pulsbeginn
INJ_NULLANTEIL = 0.05  # Anteil der Plätze ohne Puls (Kontrolle: zufällige Treffer im Fenster)
SEED = 20251101
REFERENZ = ("50-300", 4, "zweiseitig", 128)  # Kombination für die Prüfung der Nullverteilung
S_KANTEN = np.round(np.arange(-10.0, 30.0 + 1e-9, 0.1), 1)


@dataclass(frozen=True)
class Konfig:
    k: int
    modus: str
    w: int
    zs: str
    band: str

    def name(self) -> str:
        return f"k={self.k}, {self.modus} {self.w * DT:.1f} s, {self.zs}, {self.band} keV"

    def kurz(self) -> dict:
        return {"k": self.k, "untergrund": self.modus, "w_s": round(self.w * DT, 3), "zeitskalen": self.zs,
                "band_kev": self.band}


KONFIGS = [Konfig(k, modus, w, zs, band) for k in KS for modus, w in HINTERGRUENDE for zs in ZS_SAETZE
           for band in BAND_SAETZE]
INDEX = {c: i for i, c in enumerate(KONFIGS)}


def alle_signifikanzen(zaehl: dict[str, np.ndarray], dauer: np.ndarray, kombis=None) -> dict:
    n = dauer.shape[1]
    voll = float(np.median(dauer[dauer > 0]))
    ce = kumulieren(dauer, RAND)
    S = {}
    for b in TRIGGERBAENDER:
        cz = kumulieren(zaehl[b], RAND)
        for modus, w in HINTERGRUENDE:
            for m in ZEITSKALEN:
                if kombis is None or (b, m, modus, w) in kombis:
                    S[(b, m, modus, w)] = signifikanz(cz, ce, n, RAND, m, w, G, modus, voll)
        del cz
    return S


def konfig_reihen(ck: dict):
    """(Index, Konfiguration, Reihe C(j)) für alle Konfigurationen; C = max über Bänder und Zeitskalen des Satzes."""
    for k in KS:
        for modus, w in HINTERGRUENDE:
            m_band: dict[str, dict[str, np.ndarray]] = {}
            for b in TRIGGERBAENDER:
                laufend, schnappschuss = None, {}
                for m in ZEITSKALEN:
                    r = ck[(b, m, modus, w)][k]
                    laufend = r.copy() if laufend is None else np.maximum(laufend, r)
                    for name, bis in ZS_SAETZE.items():
                        if bis == m:
                            schnappschuss[name] = laufend.copy()
                m_band[b] = schnappschuss
            for zs in ZS_SAETZE:
                for band, glieder in BAND_SAETZE.items():
                    c = m_band[glieder[0]][zs]
                    for g in glieder[1:]:
                        c = np.maximum(c, m_band[g][zs])
                    cfg = Konfig(k, modus, w, zs, band)
                    yield INDEX[cfg], cfg, c


def _hist(werte: np.ndarray) -> np.ndarray:
    h, _ = np.histogram(np.clip(werte, KANTEN[0], KANTEN[-1] - 1e-6), bins=KANTEN)
    return h.astype(np.int32)


def _verschieben(s: np.ndarray, v: int) -> np.ndarray:
    """Durchlauf v (1..N): Detektor d um round(d * phi_v * n) Zellen zyklisch verschoben."""
    n = s.shape[1]
    return np.stack([np.roll(s[d], int(round(d * PHI[v - 1] * n)) % n) for d in range(s.shape[0])])


def _fenster_max(c: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """max(c[lo_i:hi_i]) für sortierte, disjunkte Fenster; leeres Fenster -> -inf."""
    out = np.full(lo.size, -np.inf, dtype=np.float32)
    voll = hi > lo
    if not np.any(voll):
        return out
    l, h = lo[voll], hi[voll]
    idx = np.empty(2 * l.size, dtype=np.int64)
    idx[0::2], idx[1::2] = l, h
    if idx[-1] >= c.size:
        idx = idx[:-1]
    out[voll] = np.maximum.reduceat(c, idx)[0::2][:l.size]
    return out


def abgedeckt(dauer: np.ndarray, t0: float, t: float, fenster=TRIGGERFENSTER, min_det: int = 6) -> bool:
    """Liegen im Fenster um t mindestens 80 % der Zellen mit Daten von mindestens min_det Detektoren?"""
    a = int(np.floor((t + fenster[0] - t0) / DT))
    e = int(np.ceil((t + fenster[1] - t0) / DT))
    if a < 0 or e > dauer.shape[1]:
        return False
    live = (dauer[:, a:e] > 0).sum(axis=0) >= min_det
    return bool(live.mean() >= 0.8)


def injektionen_planen(tag: einlesen.Tag, alle_trigger_met: np.ndarray, rng: np.random.Generator) -> dict:
    """Plätze alle 200 s (±20 s), ohne Katalog-Trigger im Abstand von 200 s und ohne Datenlücke im Umfeld."""
    t_end = tag.t0 + tag.n * DT
    plaetze = tag.t0 + 150.0 + INJ_ABSTAND * np.arange(int((t_end - tag.t0 - 300) // INJ_ABSTAND))
    plaetze = plaetze + rng.uniform(-20, 20, plaetze.size)
    alle_live = (tag.dauer > 0).all(axis=0)
    cl = np.concatenate([[0], np.cumsum(~alle_live)])
    a = np.floor((plaetze - 60 - tag.t0) / DT).astype(int)
    e = np.ceil((plaetze + 100 - tag.t0) / DT).astype(int)
    ok = (a >= 0) & (e <= tag.n)
    ok[ok] &= (cl[e[ok]] - cl[a[ok]]) == 0
    if alle_trigger_met.size:
        i = np.searchsorted(alle_trigger_met, plaetze)
        nah = np.zeros(plaetze.size, bool)
        for d in (i - 1, i):
            d = np.clip(d, 0, alle_trigger_met.size - 1)
            nah |= np.abs(alle_trigger_met[d] - plaetze) < 200
        ok &= ~nah
    t = plaetze[ok]
    n = t.size
    klasse = (rng.random(n) < 0.5).astype(np.int8)  # 0 = kurz, 1 = lang
    namen = ["kurz", "lang"]
    amp = np.empty(n)
    for i, kl in enumerate(klasse):
        k = simulation.KLASSEN[namen[kl]]
        amp[i] = np.exp(rng.uniform(np.log(k.a_min), np.log(k.a_max)))
    amp[rng.random(n) < INJ_NULLANTEIL] = 0.0
    theta = np.degrees(np.arccos(rng.uniform(-1, 1, n)))
    phi = rng.uniform(0, 360, n)
    faktor = np.stack([simulation.antwort(th, ph) for th, ph in zip(theta, phi)]) if n else np.zeros((0, 12))
    return {"met": t, "klasse": klasse, "amplitude": amp, "theta": theta, "phi": phi, "faktor": faktor}


def tag_rechnen(aufgabe: tuple) -> dict:
    raw, datum, eintraege, tag_trigger, alle_trigger_met, tag_index = aufgabe
    t_start = time.time()
    tag = einlesen.tag_lesen(Path(raw), datum, eintraege, mit_position=False)
    n = tag.n
    zeit = tag.t0 + DT * np.arange(n + 1)  # Ende des Signalfensters je Index j
    nk, nb = len(KONFIGS), KANTEN.size - 1
    erg: dict = {"datum": datum, "tag_index": tag_index, "t0": tag.t0, "n_zellen": n,
                 "lebenszeit_s": float(((tag.dauer > 0).sum(axis=0) >= 6).sum() * DT)}
    trig_met = np.array([t["met"] for t in tag_trigger], float)
    alle = np.sort(np.asarray(alle_trigger_met, float))

    # 1. Nulllauf
    S = alle_signifikanzen(tag.zaehl, tag.dauer)
    ref = S[REFERENZ]
    erg["s_hist"] = np.stack([np.histogram(r[np.isfinite(r)], bins=S_KANTEN)[0] for r in ref]).astype(np.int64)
    ck = {kombi: ordnungsstatistik(s, KS) for kombi, s in S.items()}
    h_unb = np.zeros((nk, nb), np.int32)
    h_bes = np.zeros((nk, nb), np.int32)
    gueltig = np.zeros(nk)
    tmax = np.full((nk, trig_met.size), -np.inf, np.float32)
    t_lo = np.searchsorted(zeit, trig_met + TRIGGERFENSTER[0])
    t_hi = np.searchsorted(zeit, trig_met + TRIGGERFENSTER[1], side="right")
    for i, cfg, c in konfig_reihen(ck):
        gueltig[i] = np.isfinite(c).sum() * DT
        a, e, _, spitze = ereignisse(c, BODEN, VERBINDEN)
        if a.size:
            lo = np.searchsorted(alle, zeit[a] - ABGLEICH[0])
            hi = np.searchsorted(alle, zeit[e] + ABGLEICH[1], side="right")
            passt = hi > lo
            h_unb[i] = _hist(spitze[~passt])
            h_bes[i] = _hist(spitze[passt])
        for t in range(trig_met.size):
            if t_hi[t] > t_lo[t]:
                tmax[i, t] = np.max(c[t_lo[t]:t_hi[t]])
    erg.update(hist_unbestaetigt=h_unb, hist_bestaetigt=h_bes, gueltig_s=gueltig, trigger_max=tmax,
               trigger=[{"name": t["trigger_name"], "met": t["met"], "typ": t["trigger_type"],
                         "abgedeckt": abgedeckt(tag.dauer, tag.t0, t["met"])} for t in tag_trigger])
    del ck

    # 2. Zeitverschiebungen
    h_ver = np.zeros((nk, nb), np.int32)
    g_ver = np.zeros(nk)
    for v in range(1, N_VERSCHIEBUNGEN + 1):
        ck = {kombi: ordnungsstatistik(_verschieben(s, v), KS) for kombi, s in S.items()}
        for i, cfg, c in konfig_reihen(ck):
            g_ver[i] += np.isfinite(c).sum() * DT
            _, _, _, spitze = ereignisse(c, BODEN, VERBINDEN)
            if spitze.size:
                h_ver[i] += _hist(spitze)
        del ck
    erg.update(hist_verschoben=h_ver, gueltig_verschoben_s=g_ver)
    del S

    # 3. Künstliche Ausbrüche
    rng = np.random.default_rng([SEED, tag_index])
    inj = injektionen_planen(tag, alle, rng)
    zaehl = {b: tag.zaehl[b].copy() for b in TRIGGERBAENDER}
    namen = ["kurz", "lang"]
    for t, kl, amp, f in zip(inj["met"], inj["klasse"], inj["amplitude"], inj["faktor"]):
        if amp <= 0:
            continue
        j0 = int(np.floor((t - tag.t0) / DT))
        simulation.einbauen(zaehl, tag.dauer, j0, t - (tag.t0 + j0 * DT), simulation.KLASSEN[namen[kl]], amp, f,
                            rng)
    S = alle_signifikanzen(zaehl, tag.dauer)
    ck = {kombi: ordnungsstatistik(s, KS) for kombi, s in S.items()}
    del S
    lo = np.searchsorted(zeit, inj["met"] + INJ_FENSTER[0])
    hi = np.searchsorted(zeit, inj["met"] + INJ_FENSTER[1], side="right")
    imax = np.full((nk, inj["met"].size), -np.inf, np.float32)
    for i, cfg, c in konfig_reihen(ck):
        imax[i] = _fenster_max(c, lo, hi)
    erg["inj"] = {k: v for k, v in inj.items() if k != "faktor"}
    erg["inj"]["faktor_max"] = inj["faktor"].max(axis=1) if inj["met"].size else np.zeros(0)
    erg["inj_max"] = imax
    erg["rechenzeit_s"] = round(time.time() - t_start, 1)
    return erg
