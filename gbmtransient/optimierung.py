"""Optimierung des Triggers: höchste Nachweiswahrscheinlichkeit bei fester Fehlalarmrate (Neyman-Pearson).

Für jede der 162 Konfigurationen wird die Schwelle s* als kleinste Schwelle bestimmt, bei der auf den
Trainingstagen beide Nebenbedingungen gelten:
  - zufällige Fehlalarme (aus den Zeitverschiebungen) höchstens FAR_MAX je Tag,
  - unbestätigte Kandidaten (Nulllauf, nicht im Triggerkatalog) höchstens RU_MAX je Tag (Prüfaufwand).
Beide Raten fallen monoton mit der Schwelle; s* ist darum exakt (auf dem Raster der Histogramme, 0,02).
Zielgröße J = mittlere Nachweiswahrscheinlichkeit der künstlichen Ausbrüche (Stärken log-gleichverteilt, beide
Pulsklassen gleich gewichtet) bei s*. Gewählt wird die Konfiguration mit dem höchsten J auf den Trainingstagen
(vollständige Suche über das Gitter); bewertet wird auf den übrigen Tagen (Testtage).
"""

from __future__ import annotations

import numpy as np

from .simulation import KLASSEN, spitze_1024
from .statistik import logistik_anpassen, logistik_quantil, poisson_intervall, wilson
from .suche import KANTEN, KONFIGS

SPITZE_1024 = np.array([spitze_1024(KLASSEN["kurz"]), spitze_1024(KLASSEN["lang"])])

FAR_MAX = 0.1  # je Tag, Annahme: ein zufälliger Fehlalarm in zehn Tagen
RU_MAX = 2.0  # je Tag, Annahme: zwei zusätzliche Kandidaten je Tag sind prüfbar
TAG_S = 86400.0


def aufteilen(tage: list[dict]) -> tuple[list[int], list[int]]:
    """Training: gerade Tagesnummer im Zeitraum (1., 3., 5. ...), Test: ungerade."""
    train = [i for i, t in enumerate(tage) if t["tag_index"] % 2 == 0]
    test = [i for i, t in enumerate(tage) if t["tag_index"] % 2 == 1]
    return train, test


def summieren(tage: list[dict], auswahl, gewichte=None) -> dict:
    w = np.ones(len(auswahl)) if gewichte is None else np.asarray(gewichte, float)
    out = {}
    for key in ("hist_unbestaetigt", "hist_bestaetigt", "hist_verschoben"):
        out[key] = sum(wi * tage[i][key].astype(np.float64) for wi, i in zip(w, auswahl))
    for key in ("gueltig_s", "gueltig_verschoben_s"):
        out[key] = sum(wi * tage[i][key] for wi, i in zip(w, auswahl))
    return out


def rate_ab(hist: np.ndarray, gueltig_s) -> np.ndarray:
    """Rate je Tag der Ereignisse mit Spitze >= Kante (für jede untere Histogrammkante), je Konfiguration."""
    ab = np.cumsum(hist[..., ::-1], axis=-1)[..., ::-1]
    g = np.asarray(gueltig_s, float)[..., None] / TAG_S
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(g > 0, ab / g, np.inf)


def schwellen(agg: dict, far_max: float = FAR_MAX, ru_max: float = RU_MAX):
    """s* je Konfiguration und welche Nebenbedingung sie festlegt ('zufall', 'pruefaufwand'); inf, wenn keine."""
    rz = rate_ab(agg["hist_verschoben"], agg["gueltig_verschoben_s"])
    ru = rate_ab(agg["hist_unbestaetigt"], agg["gueltig_s"])
    kanten = KANTEN[:-1]

    def erste(ok):
        # erste Kante, ab der die Bedingung für alle höheren Kanten gilt
        schlecht = ~ok
        letzte_schlechte = np.where(schlecht.any(axis=-1), ok.shape[-1] - 1 - np.argmax(schlecht[..., ::-1], axis=-1),
                                    -1)
        return letzte_schlechte + 1

    iz, iu = erste(rz <= far_max), erste(ru <= ru_max)
    i = np.maximum(iz, iu)
    s = np.where(i < kanten.size, kanten[np.minimum(i, kanten.size - 1)], np.inf)
    bindend = np.where(iz >= iu, "zufall", "pruefaufwand")
    return s, bindend


def injektionen(tage: list[dict], auswahl) -> dict:
    """Künstliche Ausbrüche der gewählten Tage; 'staerke' = erwartete höchste Überschussrate über 1,024 s im
    hellsten Detektor, 50-300 keV (cts/s), dieselbe Größe, die für echte Ausbrüche gemessen wird."""
    teile = [tage[i]["inj"] for i in auswahl]
    out = {k: np.concatenate([t[k] for t in teile]) for k in ("met", "klasse", "amplitude", "theta", "phi",
                                                                 "faktor_max")}
    out["max"] = np.concatenate([tage[i]["inj_max"] for i in auswahl], axis=1)
    out["tag"] = np.concatenate([np.full(tage[i]["inj"]["met"].size, i) for i in auswahl])
    out["staerke"] = out["amplitude"] * out["faktor_max"] * SPITZE_1024[out["klasse"]]
    return out


def ziel(inj: dict, s) -> np.ndarray:
    """J je Konfiguration: Mittel der Nachweisquoten beider Pulsklassen (ohne Nullplätze)."""
    s = np.asarray(s, float)
    werte = []
    for kl in (0, 1):
        sel = (inj["klasse"] == kl) & (inj["amplitude"] > 0)
        werte.append((inj["max"][:, sel] >= s[:, None]).mean(axis=1))
    return 0.5 * (werte[0] + werte[1])


def waehlen(tage: list[dict], train, far_max=FAR_MAX, ru_max=RU_MAX) -> dict:
    agg = summieren(tage, train)
    s, bindend = schwellen(agg, far_max, ru_max)
    inj = injektionen(tage, train)
    j = ziel(inj, s)
    j = np.where(np.isfinite(s), j, 0.0)
    ru = rate_ab(agg["hist_unbestaetigt"], agg["gueltig_s"])
    ru_s = np.array([ru[i, np.searchsorted(KANTEN[:-1], s[i])] if np.isfinite(s[i]) else np.nan
                     for i in range(len(KONFIGS))])
    reihenfolge = np.lexsort((ru_s, -j))  # höchstes J, bei Gleichstand weniger Kandidaten
    return {"s": s, "bindend": bindend, "j": j, "ru": ru_s, "beste": int(reihenfolge[0]),
            "reihenfolge": reihenfolge}


def stabilitaet(tage: list[dict], train, n: int = 200, seed: int = 7, far_max=FAR_MAX, ru_max=RU_MAX) -> dict:
    """Bootstrap über die Trainingstage: wie oft wird dieselbe Konfiguration gewählt?"""
    rng = np.random.default_rng(seed)
    inj_alle = injektionen(tage, train)
    tag_pos = {t: k for k, t in enumerate(train)}
    gewinner = []
    for _ in range(n):
        zug = rng.integers(0, len(train), len(train))
        mult = np.bincount(zug, minlength=len(train))
        agg = summieren(tage, train, mult)
        s, _ = schwellen(agg, far_max, ru_max)
        gew = mult[[tag_pos[t] for t in inj_alle["tag"]]]
        werte = []
        for kl in (0, 1):
            sel = (inj_alle["klasse"] == kl) & (inj_alle["amplitude"] > 0) & (gew > 0)
            treffer = (inj_alle["max"][:, sel] >= s[:, None]) * gew[sel][None, :]
            werte.append(treffer.sum(axis=1) / gew[sel].sum())
        j = np.where(np.isfinite(s), 0.5 * (werte[0] + werte[1]), 0.0)
        gewinner.append(int(np.argmax(j)))
    werte, anzahl = np.unique(gewinner, return_counts=True)
    o = np.argsort(-anzahl)
    return {"ziehungen": n, "haeufigste": [{"konfig": int(werte[i]), "anteil": float(anzahl[i] / n)} for i in o[:5]]}


def bewerten(tage: list[dict], auswahl, i_cfg: int, s: float, trigger_typen: dict, n_boot: int = 1000,
             seed: int = 11) -> dict:
    """Kennzahlen einer Konfiguration mit fester Schwelle s auf den gewählten Tagen, mit Unsicherheiten."""
    agg = summieren(tage, auswahl)
    kz = int(round(agg["hist_verschoben"][i_cfg][np.searchsorted(KANTEN[:-1], s):].sum()))
    zeit_z = agg["gueltig_verschoben_s"][i_cfg] / TAG_S
    lo, hi = poisson_intervall(kz)
    # unbestätigte Kandidaten: Zahl je Tag, Bootstrap über Tage (Tage sind nicht gleich ruhig)
    b = np.searchsorted(KANTEN[:-1], s)
    je_tag = np.array([tage[i]["hist_unbestaetigt"][i_cfg][b:].sum() for i in auswahl], float)
    zeit_tag = np.array([tage[i]["gueltig_s"][i_cfg] for i in auswahl]) / TAG_S
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        z = rng.integers(0, len(auswahl), len(auswahl))
        boot.append(je_tag[z].sum() / zeit_tag[z].sum())
    inj = injektionen(tage, auswahl)
    klassen = {}
    for kl, name in ((0, "kurz"), (1, "lang")):
        sel = (inj["klasse"] == kl) & (inj["amplitude"] > 0)
        treffer = inj["max"][i_cfg, sel] >= s
        x = np.log10(inj["staerke"][sel])
        a, bb, cov = logistik_anpassen(x, treffer)
        q50, se50 = logistik_quantil(a, bb, cov, 0.5)
        q90, se90 = logistik_quantil(a, bb, cov, 0.9)
        wl, wh = wilson(treffer.sum(), treffer.size)
        klassen[name] = {"n": int(treffer.size), "erkannt": int(treffer.sum()), "quote": float(treffer.mean()),
                         "ki95": [float(wl), float(wh)],
                         "staerke50_cts_s": float(10 ** q50), "staerke50_ki95": [float(10 ** (q50 - 1.96 * se50)),
                                                                                float(10 ** (q50 + 1.96 * se50))],
                         "staerke90_cts_s": float(10 ** q90), "staerke90_ki95": [float(10 ** (q90 - 1.96 * se90)),
                                                                                float(10 ** (q90 + 1.96 * se90))],
                         "logistik": [a, bb]}
    null = inj["amplitude"] == 0
    katalog = {}
    for typ in sorted(set(trigger_typen.values())):
        tr, ab = [], 0
        for i in auswahl:
            for t, m in zip(tage[i]["trigger"], tage[i]["trigger_max"][i_cfg]):
                if t["typ"] == typ and t["abgedeckt"]:
                    ab += 1
                    tr.append(bool(m >= s))
        if ab:
            wl, wh = wilson(sum(tr), ab)
            katalog[typ] = {"abgedeckt": ab, "erkannt": int(sum(tr)), "quote": float(np.mean(tr)),
                            "ki95": [float(wl), float(wh)]}
    return {
        "schwelle": float(s), "tage": len(auswahl),
        "zufall": {"anzahl": kz, "tage_verschoben": float(zeit_z), "je_tag": float(kz / zeit_z),
                   "ki95": [float(lo / zeit_z), float(hi / zeit_z)]},
        "unbestaetigt": {"anzahl": int(je_tag.sum()), "je_tag": float(je_tag.sum() / zeit_tag.sum()),
                         "ki95_bootstrap": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                         "max_an_einem_tag": int(je_tag.max()) if je_tag.size else 0},
        "injektionen": klassen, "j": float(0.5 * (klassen["kurz"]["quote"] + klassen["lang"]["quote"])),
        "nullplaetze": {"n": int(null.sum()), "treffer": int((inj["max"][i_cfg, null] >= s).sum())},
        "katalog": katalog,
    }


def pareto(tage: list[dict], train, test, grund: int, budgets=(0.5, 1, 2, 5, 10, 20), far_max=FAR_MAX) -> list:
    """Für mehrere Prüfbudgets: beste Konfiguration (Training) und J auf den Testtagen, neben der Grundeinstellung."""
    out = []
    inj_test = injektionen(tage, test)
    for ru_max in budgets:
        w = waehlen(tage, train, far_max, ru_max)
        b = w["beste"]
        j_test = ziel(inj_test, w["s"])
        out.append({"budget_je_tag": ru_max, "beste": b, "schwelle": float(w["s"][b]),
                    "j_train": float(w["j"][b]), "j_test": float(j_test[b]),
                    "grund_schwelle": float(w["s"][grund]), "grund_j_test": float(j_test[grund])})
    return out
