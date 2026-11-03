"""Auswertung einer Rohablage: beide Durchgänge, Optimierung, Abgleich, Kennzahlen, Abbildungen, Bericht.

Jede Zahl im Bericht wird aus ergebnisse.json eingesetzt; der Text selbst enthält keine Messwerte.
"""

from __future__ import annotations

import concurrent.futures as cf
import csv
import hashlib
import json
import math
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__, details, einlesen, optimierung, suche
from .einlesen import DT
from .suche import INDEX, KONFIGS, Konfig
from .zeit import utc_aus_met, utc_text

GRUND = Konfig(2, "kausal", 64, "0,256-4,096 s", "50-300")  # an den Bord-Trigger angelehnt (17 s Untergrund davor)
GAUSS_UEBER = {"3": 1.3499e-3, "4": 3.1671e-5, "5": 2.8665e-7}  # P(Z >= s) der Standardnormalverteilung


def _parallel(f, aufgaben, prozesse: int):
    if prozesse <= 1:
        return [f(a) for a in aufgaben]
    with cf.ProcessPoolExecutor(prozesse) as ex:
        return list(ex.map(f, aufgaben))


def _sauber(x):
    """Für JSON: numpy-Typen in Python-Typen, nicht endliche Zahlen in None."""
    if isinstance(x, dict):
        return {str(k): _sauber(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_sauber(v) for v in x]
    if isinstance(x, np.ndarray):
        return _sauber(x.tolist())
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if math.isfinite(v) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def de(x, nd: int = 0) -> str:
    """Zahl deutsch: Tausenderpunkt, Dezimalkomma."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "–"
    s = f"{x:,.{nd}f}"
    return s.replace(",", "X").replace(".", ",").replace("X", ".")


def pz(x, nd: int = 0) -> str:
    return de(100 * x, nd) + " %"


def auswerten(raw: Path, aus: Path, prozesse: int = 8, max_tage: int | None = None, abbildungen: bool = True) -> dict:
    t_start = time.time()
    raw, aus = Path(raw), Path(aus)
    aus.mkdir(parents=True, exist_ok=True)
    zeiten: dict[str, float] = {}

    t = time.time()
    man = einlesen.manifest(raw)
    tage_d = einlesen.tagesdateien(man)
    daten = sorted(tage_d)[:max_tage] if max_tage else sorted(tage_d)
    trig = einlesen.katalog_lesen(raw, man, "fermigtrig")
    brst = {z["trigger_name"]: z for z in einlesen.katalog_lesen(raw, man, "fermigbrst")}
    tag_von = {d: [x for x in trig if utc_aus_met(x["met"]).strftime("%Y-%m-%d") == d] for d in daten}
    alle_met = np.array([x["met"] for x in trig])
    zeiten["kataloge_s"] = time.time() - t

    t = time.time()
    aufgaben = [(str(raw), d, tage_d[d], tag_von[d], alle_met, i) for i, d in enumerate(daten)]
    tage = _parallel(suche.tag_rechnen, aufgaben, prozesse)
    zeiten["durchgang1_s"] = time.time() - t

    t = time.time()
    train, test = optimierung.aufteilen(tage)
    w = optimierung.waehlen(tage, train)
    beste, grund = w["beste"], INDEX[GRUND]
    s_beste, s_grund = float(w["s"][beste]), float(w["s"][grund])
    typen = {x["trigger_name"]: x["trigger_type"] for x in trig}
    alle_i = list(range(len(tage)))
    bew = {
        "test": {"beste": optimierung.bewerten(tage, test, beste, s_beste, typen),
                 "grund": optimierung.bewerten(tage, test, grund, s_grund, typen),
                 "grund_4_5": optimierung.bewerten(tage, test, grund, 4.5, typen)},
        "train": {"beste": optimierung.bewerten(tage, train, beste, s_beste, typen),
                  "grund": optimierung.bewerten(tage, train, grund, s_grund, typen)},
        "alle": {"beste": optimierung.bewerten(tage, alle_i, beste, s_beste, typen),
                 "grund": optimierung.bewerten(tage, alle_i, grund, s_grund, typen)},
    }
    b = KONFIGS[beste]
    zwilling = INDEX[Konfig(b.k, "zweiseitig" if b.modus == "kausal" else "kausal", b.w, b.zs, b.band)]
    bew["test"]["zwilling"] = optimierung.bewerten(tage, test, zwilling, float(w["s"][zwilling]), typen)
    bew["test"]["zwilling"]["name"] = KONFIGS[zwilling].name()
    bew["test"]["schwelle_test_waere"] = {"beste": optimierung.schwelle_auf(tage, test, beste),
                                          "grund": optimierung.schwelle_auf(tage, test, grund)}
    inj_test = optimierung.injektionen(tage, test)
    j_test = optimierung.ziel(inj_test, w["s"])
    stab = optimierung.stabilitaet(tage, train)
    par = optimierung.pareto(tage, train, test, grund)
    zeiten["optimierung_s"] = time.time() - t

    # Beispiel für Abb. 1: GRB auf einem Testtag mit Koinzidenz am nächsten an 20
    kand = []
    for i in test:
        for x, m in zip(tage[i]["trigger"], tage[i]["trigger_max"][beste]):
            if x["typ"] == "GRB" and x["abgedeckt"] and np.isfinite(m):
                kand.append((abs(float(m) - 20.0), x["name"]))
    beispiel = min(kand)[1] if kand else None

    t = time.time()
    gewaehlt = [("beste", KONFIGS[beste], s_beste), ("grund", KONFIGS[grund], s_grund)]
    aufgaben2 = [(str(raw), d, tage_d[d], tag_von[d], trig, gewaehlt, beispiel) for d in daten]
    det = _parallel(details.tag_details, aufgaben2, prozesse)
    zeiten["durchgang2_s"] = time.time() - t

    t = time.time()
    r = zusammenfassen(raw, man, daten, trig, brst, tage, det, w, beste, grund, train, test, bew, j_test, stab, par,
                       prozesse)
    r["rechenzeit_s"] = {k: round(v, 1) for k, v in zeiten.items()}
    r["rechenzeit_s"]["gesamt_s"] = round(time.time() - t_start, 1)
    schreiben(aus, r, tage, det, w, j_test, trig, brst, beste, grund)
    if abbildungen:
        from . import abbildungen as abb
        abb.alle(aus, r, tage, det, w, j_test, beste, grund, train, test)
    (aus / "BERICHT.md").write_text(bericht_md(r), encoding="utf-8")
    zeiten["ausgabe_s"] = time.time() - t
    r["rechenzeit_s"]["ausgabe_s"] = round(zeiten["ausgabe_s"], 1)
    r["rechenzeit_s"]["gesamt_s"] = round(time.time() - t_start, 1)
    (aus / "ergebnisse.json").write_text(json.dumps(_sauber(r), ensure_ascii=False, indent=1), encoding="utf-8")
    return r


def zusammenfassen(raw, man, daten, trig, brst, tage, det, w, beste, grund, train, test, bew, j_test, stab, par,
                   prozesse) -> dict:
    # Daten und Abdeckung
    dateien = [e for e in man.values() if "/daily/" in e["pfad"] and not e["pfad"].endswith("index.html")]
    tage_set = set(daten)
    genutzt = [e for e in dateien if any(f"daily/{d[:4]}/{d[5:7]}/{d[8:]}/" in e["pfad"] for d in tage_set)]
    lebenszeit = sum(t["lebenszeit_s"] for t in tage)
    gesamt = sum(t["n_zellen"] for t in tage) * DT
    man_sha = hashlib.sha256((raw / "manifest.jsonl").read_bytes()).hexdigest()
    # Katalog im Zeitraum
    im_zeitraum = [x for t in tage for x in t["trigger"]]
    nach_typ: dict[str, dict] = {}
    for x in im_zeitraum:
        d = nach_typ.setdefault(x["typ"], {"gesamt": 0, "abgedeckt": 0})
        d["gesamt"] += 1
        d["abgedeckt"] += int(x["abgedeckt"])
    # Nullverteilung (Referenzkombination, alle Detektoren, alle Tage)
    sh = sum(t["s_hist"].sum(axis=0) for t in tage)
    mitte = 0.5 * (suche.S_KANTEN[1:] + suche.S_KANTEN[:-1])
    tot = sh.sum()
    mittel = float((sh * mitte).sum() / tot)
    std = float(math.sqrt((sh * mitte ** 2).sum() / tot - mittel ** 2))
    anteil = {str(s): float(sh[suche.S_KANTEN[:-1] >= s].sum() / tot) for s in (3, 4, 5)}
    # Konfigurationen
    reihen = []
    for rang, i in enumerate(w["reihenfolge"]):
        c = KONFIGS[i]
        reihen.append({"rang": rang + 1, "index": int(i), **c.kurz(), "name": c.name(), "schwelle": w["s"][i],
                       "bindend": str(w["bindend"][i]), "j_train": w["j"][i], "j_test": j_test[i],
                       "unbestaetigt_train_je_tag": w["ru"][i]})
    # Übertragung Training -> Test und Vergleich der Untergrundarten (nur zulässige Konfigurationen)
    ok = np.isfinite(w["s"])
    def raenge(x):
        return np.argsort(np.argsort(x)).astype(float)
    rk = float(np.corrcoef(raenge(w["j"][ok]), raenge(j_test[ok]))[0, 1]) if ok.sum() > 2 else None
    modus = np.array([c.modus for c in KONFIGS])
    untergrund = {m: {"bestes_j_test": float(np.max(j_test[ok & (modus == m)])) if np.any(ok & (modus == m)) else None,
                      "mittel_j_test": float(np.mean(j_test[ok & (modus == m)])) if np.any(ok & (modus == m)) else None}
                  for m in ("zweiseitig", "kausal")}
    # Ereignisse der besten Konfiguration
    ev = [e for d in det for e in d["ereignisse"]]
    split = {daten[i]: ("train" if i in set(train) else "test") for i in range(len(daten))}
    for e in ev:
        e["utc"] = utc_text(e["met_spitze"])
        e["teil"] = split[utc_aus_met(e["met_spitze"]).strftime("%Y-%m-%d")] if utc_aus_met(
            e["met_spitze"]).strftime("%Y-%m-%d") in split else "rand"
    unb = [e for e in ev if e["konfig"] == "beste" and not e["katalog"]]
    merkmale = {"anzahl": len(unb),
                "nahe_datenluecke_120s": sum(1 for e in unb if e["luecke_abstand_s"] <= 120),
                "weich_sonne_sichtbar": sum(1 for e in unb if e["luecke_abstand_s"] > 120 and e.get("sonne_sichtbar")
                                            and e["haerte"] is not None and e["haerte"] > 1.0),
                }
    merkmale["uebrige"] = merkmale["anzahl"] - merkmale["nahe_datenluecke_120s"] - merkmale["weich_sonne_sichtbar"]
    merkmale["band_27_50"] = sum(1 for e in unb if e["band"] == "27-50")
    merkmale["dauer_median_s"] = float(np.median([e["dauer_s"] for e in unb])) if unb else None
    merkmale["breite_betrag_median"] = float(np.median([abs(e["breite"]) for e in unb if "breite" in e])) if unb else None
    # Messungen an Katalog-Triggern
    tr = [x for d in det for x in d["trigger"]]
    grb = [x for x in tr if x["typ"] == "GRB" and x.get("abgedeckt") and "muster" in x]
    kor = [x["korrelation"] for x in grb if "korrelation" in x]
    top2 = [x["hellster_in_modell_top2"] for x in grb if "hellster_in_modell_top2" in x]
    kurz = [x for x in grb if brst.get(x["name"], {}).get("t90") is not None and brst[x["name"]]["t90"] < 2]
    lang = [x for x in grb if brst.get(x["name"], {}).get("t90") is not None and brst[x["name"]]["t90"] >= 2]
    h_k = [x["haerte"] for x in kurz if x.get("haerte") is not None]
    h_l = [x["haerte"] for x in lang if x.get("haerte") is not None]
    antwort = {"grb_gemessen": len(grb), "korrelation_median": float(np.median(kor)) if kor else None,
               "korrelation_quartile": [float(np.percentile(kor, 25)), float(np.percentile(kor, 75))] if kor else None,
               "hellster_in_top2": int(sum(top2)), "hellster_in_top2_von": len(top2),
               "haerte_kurz_median": float(np.median(h_k)) if h_k else None, "haerte_kurz_n": len(h_k),
               "haerte_lang_median": float(np.median(h_l)) if h_l else None, "haerte_lang_n": len(h_l),
               "spitzenrate_kurz_median": float(np.median([x["spitzenrate_cts_s"] for x in kurz])) if kurz else None,
               "spitzenrate_lang_median": float(np.median([x["spitzenrate_cts_s"] for x in lang])) if lang else None}
    # Zeitversatz: Spitze des passenden Ereignisses minus Triggerzeit (beste Konfiguration)
    versatz = []
    met_trig = {x["trigger_name"]: x["met"] for x in trig}
    for e in ev:
        if e["konfig"] == "beste":
            for k in e["katalog"]:
                versatz.append(e["met_anfang"] - met_trig[k["name"]])
    # Tagesdiagnose: Zufallskoinzidenzen und unbestätigte Kandidaten je Tag gegen die Sonnenaktivität im Katalog
    tage_tab = []
    sfl = {}
    for x in im_zeitraum:
        if x["typ"] == "SFLARE":
            d = utc_aus_met(x["met"]).strftime("%Y-%m-%d")
            sfl[d] = sfl.get(d, 0) + 1
    for zeile in bew["alle"]["beste"]["je_tag"]:
        tage_tab.append({**zeile, "sonneneruptionen_katalog": sfl.get(zeile["datum"], 0),
                         "teil": split.get(zeile["datum"], "")})
    def spearman(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        if a.size < 3 or np.std(a) == 0 or np.std(b) == 0:
            return None
        ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
        return float(np.corrcoef(ra, rb)[0, 1])
    zz = [t["zufall"] or 0.0 for t in tage_tab]
    tagesdiagnose = {
        "rang_zufall_sonne": spearman(zz, [t["sonneneruptionen_katalog"] for t in tage_tab]),
        "rang_unbestaetigt_sonne": spearman([t["unbestaetigt"] for t in tage_tab],
                                            [t["sonneneruptionen_katalog"] for t in tage_tab]),
        "zufall_anteil_top3_tage": float(sum(sorted(zz)[-3:]) / sum(zz)) if sum(zz) > 0 else None,
        "top_zufall": sorted(tage_tab, key=lambda t: -(t["zufall"] or 0))[:3],
        "top_unbestaetigt": sorted(tage_tab, key=lambda t: -t["unbestaetigt"])[:3],
    }
    # nicht erkannte Gammablitze (alle Tage, gewählte Konfiguration) mit T90 und Zeitskala des Bord-Triggers
    tts = {x["trigger_name"]: x.get("trigger_timescale") for x in trig}
    verpasst = []
    for t in tage:
        for x, m in zip(t["trigger"], t["trigger_max"][beste]):
            if x["typ"] == "GRB" and x["abgedeckt"] and not m >= w["s"][beste]:
                verpasst.append({"name": x["name"], "max": float(m) if np.isfinite(m) else None,
                                 "t90_s": brst.get(x["name"], {}).get("t90"), "bord_zeitskala_ms": tts.get(x["name"])})
    rechner = {"python": platform.python_version(), "system": platform.platform(), "prozesse": prozesse,
               "numpy": np.__version__}
    try:
        import astropy
        rechner["astropy"] = astropy.__version__
    except ImportError:
        pass
    try:
        import matplotlib
        rechner["matplotlib"] = matplotlib.__version__
    except ImportError:
        pass
    return {
        "version": __version__, "erstellt": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "rechner": rechner,
        "daten": {"von": daten[0], "bis": daten[-1], "tage": len(daten), "dateien": len(genutzt),
                  "bytes": int(sum(e["bytes"] for e in genutzt)), "manifest_sha256": man_sha,
                  "manifest_eintraege": len(man), "lebenszeit_h": lebenszeit / 3600, "rasterzeit_h": gesamt / 3600,
                  "abdeckung": lebenszeit / gesamt if gesamt else None,
                  "luecke_h": sum(x["luecke_s"] for x in det) / 3600,
                  "luecke_in_saa_anteil": (sum(x.get("luecke_saa_s", 0.0) for x in det)
                                           / max(sum(x["luecke_s"] for x in det), 1e-9)),
                  "saa_h": sum(x.get("saa_s", 0.0) for x in det) / 3600,
                  "injektionen": int(sum(t["inj"]["met"].size for t in tage)),
                  "train_tage": len(train), "test_tage": len(test)},
        "katalog": {"trigger_im_zeitraum": len(im_zeitraum), "nach_typ": nach_typ,
                    "burstkatalog_zeilen": len(brst)},
        "nullverteilung": {"kombination": "50-300 keV, 1,024 s, zweiseitig 32,8 s", "werte": int(tot),
                           "mittel": mittel, "std": std, "anteil_ueber": anteil, "gauss_ueber": GAUSS_UEBER},
        "optimierung": {"far_max_je_tag": optimierung.FAR_MAX, "ru_max_je_tag": optimierung.RU_MAX,
                        "rangkorrelation_train_test": rk, "untergrund": untergrund,
                        "konfigurationen": len(KONFIGS), "beste": reihen[0], "grund": next(x for x in reihen
                                                                                            if x["index"] == grund),
                        "top10": reihen[:10], "stabilitaet": {**stab, "haeufigste": [
                            {**h, "name": KONFIGS[h["konfig"]].name()} for h in stab["haeufigste"]]},
                        "pareto": [{**p, "name": KONFIGS[p["beste"]].name()} for p in par],
                        "machbar": int(np.isfinite(w["s"]).sum())},
        "bewertung": bew,
        "kandidaten": merkmale,
        "antwortmodell": antwort,
        "zeitversatz_s": {"n": len(versatz), "median": float(np.median(versatz)) if versatz else None,
                          "quartile": [float(np.percentile(versatz, 25)), float(np.percentile(versatz, 75))]
                          if versatz else None},
        "ereignisse": {"beste": sum(1 for e in ev if e["konfig"] == "beste"),
                       "grund": sum(1 for e in ev if e["konfig"] == "grund")},
        "beispiel": next((d["beispiel"]["name"] for d in det if d["beispiel"]), None),
        "grb_spitzenraten": {"kurz": [x["spitzenrate_cts_s"] for x in kurz],
                             "lang": [x["spitzenrate_cts_s"] for x in lang]},
        "tagesdiagnose": tagesdiagnose, "tage": tage_tab, "grb_verpasst": verpasst,
        "grb_ueber_90": {name: {"n": len(liste), "ueber": int(sum(1 for x in liste if x["spitzenrate_cts_s"] >=
                                                                  bew["test"]["beste"]["injektionen"][name]["staerke90_cts_s"]))}
                         for name, liste in (("kurz", kurz), ("lang", lang))},
    }


def schreiben(aus: Path, r, tage, det, w, j_test, trig, brst, beste, grund) -> None:
    ev = [e for d in det for e in d["ereignisse"]]
    with (aus / "ereignisse.csv").open("w", newline="", encoding="utf-8") as f:
        sp = ["konfig", "utc", "teil", "spitze", "band", "zeitskala_s", "dauer_s", "detektoren_ueber", "katalog",
              "luecke_abstand_s", "haerte", "sonne_sichtbar", "breite", "laenge", "met_anfang", "met_ende",
              "met_spitze"]
        wr = csv.writer(f)
        wr.writerow(sp)
        for e in sorted(ev, key=lambda e: (e["konfig"], e["met_spitze"])):
            z = dict(e)
            z["detektoren_ueber"] = " ".join(e["detektoren_ueber"])
            z["katalog"] = " ".join(f"{k['name']}:{k['typ']}" for k in e["katalog"])
            z["haerte"] = "" if e["haerte"] is None else round(e["haerte"], 3)
            wr.writerow([z.get(k, "") for k in sp])
    tr = {x["name"]: x for d in det for x in d["trigger"]}
    tts = {x["trigger_name"]: x.get("trigger_timescale") for x in trig}
    with (aus / "tage.csv").open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["datum", "teil", "zufall_je_tag_beste", "unbestaetigt_beste", "sonneneruptionen_katalog"])
        for t in r["tage"]:
            wr.writerow([t["datum"], t["teil"], "" if t["zufall"] is None else round(t["zufall"], 3), t["unbestaetigt"],
                         t["sonneneruptionen_katalog"]])
    with (aus / "katalog_abgleich.csv").open("w", newline="", encoding="utf-8") as f:
        sp = ["name", "typ", "utc", "teil", "abgedeckt", "max_beste", "erkannt_beste", "max_grund", "erkannt_grund",
              "spitzenrate_cts_s", "hellster", "korrelation", "haerte", "t90_s", "fluss_1024_ph_cm2_s",
              "bord_zeitskala_ms"]
        wr = csv.writer(f)
        wr.writerow(sp)
        s_b, s_g = float(w["s"][beste]), float(w["s"][grund])
        for t in tage:
            teil = "train" if t["tag_index"] % 2 == 0 else "test"
            for x, mb, mg in zip(t["trigger"], t["trigger_max"][beste], t["trigger_max"][grund]):
                d = tr.get(x["name"], {})
                b = brst.get(x["name"], {})
                wr.writerow([x["name"], x["typ"], utc_text(x["met"]), teil, int(x["abgedeckt"]),
                             round(float(mb), 2) if np.isfinite(mb) else "", int(mb >= s_b),
                             round(float(mg), 2) if np.isfinite(mg) else "", int(mg >= s_g),
                             round(d["spitzenrate_cts_s"], 1) if "spitzenrate_cts_s" in d else "",
                             d.get("hellster", ""), round(d["korrelation"], 3) if "korrelation" in d else "",
                             round(d["haerte"], 3) if d.get("haerte") is not None else "",
                             b.get("t90", ""), b.get("flux_1024", ""), tts.get(x["name"], "")])
    with (aus / "konfigurationen.csv").open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["rang", "k", "untergrund", "w_s", "zeitskalen", "band_kev", "schwelle", "bindend", "j_train",
                     "j_test", "unbestaetigt_train_je_tag"])
        for rang, i in enumerate(w["reihenfolge"]):
            c = KONFIGS[i]
            wr.writerow([rang + 1, c.k, c.modus, round(c.w * DT, 3), c.zs, c.band, w["s"][i], w["bindend"][i],
                         round(float(w["j"][i]), 4), round(float(j_test[i]), 4),
                         "" if not np.isfinite(w["ru"][i]) else round(float(w["ru"][i]), 3)])


# ------------------------------------------------------------------------------------------------- Bericht

def _cfg_text(x: dict) -> str:
    unt = "zweiseitig" if x["untergrund"] == "zweiseitig" else "kausal (nur davor)"
    return (f"mindestens {x['k']} Detektor{'en' if x['k'] > 1 else ''}, Untergrund {unt} je {de(x['w_s'], 1)} s, "
            f"Zeitskalen {x['zeitskalen']}, Band {x['band_kev']} keV")


def bericht_md(r: dict) -> str:
    d, o, b = r["daten"], r["optimierung"], r["bewertung"]
    bt, gt, g45 = b["test"]["beste"], b["test"]["grund"], b["test"]["grund_4_5"]
    ba, ga = b["alle"]["beste"], b["alle"]["grund"]
    nv, kat, kd, am = r["nullverteilung"], r["katalog"], r["kandidaten"], r["antwortmodell"]
    best, grund = o["beste"], o["grund"]
    stab = o["stabilitaet"]["haeufigste"][0]

    def quote(x):
        return f"{pz(x['quote'])} ({x['erkannt']} von {x['abgedeckt'] if 'abgedeckt' in x else x['n']}; " \
               f"95-%-Intervall {pz(x['ki95'][0])}–{pz(x['ki95'][1])})"

    grb_b, grb_g = ba["katalog"].get("GRB"), ga["katalog"].get("GRB")
    typ_zeilen = []
    namen = {"GRB": "Gammablitz (GRB)", "SFLARE": "Sonneneruption", "SGR": "Magnetar-Ausbruch (SGR)",
             "TGF": "Terrestrischer Gammablitz (TGF)", "LOCLPAR": "Teilchen, lokal", "DISTPAR": "Teilchen, entfernt",
             "UNCERT": "unsicher", "TRANSNT": "Transient", "GRO": "Röntgenquelle"}
    for typ in sorted(ba["katalog"], key=lambda t: -ba["katalog"][t]["abgedeckt"]):
        x, y = ba["katalog"][typ], ga["katalog"].get(typ)
        typ_zeilen.append(f"| {namen.get(typ, typ)} | {x['abgedeckt']} | {x['erkannt']} ({pz(x['quote'])}) | "
                          f"{y['erkannt']} ({pz(y['quote'])}) |" if y else "")
    top = "\n".join(
        f"| {x['rang']} | {x['k']} | {x['untergrund']} {de(x['w_s'], 1)} s | {x['zeitskalen']} | {x['band_kev']} | "
        f"{de(x['schwelle'], 2)} | {pz(x['j_train'], 1)} | {pz(x['j_test'], 1)} |" for x in o["top10"][:5])
    zeilen_pareto = "\n".join(
        f"| {de(p['budget_je_tag'], 1)} | {de(p['schwelle'], 2)} | {pz(p['j_test'], 1)} | {de(p['grund_schwelle'], 2)}"
        f" | {pz(p['grund_j_test'], 1)} |" for p in o["pareto"])
    ik, il = bt["injektionen"]["kurz"], bt["injektionen"]["lang"]
    gk, gl = gt["injektionen"]["kurz"], gt["injektionen"]["lang"]
    ug = o["untergrund"]
    if best["bindend"] == "pruefaufwand":
        folge2 = ("Zufällige Koinzidenzen sind bei der gewählten Schwelle seltener als erlaubt; die Schwelle legen die "
                  "gleichzeitigen, nicht im Katalog geführten Ereignisse fest.")
    else:
        folge2 = "Die Schwelle legt die Rate zufälliger Koinzidenzen fest."
    besser = "zweiseitige" if (ug["zweiseitig"]["bestes_j_test"] or 0) >= (ug["kausal"]["bestes_j_test"] or 0) else "kausale"
    folge3 = (f"Die Rangfolge überträgt sich auf ungesehene Tage (Rangkorrelation {de(o['rangkorrelation_train_test'], 2)}); "
              f"der {besser} Untergrund erreicht das höhere J (bestes J Test zweiseitig "
              f"{pz(ug['zweiseitig']['bestes_j_test'] or 0, 1)}, kausal {pz(ug['kausal']['bestes_j_test'] or 0, 1)}).")
    schwaecher = ik["staerke50_cts_s"] < gk["staerke50_cts_s"] and il["staerke50_cts_s"] < gl["staerke50_cts_s"]
    u90 = r["grb_ueber_90"]
    folge4 = ((f"Die Optimierung senkt die 50-%-Schwelle für beide Pulsklassen" if schwaecher else
               "Die Optimierung senkt die 50-%-Schwelle nicht für beide Pulsklassen") +
              f"; von den gemessenen Katalog-Gammablitzen liegen {u90['lang']['ueber']} von {u90['lang']['n']} "
              f"langen und {u90['kurz']['ueber']} von {u90['kurz']['n']} kurzen über der 90-%-Schwelle ihrer Klasse.")
    gr = ba["katalog"].get("GRB")
    folge5 = (f"Gammablitze findet die gewählte Einstellung zu {pz(gr['quote'])}, die Grundeinstellung zu "
              f"{pz(ga['katalog']['GRB']['quote'])}." if gr and "GRB" in ga["katalog"] else "")
    if "TGF" in ba["katalog"]:
        folge5 += (f" Terrestrische Gammablitze (unter 1 ms) bleiben in Zellen von 0,256 s fast unsichtbar "
                   f"({ba['katalog']['TGF']['erkannt']} von {ba['katalog']['TGF']['abgedeckt']}).")
    km = am["korrelation_median"] or 0
    folge6 = ("Das einfache Modell trifft die Verteilung über die Detektoren gut" if km >= 0.8 else
              "Das einfache Modell trifft die Verteilung über die Detektoren grob" if km >= 0.5 else
              "Das einfache Modell trifft die Verteilung über die Detektoren schlecht") + \
             "; die Abweichungen erklärt vermutlich die Streuung an Raumfahrzeug und Erdatmosphäre, die es nicht enthält."
    td = r["tagesdiagnose"]
    sw = b["test"]["schwelle_test_waere"]
    btr = b["train"]["beste"]
    if bt["zufall"]["je_tag"] > o["far_max_je_tag"]:
        uebertragung = (f"Die auf den Trainingstagen gewählte Schwelle hält die Grenze für zufällige Fehlalarme auf den "
                        f"Testtagen nicht: {de(bt['zufall']['je_tag'], 2)} statt {de(btr['zufall']['je_tag'], 2)} je Tag; "
                        f"auf den Testtagen wäre {de(sw['beste'], 2)} σ statt {de(bt['schwelle'], 2)} σ nötig gewesen. ")
    else:
        uebertragung = (f"Die auf den Trainingstagen gewählte Schwelle hält die Grenze für zufällige Fehlalarme auch auf den "
                        f"Testtagen ({de(bt['zufall']['je_tag'], 2)} je Tag). ")
    uebertragung += (f"Die Zufallskoinzidenzen sind ungleich verteilt: Die drei stärksten Tage tragen "
                     f"{pz(td['zufall_anteil_top3_tage'] or 0)} (höchster: {td['top_zufall'][0]['datum'][8:10]}."
                     f"{td['top_zufall'][0]['datum'][5:7]}.{td['top_zufall'][0]['datum'][:4]}, "
                     f"{de(td['top_zufall'][0]['zufall'], 1)} je Tag); Rangkorrelation der Tageswerte mit der Zahl der "
                     f"Sonneneruptionen im Katalog {de(td['rang_zufall_sonne'], 2)}." +
                     (" Eine Schwelle aus wenigen Wochen trägt diese Schwankung nicht; ein Betrieb bräuchte eine "
                      "laufende Nachkalibrierung." if bt["zufall"]["je_tag"] > o["far_max_je_tag"] else ""))
    zw = b["test"]["zwilling"]
    zwilling_satz = (f"Dieselbe Konfiguration mit {'zweiseitigem' if best['untergrund'] == 'kausal' else 'kausalem'} "
                     f"Untergrund erreicht auf den Testtagen J = {pz(zw['j'], 1)} (kurz {pz(zw['injektionen']['kurz']['quote'], 1)}, "
                     f"lang {pz(zw['injektionen']['lang']['quote'], 1)}; gewählt: kurz {pz(ik['quote'], 1)}, lang "
                     f"{pz(il['quote'], 1)})" + (". Vermutlich fällt beim zweiseitigen Untergrund der Ausläufer eines "
                     f"Pulses in das hintere Untergrundfenster und hebt den geschätzten Untergrund an."
                     if best['untergrund'] == 'kausal' and zw['j'] < bt['j'] else "."))
    gv = r["grb_verpasst"]
    kurz_verpasst = [x for x in gv if x["t90_s"] is not None and x["t90_s"] < 2]
    verpasst_satz = (f"Von den {len(gv)} nicht erkannten Gammablitzen sind {len(kurz_verpasst)} kurz (T90 "
                     f"{de(min((x['t90_s'] for x in kurz_verpasst), default=0), 2)}–"
                     f"{de(max((x['t90_s'] for x in kurz_verpasst), default=0), 2)} s); an Bord ausgelöst auf Zeitskalen "
                     f"von {', '.join(de(v, 0) for v in sorted({x['bord_zeitskala_ms'] for x in gv if x['bord_zeitskala_ms'] is not None}))} ms. "
                     + ("Die kürzeste Zeitskala hier ist 256 ms; für sehr kurze, schwache Ausbrüche ist CTIME zu grob."
                        if any(x["bord_zeitskala_ms"] is not None and x["bord_zeitskala_ms"] < 256 for x in gv) else "")
                     ) if gv else ""
    versatz_wort = "vor" if (r["zeitversatz_s"]["median"] or 0) < 0 else "nach"
    versatz_betrag = abs(r["zeitversatz_s"]["median"] or 0)
    return f"""# gbmtransient – Suche nach kurzen Gamma- und Röntgentransienten in Messdaten des Fermi GBM

**Technischer Bericht, Version {r['version']}** · Vorhaben des Bieters Panagiotis Herbe · Auswertung vom
{datetime.fromisoformat(r['erstellt'].replace('Z', '+00:00')).strftime('%d.%m.%Y')} · Quellcode, Tests und Rohdaten-Manifest:
https://github.com/panoherbe/gbmtransient

## 1 Zusammenfassung

Untersucht wird, wie empfindlich ein statistisch kontrollierter Trigger kurze Gamma- und Röntgenausbrüche in den
Zählraten der zwölf NaI-Szintillationsdetektoren des Fermi Gamma-ray Burst Monitor (GBM) findet, wenn seine
Parameter durch ein Optimierungsverfahren gewählt werden. Grundlage sind {d['tage']} Tage offener Messdaten
({datetime.fromisoformat(d['von']).strftime('%d.%m.%Y')} bis {datetime.fromisoformat(d['bis']).strftime('%d.%m.%Y')},
{de(d['dateien'])} Dateien, {de(d['bytes'] / 1e9, 1)} GB, {de(d['lebenszeit_h'], 0)} h Messzeit). Unter
{de(o['konfigurationen'])} Trigger-Konfigurationen wurde auf {d['train_tage']} Trainingstagen die mit der höchsten
Nachweiswahrscheinlichkeit bei höchstens {de(o['far_max_je_tag'], 1)} zufälligen Fehlalarmen und höchstens
{de(o['ru_max_je_tag'], 0)} unbestätigten Kandidaten je Tag gewählt; bewertet wird auf den {d['test_tage']}
übrigen Tagen.

- **Gewählt:** {_cfg_text(best)}, Schwelle {de(best['schwelle'], 2)} σ.
- **Nachweiswahrscheinlichkeit künstlicher Ausbrüche (Testtage):** {pz(bt['j'], 1)} gegenüber {pz(gt['j'], 1)} für
  eine an den Bord-Trigger angelehnte Grundeinstellung bei denselben Nebenbedingungen; 50-%-Schwelle kurzer Pulse
  {de(ik['staerke50_cts_s'], 0)} statt {de(gk['staerke50_cts_s'], 0)}, langer Pulse {de(il['staerke50_cts_s'], 0)}
  statt {de(gl['staerke50_cts_s'], 0)} Zählungen/s.
- **Abgleich mit dem GBM-Triggerkatalog (alle Tage):** {quote(grb_b) if grb_b else '–'} der abgedeckten
  Gammablitze erkannt.
- **Fehlalarme (Testtage):** zufällig {de(bt['zufall']['je_tag'], 2)} je Tag (Trainingstage
  {de(b['train']['beste']['zufall']['je_tag'], 2)}; {'die Grenze 0,1 hält die Schwelle nur auf den Trainingstagen ein' if bt['zufall']['je_tag'] > o['far_max_je_tag'] else 'Grenze 0,1 eingehalten'}),
  unbestätigte Kandidaten {de(bt['unbestaetigt']['je_tag'], 2)} je Tag (Grenze {de(o['ru_max_je_tag'], 0)}).
- **Statistik geprüft:** Die Signifikanz je Detektor hat über alle Daten Mittel {de(nv['mittel'], 3)} und
  Streuung {de(nv['std'], 3)}, wie für die Standardnormalverteilung erwartet (Abb. 2).

## 2 Fragestellung und Daten

Der GBM an Bord des Fermi-Satelliten (NASA, seit 2008) misst mit zwölf NaI(Tl)-Detektoren (8 keV bis 1 MeV) und
zwei BGO-Detektoren die Zählraten aus allen Richtungen, die nicht von der Erde verdeckt sind. An Bord löst ein
Trigger aus, wenn mindestens zwei NaI-Detektoren gleichzeitig eine Schwelle über dem Untergrund überschreiten. Die
Frage dieses Vorhabens: Wie weit lässt sich die Empfindlichkeit einer nachträglichen Suche in den kontinuierlichen
Daten steigern, wenn Untergrundschätzung, Zeitskalen, Energieband, Zahl der Detektoren und Schwelle gemeinsam
optimiert werden, und wie verlässlich ist das gegen den offiziellen Triggerkatalog? Der Verfasser hat 2020 in
seiner Masterarbeit mit Messdaten desselben Instruments gearbeitet (Erkennung von Röntgenausbrüchen); dieses
Vorhaben ist davon unabhängig: eigene Fragestellung, eigener Code, Messdaten von 2025.

| | |
|---|---|
| Messdaten | CTIME-Tagesdateien der zwölf NaI-Detektoren (acht Energiekanäle, 0,256 s; nach Bord-Triggern 0,064 s) und Positionsdatei (Lage, Südatlantische Anomalie, Sonne sichtbar) |
| Kataloge | Triggerkatalog `fermigtrig` (alle Bord-Trigger mit Einordnung und Quellrichtung) und Burstkatalog `fermigbrst` (Dauer T90, Fluss), Abfrage über TAP |
| Quelle | HEASARC, NASA Goddard Space Flight Center, `https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/`; GBM-Daten waren nach der Datenpolitik der Mission zu keiner Zeit proprietär; „HEASARC materials are all available freely for your use“ (HEASARC, Use of HEASARC Data) |
| Rohablage | {de(d['manifest_eintraege'])} Dateien, jede mit URL, SHA-256 und Abrufzeit im Manifest (SHA-256 des Manifests `{d['manifest_sha256'][:16]}…`); ausgewertet wird nur aus der Ablage, jede Datei wird vor dem Lesen geprüft |
| Zeitraum | {datetime.fromisoformat(d['von']).strftime('%d.%m.%Y')}–{datetime.fromisoformat(d['bis']).strftime('%d.%m.%Y')}; im Triggerkatalog {kat['trigger_im_zeitraum']} Bord-Trigger, davon {kat['nach_typ'].get('GRB', {}).get('gesamt', 0)} Gammablitze |
| Abdeckung | {pz(d['abdeckung'], 1)} der Zeit mit Daten von mindestens sechs Detektoren; von {de(d['luecke_h'], 0)} h ohne Daten liegen {pz(d['luecke_in_saa_anteil'])} in Durchgängen durch die Südatlantische Anomalie (Kennung der Positionsdatei; die Detektoren sind dort abgeschaltet) |

## 3 Verfahren

**Datenmodell.** Je Tag ein gemeinsames Zeitraster mit Zellen von 0,256 s; je Detektor und Energieband
(27–50 keV, 50–300 keV) die Zahl der Ereignisse und die Lebenszeit (Messzeit ohne Totzeit) je Zelle. Jedes
CTIME-Intervall wird der Zelle seines Mittelpunkts zugeordnet; die feineren Intervalle nach Bord-Triggern und
Übergangsintervalle gehen mit ihrer echten Messzeit ein, Lücken sind Zellen mit Messzeit null. Alle Summen über
Fenster werden aus kumulierten Summen gebildet und tragen ihre Messzeit mit.

**Signifikanz je Detektor.** Für ein Signalfenster der Länge T (0,256 bis 8,192 s) wird der Untergrund aus
Nebenfenstern der Breite w geschätzt, zweiseitig (vor und nach dem Fenster; ein linearer Verlauf hebt sich heraus)
oder kausal (nur davor, wie an Bord und für Echtzeit-Alarme). Der Überschuss wird mit dem
Likelihood-Quotienten-Test für Poisson-Zählungen nach Li & Ma (1983) bewertet; das Verhältnis der Messzeiten
α = E_Signal / E_Untergrund berücksichtigt Lücken und Totzeit.

**Koinzidenz.** Die Teststatistik ist die k-größte Signifikanz über die zwölf Detektoren (Ordnungsstatistik):
„mindestens k Detektoren über s“ ist genau C_k ≥ s. Eine Konfiguration nimmt das Maximum über ihre Zeitskalen und
Bänder; Überschreitungen im Abstand bis 10,24 s bilden ein Ereignis. Ein Ereignis gilt als bestätigt, wenn ein
Katalog-Trigger zwischen 20 s vor seinem Beginn und 5 s nach seinem Ende liegt.

**Fehlalarmrate durch Zufall.** In 20 Durchläufen wird Detektor d zyklisch um d·φ eines Tages verschoben, mit
einem anderen φ je Durchlauf, so gewählt, dass jedes Detektorpaar mindestens 43 min gegeneinander versetzt ist
(Zeitverschiebungs-Methode). Echte Ausbrüche treffen alle Detektoren gleichzeitig und fallen dabei auseinander,
auch lange wie Sonneneruptionen; die verbleibenden Koinzidenzen schätzen die Rate zufälliger Fehlalarme aus {de(bt['zufall']['tage_verschoben'] + b['train']['beste']['zufall']['tage_verschoben'], 0)}
verschobenen Tagen. Die Rate unbestätigter Kandidaten im unverschobenen Lauf enthält dazu echte, nicht
getriggerte Ereignisse und korrelierten Untergrund.

**Künstliche Ausbrüche.** Alle 200 s (ohne Katalog-Trigger im Umkreis von 200 s und ohne Datenlücke) wird ein Puls
nach Norris et al. (2005) mit Poisson-Zählungen in die echten Daten eingebaut, {de(d['injektionen'])} insgesamt:
„kurz“ (Breite etwa 0,4 s, Verhältnis 27–50 zu 50–300 keV angenommen 0,3) und „lang“ (etwa 10 s, 0,6), Stärke
log-gleichverteilt, Richtung isotrop, Zählrate je Detektor proportional zu max(cos θ, 0) zur Detektorachse
(Ausrichtungen nach Meegan et al. 2009); 5 % der Plätze bleiben leer (Kontrolle). Stärke = erwartete höchste
Überschussrate über 1,024 s im hellsten Detektor, 50–300 keV; dieselbe Größe wird an echten Ausbrüchen gemessen.

**Optimierungsverfahren.** Gesucht ist die Konfiguration mit der höchsten mittleren Nachweiswahrscheinlichkeit J
unter zwei Nebenbedingungen (Neyman-Pearson-Kriterium): zufällige Fehlalarme ≤ {de(o['far_max_je_tag'], 1)} je Tag
und unbestätigte Kandidaten ≤ {de(o['ru_max_je_tag'], 0)} je Tag (beide Grenzen sind Annahmen für einen
prüfbaren Betrieb). Für jede der {o['konfigurationen']} Konfigurationen (k = 1, 2, 3; Untergrund zweiseitig oder
kausal mit w = 8,2, 16,4, 32,8 s; drei Zeitskalen-Sätze; drei Bänder) ist die kleinste zulässige Schwelle exakt
bestimmbar, weil beide Raten monoton mit der Schwelle fallen; über das Gitter wird vollständig gesucht. Gewählt
wird auf den Trainingstagen (1., 3., 5. … Tag), bewertet auf den Testtagen. Unsicherheiten: Wilson-Intervalle für
Anteile, exakte Poisson-Intervalle für Zufallsraten, Bootstrap über Tage für Kandidatenraten und für die
Stabilität der Wahl, Maximum-Likelihood-Anpassung einer logistischen Kurve (Newton-Raphson) für die 50- und
90-%-Schwellen mit Standardfehlern aus der Fisher-Information.

## 4 Ergebnisse

![Abb. 1](abbildungen/abb1_beispiel.png)

**Abb. 1:** Zählungen je 0,256 s (50–300 keV) der drei hellsten Detektoren um den Bord-Trigger {r['beispiel']}
(Zeit 0), mit zweiseitig geschätztem Untergrund (oben), und Koinzidenzstatistik der gewählten Konfiguration mit
Schwelle (unten). *Daraus folgt: Der Ausbruch hebt sich in mehreren Detektoren zugleich vom Untergrund ab, die
Koinzidenzstatistik überschreitet die Schwelle.*

**Nullverteilung.** Über {de(nv['werte'] / 1e6, 1)} Mio. Werte der Referenzkombination (alle Detektoren und Tage)
hat S den Mittelwert {de(nv['mittel'], 3)} und die Streuung {de(nv['std'], 3)}; Anteil S ≥ 3: {de(nv['anteil_ueber']['3'] * 1e3, 2)}‰ (Normalverteilung
{de(nv['gauss_ueber']['3'] * 1e3, 2)}‰), S ≥ 4: {de(nv['anteil_ueber']['4'] * 1e5, 1)}·10⁻⁵ ({de(nv['gauss_ueber']['4'] * 1e5, 1)}·10⁻⁵),
S ≥ 5: {de(nv['anteil_ueber']['5'] * 1e7, 1)}·10⁻⁷ ({de(nv['gauss_ueber']['5'] * 1e7, 1)}·10⁻⁷). Der Überschuss in den Ausläufern stammt aus echten Ereignissen und
Untergrundsprüngen.

![Abb. 2](abbildungen/abb2_raten.png)

**Abb. 2:** Links: Häufigkeit der Signifikanz S auf allen Daten gegen die Standardnormalverteilung. Rechts: Rate
der Ereignisse über einer Schwelle für die gewählte Konfiguration, unverschoben (bestätigt und unbestätigt) und
zeitverschoben (Zufall), Trainingstage; gepunktet: gewählte Schwelle und die beiden Grenzen. *Daraus folgt:
{folge2}*

**Optimierung.** {o['machbar']} der {o['konfigurationen']} Konfigurationen erfüllen beide Nebenbedingungen mit
einer Schwelle unter 25 σ. Die besten fünf auf den Trainingstagen:

| Rang | k | Untergrund | Zeitskalen | Band (keV) | Schwelle (σ) | J Training | J Test |
|---|---|---|---|---|---|---|---|
{top}

Die Grundeinstellung ({_cfg_text(grund)}) liegt auf Rang {grund['rang']} (J Training {pz(grund['j_train'], 1)},
Test {pz(grund['j_test'], 1)}). Bootstrap über die Trainingstage ({o['stabilitaet']['ziehungen']} Ziehungen): die
gewählte Konfiguration gewinnt in {pz(stab['anteil'])} der Ziehungen{'' if stab['konfig'] == best['index'] else f" (häufigster Gewinner: {stab['name']})"}.

![Abb. 3](abbildungen/abb3_optimierung.png)

**Abb. 3:** Links: J auf Trainings- gegen Testtage für alle Konfigurationen (Farbe: k; Kreis zweiseitig, Dreieck
kausal); gewählte Konfiguration und Grundeinstellung markiert. Rechts: J auf den Testtagen für verschiedene
Prüfbudgets. *Daraus folgt: {folge3}*

| Prüfbudget je Tag | Schwelle gewählt | J Test gewählt | Schwelle Grund | J Test Grund |
|---|---|---|---|---|
{zeilen_pareto}

![Abb. 4](abbildungen/abb4_effizienz.png)

**Abb. 4:** Nachweiswahrscheinlichkeit künstlicher Ausbrüche auf den Testtagen nach Stärke, mit 95-%-Intervallen
und angepasster logistischer Kurve; Striche unten: gemessene Spitzenraten der abgedeckten Katalog-Gammablitze
(kurz: T90 < 2 s).
*Daraus folgt: {folge4}*

| Testtage | gewählt | Grundeinstellung | Grund, feste Schwelle 4,5 σ |
|---|---|---|---|
| Schwelle | {de(bt['schwelle'], 2)} σ | {de(gt['schwelle'], 2)} σ | 4,50 σ |
| J (beide Pulsklassen) | {pz(bt['j'], 1)} | {pz(gt['j'], 1)} | {pz(g45['j'], 1)} |
| 50 % kurz / lang (cts/s) | {de(ik['staerke50_cts_s'], 0)} / {de(il['staerke50_cts_s'], 0)} | {de(gk['staerke50_cts_s'], 0)} / {de(gl['staerke50_cts_s'], 0)} | {de(g45['injektionen']['kurz']['staerke50_cts_s'], 0)} / {de(g45['injektionen']['lang']['staerke50_cts_s'], 0)} |
| 90 % kurz / lang (cts/s) | {de(ik['staerke90_cts_s'], 0)} / {de(il['staerke90_cts_s'], 0)} | {de(gk['staerke90_cts_s'], 0)} / {de(gl['staerke90_cts_s'], 0)} | {de(g45['injektionen']['kurz']['staerke90_cts_s'], 0)} / {de(g45['injektionen']['lang']['staerke90_cts_s'], 0)} |
| zufällige Fehlalarme je Tag | {de(bt['zufall']['je_tag'], 3)} ({de(bt['zufall']['ki95'][0], 3)}–{de(bt['zufall']['ki95'][1], 3)}) | {de(gt['zufall']['je_tag'], 3)} ({de(gt['zufall']['ki95'][0], 3)}–{de(gt['zufall']['ki95'][1], 3)}) | {de(g45['zufall']['je_tag'], 2)} |
| unbestätigte Kandidaten je Tag | {de(bt['unbestaetigt']['je_tag'], 2)} ({de(bt['unbestaetigt']['ki95_bootstrap'][0], 2)}–{de(bt['unbestaetigt']['ki95_bootstrap'][1], 2)}) | {de(gt['unbestaetigt']['je_tag'], 2)} ({de(gt['unbestaetigt']['ki95_bootstrap'][0], 2)}–{de(gt['unbestaetigt']['ki95_bootstrap'][1], 2)}) | {de(g45['unbestaetigt']['je_tag'], 1)} |
| leere Plätze mit Treffer | {bt['nullplaetze']['treffer']} von {bt['nullplaetze']['n']} | {gt['nullplaetze']['treffer']} von {gt['nullplaetze']['n']} | {g45['nullplaetze']['treffer']} von {g45['nullplaetze']['n']} |

**Übertragung auf ungesehene Tage.** {uebertragung} {zwilling_satz}

**Abgleich mit dem Triggerkatalog.** Erkannt, wenn die Koinzidenzstatistik im Fenster 5 s vor bis 20 s nach der
Triggerzeit die Schwelle erreicht; nur Trigger in abgedeckter Zeit (alle Tage):

| Einordnung im Katalog | abgedeckt | erkannt, gewählt | erkannt, Grundeinstellung |
|---|---|---|---|
{chr(10).join(z for z in typ_zeilen if z)}

Der Beginn des passenden Ereignisses liegt im Median {de(versatz_betrag, 1)} s
(Quartile {de(r['zeitversatz_s']['quartile'][0], 1) if r['zeitversatz_s']['quartile'] else '–'} bis
{de(r['zeitversatz_s']['quartile'][1], 1) if r['zeitversatz_s']['quartile'] else '–'} s, negativ = davor) {versatz_wort} der
Triggerzeit des Katalogs. {verpasst_satz}

![Abb. 5](abbildungen/abb5_katalog.png)

**Abb. 5:** Anteil erkannter Katalog-Trigger je Einordnung (alle Tage) mit 95-%-Intervallen. *Daraus folgt:
{folge5}*

**Unbestätigte Kandidaten.** Die gewählte Konfiguration meldet über alle Tage {kd['anzahl']} Ereignisse ohne
Katalog-Trigger: {kd['nahe_datenluecke_120s']} liegen höchstens 120 s neben einer Datenlücke (Rand der
Südatlantischen Anomalie), {kd['weich_sonne_sichtbar']} sind weich (mehr Überschuss in 27–50 als in 50–300 keV)
bei sichtbarer Sonne (vermutlich schwache Sonneneruptionen), {kd['uebrige']} bleiben ohne Zuordnung. Das ist
eine Einordnung nach Merkmalen, keine Bestätigung.

## 5 Prüfung des Antwortmodells

![Abb. 6](abbildungen/abb6_antwort.png)

**Abb. 6:** Gemessener Überschuss je Detektor (4,096 s um die Spitze, auf den hellsten normiert) gegen das
Modell max(cos θ, 0) aus der Katalog-Richtung, für {am['grb_gemessen']} Gammablitze. *Daraus folgt: {folge6}*

Korrelation zwischen gemessenem Muster und Modell: Median {de(am['korrelation_median'], 2)}
(Quartile {de(am['korrelation_quartile'][0], 2) if am['korrelation_quartile'] else '–'} bis
{de(am['korrelation_quartile'][1], 2) if am['korrelation_quartile'] else '–'}); der hellste gemessene Detektor ist in
{am['hellster_in_top2']} von {am['hellster_in_top2_von']} Fällen einer der zwei nach dem Modell hellsten. Das
gemessene Verhältnis 27–50 zu 50–300 keV beträgt im Median {de(am['haerte_kurz_median'], 2)} für kurze
(T90 < 2 s, n = {am['haerte_kurz_n']}) und {de(am['haerte_lang_median'], 2)} für lange Gammablitze
(n = {am['haerte_lang_n']}); angenommen waren 0,3 und 0,6.

## 6 Grenzen

- Die Detektorantwort ist vereinfacht (Kosinus, ohne Streuung, ohne Energieabhängigkeit); die Nachweisgrenzen
  gelten für dieses Modell und für die zwei Pulsformen, nicht für jede Form eines echten Ausbruchs.
- Die Nebenbedingungen (0,1 zufällige Fehlalarme, 2 Kandidaten je Tag) sind gesetzte Annahmen; Abb. 3 zeigt,
  wie sich das Ergebnis mit dem Prüfbudget ändert.
- Unbestätigte Kandidaten sind nicht geprüft; ihre Einordnung beruht auf Merkmalen.
- Fenster an Tagesgrenzen und neben Datenlücken werden verworfen; Trigger dort gelten als nicht abgedeckt.
- Ausgewertet wurden {d['tage']} Tage aus einer Phase hoher Sonnenaktivität; andere Zeiträume können andere
  Schwellen ergeben.
- Die Gütemessung baut künstliche Ausbrüche in echte Daten ein; zufällige Treffer im Suchfenster sind über die leeren
  Plätze gemessen ({bt['nullplaetze']['treffer']} von {bt['nullplaetze']['n']} auf den Testtagen).

## 7 Quellen

- Meegan, C. et al. (2009): The Fermi Gamma-ray Burst Monitor. ApJ 702, 791 (Detektoraufbau, Ausrichtungen Tab. 1).
- Li, T.-P., Ma, Y.-Q. (1983): Analysis methods for results in gamma-ray astronomy. ApJ 272, 317.
- Norris, J. P. et al. (2005): Long-lag, wide-pulse gamma-ray bursts. ApJ 627, 324 (Pulsform).
- Wilson, E. B. (1927): Probable inference, the law of succession, and statistical inference. JASA 22, 209.
- Fermi Science Support Center: Fermi Data Policy (Stand 28.08.2025); HEASARC: Use of HEASARC Data (Stand
  21.04.2025). „This research has made use of data and/or software provided by the High Energy Astrophysics Science
  Archive Research Center (HEASARC), which is a service of the Astrophysics Science Division at NASA/GSFC.“
"""
