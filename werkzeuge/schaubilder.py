"""Schaubilder für die README und den Rundgang, aus der Rohablage und den Ergebnissen eines Laufs.

    python werkzeuge/schaubilder.py --raw daten/raw/2026-11-03 --nachweis nachweis/2026-11-03 --aus medien \
        --bilder <ordner-für-animationsbilder>

Erzeugt in --aus: himmel.png (die Gammablitze des Katalogs am Himmel, gefunden oder nicht), detektoren.png (ein
Ausbruch in allen zwölf Detektoren), treffer.png (vom Rohsignal zum Treffer), orbit.png (ein Tag Messdaten mit
Bahn und Südatlantischer Anomalie). In --bilder: Einzelbilder für zwei Animationen des Rundgangs. Jede Zahl kommt
aus ergebnisse.json, katalog_abgleich.csv, ereignisse.csv oder unmittelbar aus den Rohdaten.
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from gbmtransient import einlesen, simulation
from gbmtransient.einlesen import DT, NAI
from gbmtransient.statistik import kumulieren, ordnungsstatistik, ueberschuss
from gbmtransient.suche import G, RAND, alle_signifikanzen
from gbmtransient.zeit import utc_aus_met

GRUND, FLAECHE, TEXT, TEXT2, GITTER = "#0b1220", "#111a2e", "#e2e8f0", "#94a3b8", "#1f2a44"
GOLD, BLAU, GRAU, ORANGE = "#fbbf24", "#60a5fa", "#64748b", "#f97316"
BEISPIEL = "bn251022940"
FENSTER = (-30.0, 130.0)  # s um die Triggerzeit des Beispiels (T90 laut Burstkatalog 99,6 s)


def stil():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": GRUND, "axes.facecolor": FLAECHE, "savefig.facecolor": GRUND, "text.color": TEXT,
        "axes.labelcolor": TEXT, "axes.edgecolor": GITTER, "xtick.color": TEXT2, "ytick.color": TEXT2,
        "axes.grid": True, "grid.color": GITTER, "grid.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "legend.frameon": False, "savefig.dpi": 150,
    })
    return plt


def de(x, nd=0):
    s = f"{x:,.{nd}f}"
    return s.replace(",", "X").replace(".", ",").replace("X", ".").replace("-", "−")


def komma(ax, achse="both"):
    from matplotlib.ticker import FuncFormatter

    f = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−"))
    if achse in ("x", "both"):
        ax.xaxis.set_major_formatter(f)
    if achse in ("y", "both"):
        ax.yaxis.set_major_formatter(f)


def lade(raw: Path, nachweis: Path) -> dict:
    r = json.loads((nachweis / "ergebnisse.json").read_text(encoding="utf-8"))
    man = einlesen.manifest(raw)
    trig = {t["trigger_name"]: t for t in einlesen.katalog_lesen(raw, man, "fermigtrig")}
    abgleich = list(csv.DictReader((nachweis / "katalog_abgleich.csv").open(encoding="utf-8")))
    ereignisse = list(csv.DictReader((nachweis / "ereignisse.csv").open(encoding="utf-8")))
    grb = []
    for z in abgleich:
        if z["typ"] == "GRB" and z["abgedeckt"] == "1":
            t = trig[z["name"]]
            grb.append({"name": z["name"], "met": t["met"], "ra": float(t["ra"]), "dec": float(t["dec"]),
                        "erkannt": z["erkannt_beste"] == "1", "utc": z["utc"]})
    grb.sort(key=lambda g: g["met"])
    return {"r": r, "man": man, "trig": trig, "grb": grb, "ereignisse": ereignisse,
            "abgleich": abgleich}


# ------------------------------------------------------------------------------------------------- Himmelskarte

def _mollweide_xy(ra, dec):
    """Rektaszension nach links steigend, 0 h in der Mitte (astronomische Konvention)."""
    ra = np.asarray(ra, float)
    x = -np.radians(((ra + 180.0) % 360.0) - 180.0)
    return x, np.radians(np.asarray(dec, float))


def himmel(plt, d: dict, ziel: Path, anzahl: int | None = None, groesse=(12.8, 7.2), titel=True):
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    grb = d["grb"] if anzahl is None else d["grb"][:anzahl]
    fig = plt.figure(figsize=groesse)
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_facecolor(FLAECHE)
    ax.grid(True, color=GITTER, linewidth=0.7)
    ax.set_xticklabels([])
    ax.tick_params(axis="y", colors=TEXT2, labelsize=9)
    # Milchstraße als Band: galaktische Ebene ±10°
    for b, a, lw in ((0, 0.55, 2.2), (10, 0.18, 1.2), (-10, 0.18, 1.2)):
        g = SkyCoord(l=np.linspace(0, 360, 1441) * u.deg, b=np.full(1441, b) * u.deg, frame="galactic").icrs
        x, y = _mollweide_xy(g.ra.deg, g.dec.deg)
        sprung = np.flatnonzero(np.abs(np.diff(x)) > 1) + 1
        for teil_x, teil_y in zip(np.split(x, sprung), np.split(y, sprung)):
            ax.plot(teil_x, teil_y, color="#8b5cf6", alpha=a, linewidth=lw, zorder=1)
    erk = [g for g in grb if g["erkannt"]]
    nicht = [g for g in grb if not g["erkannt"]]
    if erk:
        x, y = _mollweide_xy([g["ra"] for g in erk], [g["dec"] for g in erk])
        ax.scatter(x, y, s=150, color=GOLD, alpha=0.18, linewidths=0, zorder=2)
        ax.scatter(x, y, s=38, color=GOLD, edgecolors="#fff7d6", linewidths=0.6, zorder=3, label="gefunden")
    if nicht:
        x, y = _mollweide_xy([g["ra"] for g in nicht], [g["dec"] for g in nicht])
        ax.scatter(x, y, s=46, facecolors="none", edgecolors=TEXT2, linewidths=1.4, zorder=3, label="nicht gefunden")
    b = next((g for g in grb if g["name"] == BEISPIEL), None)
    if b:
        x, y = _mollweide_xy([b["ra"]], [b["dec"]])
        ax.scatter(x, y, s=420, facecolors="none", edgecolors=GOLD, linewidths=1.6, zorder=4)
        ax.annotate("22.10.2025, 22:34 UTC", (x[0], y[0]), xytext=(18, 14), textcoords="offset points",
                    color=GOLD, fontsize=10, fontweight="bold")
    alle = d["grb"]
    n_erk = sum(g["erkannt"] for g in grb)
    if titel:
        fig.text(0.06, 0.93, "Gammablitze am Himmel, 01.09.–31.10.2025", fontsize=17, fontweight="bold", color=TEXT)
        t = d["r"]["bewertung"]["test"]["beste"]["katalog"]["GRB"]  # Prüftage, nicht zur Wahl benutzt
        zeile = (f"{len(alle)} Gammablitze im offiziellen Katalog · {n_erk} gefunden, auf den Prüftagen "
                 f"{t['erkannt']} von {t['abgedeckt']}" if anzahl is None else
                 f"{len(grb)} von {len(alle)} Gammablitzen, zeitlich geordnet · {n_erk} gefunden")
        fig.text(0.06, 0.885, zeile + " · Band: Milchstraße", fontsize=12, color=TEXT2)
    leg = ax.legend(loc="lower right", bbox_to_anchor=(1.02, -0.08), fontsize=10, labelcolor=TEXT)
    for t in leg.get_texts():
        t.set_color(TEXT)
    fig.savefig(ziel, dpi=125, facecolor=GRUND if titel else "none")  # ohne Titel: für das Titelbild, Grund durchsichtig
    plt.close(fig)


# ------------------------------------------------------------------------------------------- ein Ausbruch, Daten

def ausbruch(d: dict, raw: Path) -> dict:
    """Zählungen, Untergrund und Signifikanzen um den Beispiel-Ausbruch; gewählte Konfiguration aus ergebnisse.json."""
    t = d["trig"][BEISPIEL]
    datum = utc_aus_met(t["met"]).strftime("%Y-%m-%d")
    tage = einlesen.tagesdateien(d["man"])
    tag = einlesen.tag_lesen(raw, datum, tage[datum], mit_position=True)
    beste = d["r"]["optimierung"]["beste"]
    w = int(round(beste["w_s"] / DT))
    bis = {"0,256-1,024 s": 4, "0,256-4,096 s": 16, "0,256-8,192 s": 32}[beste["zeitskalen"]]
    zs = [m for m in (1, 2, 4, 8, 16, 32) if m <= bis]
    kombis = {(beste["band_kev"], m, beste["untergrund"], w) for m in zs}
    S = alle_signifikanzen(tag.zaehl, tag.dauer, kombis)
    c = np.max(np.stack([ordnungsstatistik(S[k], (beste["k"],))[beste["k"]] for k in sorted(kombis)]), axis=0)
    zeit = tag.t0 + DT * np.arange(tag.n + 1)
    lo = int(np.searchsorted(zeit, t["met"] + FENSTER[0]))
    hi = int(np.searchsorted(zeit, t["met"] + FENSTER[1]))
    voll = float(np.median(tag.dauer[tag.dauer > 0]))
    cz = kumulieren(tag.zaehl["50-300"], RAND)
    ce = kumulieren(tag.dauer, RAND)
    x1, _ = ueberschuss(cz, ce, tag.n, RAND, 1, 128, G, voll)  # zweiseitiger Untergrund je Zelle (nur zur Anzeige)
    zell = tag.zaehl["50-300"][:, lo:hi].astype(float)
    untergrund = zell - x1[:, lo + 1:hi + 1]
    winkel = np.degrees(np.arccos(np.clip(simulation.ACHSEN @ simulation.einheitsvektor(float(t["phi"]),
                                                                                       float(t["theta"])), -1, 1)))
    s4 = S[(beste["band_kev"], 4, beste["untergrund"], w)][:, lo:hi]
    return {"t_rel": zeit[lo:hi] - t["met"], "zaehl": zell, "untergrund": untergrund, "winkel": winkel,
            "s4": s4, "c": c[lo:hi], "schwelle": float(beste["schwelle"]), "tag": tag, "met": t["met"],
            "datum": datum, "k": beste["k"]}


def detektoren(plt, a: dict, ziel: Path):
    reihen = np.argsort(a["winkel"])
    fig, achsen = plt.subplots(3, 4, figsize=(12.8, 7.2), sharex=True)
    oben = float(np.nanmax(a["zaehl"])) * 1.05
    for ax, d in zip(achsen.ravel(), reihen):
        hell = a["winkel"][d] < 60
        ax.plot(a["t_rel"], a["zaehl"][d], color=GOLD if hell else BLAU, linewidth=0.8, alpha=0.95)
        ax.plot(a["t_rel"], a["untergrund"][d], color=TEXT2, linewidth=1.0, linestyle="--", alpha=0.8)
        ax.set_ylim(0, oben)
        ax.set_title(f"Detektor {NAI[d]} · {de(a['winkel'][d])}° zur Quelle", fontsize=10.5,
                     color=GOLD if hell else TEXT)
        ax.tick_params(labelsize=8.5)
        komma(ax, "x")
    for ax in achsen[-1]:
        ax.set_xlabel("Sekunden seit dem Bord-Trigger", fontsize=9.5)
    for ax in achsen[:, 0]:
        ax.set_ylabel("Zählungen je 0,256 s", fontsize=9.5)
    fig.suptitle("Ein Gammablitz in allen zwölf Detektoren (22.10.2025, 50–300 keV)", x=0.06, ha="left",
                 fontsize=16, fontweight="bold", color=TEXT)
    fig.text(0.06, 0.905, "sortiert nach dem Winkel zur Quelle; gold: Detektoren, die zur Quelle zeigen (< 60°); "
             "gestrichelt: Untergrund", fontsize=11, color=TEXT2)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(ziel, dpi=125)
    plt.close(fig)


def treffer(plt, a: dict, ziel: Path, bis_s: float | None = None):
    """Drei Stufen: Summe der Rohzählungen, Signifikanz je Detektor, Koinzidenz mit Schwelle. bis_s: Animation."""
    from matplotlib.colors import LinearSegmentedColormap

    t = a["t_rel"]
    sicht = np.ones_like(t, bool) if bis_s is None else t <= bis_s
    reihen = np.argsort(a["winkel"])
    fig, (o, m, u) = plt.subplots(3, 1, figsize=(12.8, 7.2), sharex=True,
                                  gridspec_kw={"height_ratios": [1.1, 1.3, 1.1]})
    summe, ug = a["zaehl"].sum(axis=0), np.nansum(a["untergrund"], axis=0)
    o.plot(t[sicht], summe[sicht], color=BLAU, linewidth=0.8)
    o.plot(t[sicht], ug[sicht], color=TEXT2, linewidth=1.1, linestyle="--")
    o.set_ylim(np.nanmin(summe) * 0.95, np.nanmax(summe) * 1.05)
    o.set_title("1 · Rohdaten: Zählungen aller zwölf Detektoren je 0,256 s (gestrichelt: Untergrund)", fontsize=11.5)
    karte = LinearSegmentedColormap.from_list("treffer", [FLAECHE, "#1e3a8a", BLAU, GOLD, "#fff7d6"])
    s = np.where(sicht[None, :], np.clip(a["s4"][reihen], 0, 10), np.nan)
    bild = m.imshow(s, aspect="auto", cmap=karte, vmin=0, vmax=10, interpolation="nearest",
                    extent=(t[0], t[-1], 11.5, -0.5))
    m.set_yticks(range(12), [NAI[d] for d in reihen], fontsize=8.5)
    m.grid(False)
    m.set_title("2 · Signifikanz je Detektor über 1,024 s, sortiert nach dem Winkel zur Quelle "
                "(Farbskala bei 10 σ gekappt)", fontsize=11.5)
    leiste = fig.colorbar(bild, ax=m, pad=0.01, fraction=0.03)
    leiste.ax.tick_params(colors=TEXT2, labelsize=8.5)
    leiste.set_label("σ", color=TEXT2)
    leiste.outline.set_edgecolor(GITTER)
    c = np.where(sicht, a["c"], np.nan)
    u.plot(t, c, color=TEXT, linewidth=1.2)
    u.axhline(a["schwelle"], color=ORANGE, linewidth=1.3, linestyle=":")
    u.text(t[0] + 1, a["schwelle"] + 0.6, f"Schwelle {de(a['schwelle'], 2)} σ", color=ORANGE, fontsize=10)
    ueber = sicht & (a["c"] >= a["schwelle"])
    if ueber.any():
        u.fill_between(t, a["schwelle"], np.where(ueber, a["c"], np.nan), color=GOLD, alpha=0.35, linewidth=0)
        erst = t[np.argmax(ueber)]
        u.annotate("Treffer", (erst, a["schwelle"]), xytext=(erst + 3, np.nanmax(a["c"]) * 0.8),
                   color=GOLD, fontsize=13, fontweight="bold", arrowprops=dict(arrowstyle="-", color=GOLD))
    u.set_ylim(-2, np.nanmax(a["c"]) * 1.12)
    u.set_title(f"3 · Koinzidenz: {a['k']}.-stärkster Detektor über alle Zeitskalen, gegen die Schwelle", fontsize=11.5)
    u.set_xlabel("Sekunden seit dem Bord-Trigger")
    u.set_xlim(t[0], t[-1])
    for ax in (o, u):
        komma(ax)
    fig.suptitle("Vom Rohsignal zum Treffer (Gammablitz vom 22.10.2025)", x=0.06, ha="left", fontsize=16,
                 fontweight="bold", color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ziel, dpi=125)
    plt.close(fig)


def orbit(plt, d: dict, a: dict, ziel: Path):
    """24 Stunden: Zählrate je Detektor (16 s), Breite der Bahn, Südatlantische Anomalie, Katalog und Treffer."""
    tag, pos = a["tag"], a["tag"].pos
    n16 = 64
    n = tag.n // n16 * n16
    z = tag.zaehl["50-300"][:, :n].reshape(12, -1, n16).sum(axis=(0, 2)).astype(float)
    e = tag.dauer[:, :n].reshape(12, -1, n16).sum(axis=(0, 2))
    rate = np.where(e > 0.5 * 12 * n16 * DT * 0.5, z / np.where(e > 0, e, 1), np.nan)
    t0 = datetime.strptime(a["datum"], "%Y-%m-%d").replace(tzinfo=timezone.utc)
    met0 = tag.t0 + (t0 - utc_aus_met(tag.t0)).total_seconds()
    stunden = (tag.t0 + DT * n16 * (np.arange(rate.size) + 0.5) - met0) / 3600
    fig, (o, u) = plt.subplots(2, 1, figsize=(12.8, 7.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    o.plot(stunden, rate, color=BLAU, linewidth=1.0)
    ph = (pos["met"] - met0) / 3600
    saa = (pos["flags"] & 2) > 0
    kanten = np.flatnonzero(np.diff(saa.astype(int))) + 1
    for teil, h in zip(np.split(saa, kanten), np.split(ph, kanten)):
        if teil[0]:
            for ax in (o, u):
                ax.axvspan(h[0], h[-1], color=ORANGE, alpha=0.16, linewidth=0)
    tagestrig = [t for t in d["trig"].values() if utc_aus_met(t["met"]).strftime("%Y-%m-%d") == a["datum"]]
    oben = np.nanmax(rate)
    for t in tagestrig:
        h = (t["met"] - met0) / 3600
        farbe = GOLD if t["trigger_type"] == "GRB" else TEXT2
        o.axvline(h, color=farbe, linewidth=1.0, alpha=0.9)
        o.text(h + 0.08, oben * 1.02, {"GRB": "Gammablitz", "SFLARE": "Sonneneruption", "TGF": "TGF",
                                       "UNCERT": "unsicher", "LOCLPAR": "Teilchen", "DISTPAR": "Teilchen",
                                       "SGR": "Magnetar"}.get(t["trigger_type"], t["trigger_type"]),
               color=farbe, fontsize=9, rotation=90, va="bottom")
    eigene = [e for e in d["ereignisse"] if e["konfig"] == "beste" and e["utc"].startswith(a["datum"])]
    if eigene:
        h = [(float(e["met_spitze"]) - met0) / 3600 for e in eigene]
        o.scatter(h, np.full(len(h), oben * 1.17), marker="v", s=46, color=GOLD, zorder=5,
                  label="Treffer von gbmtransient")
        leg = o.legend(loc="upper left", fontsize=9.5)
        for t in leg.get_texts():
            t.set_color(TEXT)
    o.set_ylim(np.nanmin(rate) * 0.9, oben * 1.25)
    o.set_ylabel("Zählrate je Detektor (1/s)")
    o.set_title("Zählrate 50–300 keV, Mittel über die Detektoren, je 16 s; orange: Südatlantische Anomalie "
                "(Detektoren aus); Linien: Bord-Trigger", fontsize=11.5)
    u.plot(ph, pos["lat"], color=TEXT, linewidth=1.0)
    u.set_ylabel("Breite der Bahn (°)")
    u.set_xlabel(f"Stunden (UTC) am {a['datum'][8:10]}.{a['datum'][5:7]}.{a['datum'][:4]}")
    u.set_xlim(0, 24)
    u.set_xticks(range(0, 25, 3))
    for ax in (o, u):
        komma(ax)
    fig.suptitle("Ein Tag im Orbit: Der Untergrund folgt der Bahn", x=0.06, ha="left", fontsize=16,
                 fontweight="bold", color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ziel, dpi=125)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", type=Path, required=True)
    ap.add_argument("--nachweis", type=Path, required=True)
    ap.add_argument("--aus", type=Path, required=True)
    ap.add_argument("--bilder", type=Path, help="Ordner für Animationsbilder (für den Rundgang)")
    a = ap.parse_args()
    plt = stil()
    a.aus.mkdir(parents=True, exist_ok=True)
    d = lade(a.raw, a.nachweis)
    himmel(plt, d, a.aus / "himmel.png")
    if a.bilder:
        a.bilder.mkdir(parents=True, exist_ok=True)
        himmel(plt, d, a.bilder / "himmel_rein.png", groesse=(11.2, 6.3), titel=False)
    b = ausbruch(d, a.raw)
    detektoren(plt, b, a.aus / "detektoren.png")
    treffer(plt, b, a.aus / "treffer.png")
    orbit(plt, d, b, a.aus / "orbit.png")
    zahlen = {"grb": len(d["grb"]), "gefunden": sum(g["erkannt"] for g in d["grb"]),
              "beispiel_spitze": float(np.nanmax(b["c"])), "schwelle": b["schwelle"],
              "winkel_min": float(np.min(b["winkel"])), "hell_unter_60": int((b["winkel"] < 60).sum())}
    if a.bilder:
        h = a.bilder / "himmel"
        h.mkdir(parents=True, exist_ok=True)
        for i in range(len(d["grb"]) + 1):
            himmel(plt, d, h / f"h{i:03d}.png", anzahl=i)
        s = a.bilder / "treffer"
        s.mkdir(parents=True, exist_ok=True)
        for i, bis in enumerate(np.arange(FENSTER[0] + 2, FENSTER[1] + 0.1, 2.0)):
            treffer(plt, b, s / f"t{i:03d}.png", bis_s=float(bis))
        zahlen["bilder_himmel"] = len(d["grb"]) + 1
        zahlen["bilder_treffer"] = i + 1
    (a.aus / "schaubilder.json").write_text(json.dumps(zahlen, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(zahlen, ensure_ascii=False))


if __name__ == "__main__":
    main()
