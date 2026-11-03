"""Abbildungen des Berichts (Matplotlib, 300 dpi). Jede Abbildung liest nur Ergebnisse der Auswertung."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import optimierung, simulation, suche
from .statistik import logistik_anpassen, wilson

BLAU, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"  # Kategorien 1-3 (geprüfte Reihenfolge)
TINTE, TINTE2, GITTER = "#0b0b0b", "#52514e", "#e4e3df"


def _stil():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": TINTE2, "axes.labelcolor": TINTE,
        "xtick.color": TINTE2, "ytick.color": TINTE2, "text.color": TINTE, "axes.grid": True,
        "grid.color": GITTER, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.dpi": 300, "figure.dpi": 100, "lines.linewidth": 1.4,
        "axes.formatter.use_locale": False,
    })
    return plt


def _komma(ax, achse="both"):
    from matplotlib.ticker import FuncFormatter

    f = FuncFormatter(lambda v, _: (f"{v:g}").replace(".", ","))
    if achse in ("x", "both"):
        ax.xaxis.set_major_formatter(f)
    if achse in ("y", "both"):
        ax.yaxis.set_major_formatter(f)


def alle(aus: Path, r, tage, det, w, j_test, beste, grund, train, test) -> None:
    plt = _stil()
    ziel = aus / "abbildungen"
    ziel.mkdir(parents=True, exist_ok=True)
    abb1(plt, ziel, det)
    abb2(plt, ziel, tage, w, beste, train)
    abb3(plt, ziel, r, w, j_test, beste, grund)
    abb4(plt, ziel, r, tage, beste, grund, test, w)
    abb5(plt, ziel, r)
    abb6(plt, ziel, det)


def abb1(plt, ziel, det):
    b = next((d["beispiel"] for d in det if d["beispiel"]), None)
    fig, (o, u) = plt.subplots(2, 1, figsize=(6.3, 3.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    if b:
        t = np.array(b["zeit_rel"])
        for i, (z, ug, d) in enumerate(zip(b["zaehl"], b["untergrund"], b["detektoren"])):
            farbe = (BLAU, ORANGE, AQUA)[i]
            o.step(t, z, where="post", color=farbe, lw=0.9, label=f"Detektor {d}")
            o.plot(t + 0.128, ug, color=farbe, lw=1.4, ls="--", alpha=0.9)
        o.set_ylabel("Zählungen je 0,256 s")
        o.legend(loc="upper right", ncol=3)
        o.set_title(f"{b['name']}: Zählungen 50–300 keV (Linie gestrichelt: Untergrund, zweiseitig)", loc="left")
        c = np.array([np.nan if v is None else v for v in b["koinzidenz"]], float)
        u.plot(t, c, color=TINTE, lw=1.0)
        u.axhline(b["schwelle"], color=ORANGE, lw=1.2, ls=":")
        u.text(t[0], b["schwelle"], f" Schwelle {b['schwelle']:.2f} σ".replace(".", ","), va="bottom", fontsize=7,
               color=TINTE2)
        u.set_ylabel("Koinzidenz (σ)")
        u.set_xlabel("Zeit seit Bord-Trigger (s)")
        _komma(u, "x")
    fig.tight_layout()
    fig.savefig(ziel / "abb1_beispiel.png")
    plt.close(fig)


def abb2(plt, ziel, tage, w, beste, train):
    fig, (l, r_) = plt.subplots(1, 2, figsize=(6.3, 2.7))
    sh = sum(t["s_hist"].sum(axis=0) for t in tage).astype(float)
    k = suche.S_KANTEN
    mitte = 0.5 * (k[1:] + k[:-1])
    breite = k[1] - k[0]
    dichte = sh / (sh.sum() * breite)
    l.semilogy(mitte, np.where(dichte > 0, dichte, np.nan), drawstyle="steps-mid", color=BLAU, label="gemessen")
    x = np.linspace(-6, 6, 400)
    l.semilogy(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color=TINTE2, ls="--", lw=1.0, label="Normalverteilung")
    l.set_xlim(-6, 14)
    l.set_ylim(1e-10, 1)
    l.set_xlabel("Signifikanz S je Detektor")
    l.set_ylabel("Dichte")
    l.legend(loc="upper right")
    l.set_title("Nullverteilung (1,024 s, 50–300 keV)", loc="left")
    agg = optimierung.summieren(tage, train)
    kanten = suche.KANTEN[:-1]
    for key, g, farbe, name in (("hist_unbestaetigt", "gueltig_s", ORANGE, "unbestätigt"),
                                ("hist_bestaetigt", "gueltig_s", BLAU, "im Katalog"),
                                ("hist_verschoben", "gueltig_verschoben_s", AQUA, "Zufall (verschoben)")):
        rate = optimierung.rate_ab(agg[key][beste], agg[g][beste])
        r_.semilogy(kanten, np.where(rate > 0, rate, np.nan), color=farbe, label=name)
    r_.axvline(w["s"][beste], color=TINTE2, ls=":", lw=1.0)
    r_.axhline(optimierung.RU_MAX, color=ORANGE, ls=":", lw=0.8)
    r_.axhline(optimierung.FAR_MAX, color=AQUA, ls=":", lw=0.8)
    r_.set_xlim(suche.BODEN, 15)
    r_.set_xlabel("Schwelle (σ)")
    r_.set_ylabel("Ereignisse je Tag")
    r_.legend(loc="upper right")
    r_.set_title("Raten der gewählten Konfiguration", loc="left")
    _komma(l, "x")
    _komma(r_, "x")
    fig.tight_layout()
    fig.savefig(ziel / "abb2_raten.png")
    plt.close(fig)


def abb3(plt, ziel, r, w, j_test, beste, grund):
    fig, (l, r_) = plt.subplots(1, 2, figsize=(6.3, 2.8))
    farben = {1: BLAU, 2: ORANGE, 3: AQUA}
    for k in (1, 2, 3):
        for modus, marke in (("zweiseitig", "o"), ("kausal", "^")):
            idx = [i for i, c in enumerate(suche.KONFIGS) if c.k == k and c.modus == modus and np.isfinite(w["s"][i])]
            l.scatter(100 * w["j"][idx], 100 * j_test[idx], s=16, marker=marke, color=farben[k], alpha=0.8,
                      edgecolors="white", linewidths=0.4, label=f"k={k}, {modus}")
    l.scatter([100 * w["j"][beste]], [100 * j_test[beste]], s=90, marker="*", color=TINTE, label="gewählt", zorder=5)
    l.scatter([100 * w["j"][grund]], [100 * j_test[grund]], s=40, marker="s", facecolors="none", edgecolors=TINTE,
              label="Grundeinstellung", zorder=5)
    lim = [0, 100]
    l.plot(lim, lim, color=GITTER, lw=0.8, zorder=0)
    l.set_xlabel("J Trainingstage (%)")
    l.set_ylabel("J Testtage (%)")
    l.legend(loc="upper left", fontsize=6.3, ncol=1, handletextpad=0.2)
    l.set_title("Alle Konfigurationen", loc="left")
    p = r["optimierung"]["pareto"]
    bud = [x["budget_je_tag"] for x in p]
    r_.semilogx(bud, [100 * x["j_test"] for x in p], marker="o", color=BLAU, label="gewählt je Budget")
    r_.semilogx(bud, [100 * x["grund_j_test"] for x in p], marker="s", color=ORANGE, label="Grundeinstellung")
    r_.axvline(optimierung.RU_MAX, color=TINTE2, ls=":", lw=0.8)
    r_.set_xlabel("Prüfbudget: unbestätigte Kandidaten je Tag")
    r_.set_ylabel("J Testtage (%)")
    r_.legend(loc="lower right")
    r_.set_title("Empfindlichkeit nach Prüfbudget", loc="left")
    _komma(r_, "x")
    fig.tight_layout()
    fig.savefig(ziel / "abb3_optimierung.png")
    plt.close(fig)


def abb4(plt, ziel, r, tage, beste, grund, test, w):
    inj = optimierung.injektionen(tage, test)
    fig, achsen = plt.subplots(1, 2, figsize=(6.3, 2.8), sharey=True)
    rug = r.get("grb_spitzenraten", {})
    for ax, kl, name in ((achsen[0], 0, "kurz"), (achsen[1], 1, "lang")):
        sel = (inj["klasse"] == kl) & (inj["amplitude"] > 0)
        x = inj["staerke"][sel]
        kanten = np.logspace(np.log10(x.min()), np.log10(x.max()), 16)
        mitte = np.sqrt(kanten[1:] * kanten[:-1])
        for i_cfg, farbe, label, versatz in ((beste, BLAU, "gewählt", 0.985), (grund, ORANGE, "Grundeinstellung",
                                                                                1.015)):
            treffer = inj["max"][i_cfg, sel] >= w["s"][i_cfg]
            b = np.digitize(x, kanten) - 1
            n = np.array([np.sum(b == i) for i in range(mitte.size)])
            k = np.array([np.sum(treffer[b == i]) for i in range(mitte.size)])
            lo, hi = wilson(k, n)
            v = n > 0
            q = k[v] / n[v]
            ax.errorbar(mitte[v] * versatz, 100 * q, yerr=[100 * np.clip(q - lo[v], 0, None),
                                                          100 * np.clip(hi[v] - q, 0, None)],
                        fmt="o", ms=3.2, color=farbe, elinewidth=0.8, capsize=0, label=label)
            a, bb, _ = logistik_anpassen(np.log10(x), treffer)
            xs = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
            ax.plot(xs, 100 / (1 + np.exp(-(a + bb * np.log10(xs)))), color=farbe, lw=1.2)
        werte = rug.get(name, [])
        if werte:
            ax.plot(werte, np.full(len(werte), -4), "|", color=TINTE2, ms=7, mew=0.8, label="echte GRB (gemessen)")
        ax.set_xscale("log")
        ax.set_ylim(-8, 104)
        ax.set_xlabel("Spitzenrate 1,024 s im hellsten Detektor (cts/s)")
        ax.set_title(f"Puls „{name}“", loc="left")
        _komma(ax, "x")
    achsen[0].set_ylabel("Nachweiswahrscheinlichkeit (%)")
    achsen[1].legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ziel / "abb4_effizienz.png")
    plt.close(fig)


def abb5(plt, ziel, r):
    kb = r["bewertung"]["alle"]["beste"]["katalog"]
    kg = r["bewertung"]["alle"]["grund"]["katalog"]
    typen = sorted(kb, key=lambda t: -kb[t]["abgedeckt"])
    fig, ax = plt.subplots(figsize=(6.3, 2.5))
    x = np.arange(len(typen))
    for versatz, kat, farbe, label in ((-0.2, kb, BLAU, "gewählt"), (0.2, kg, ORANGE, "Grundeinstellung")):
        q = np.array([100 * kat[t]["quote"] for t in typen])
        lo = np.array([100 * kat[t]["ki95"][0] for t in typen])
        hi = np.array([100 * kat[t]["ki95"][1] for t in typen])
        ax.bar(x + versatz, q, width=0.38, color=farbe, label=label, edgecolor="white", linewidth=1.0)
        ax.errorbar(x + versatz, q, yerr=[np.clip(q - lo, 0, None), np.clip(hi - q, 0, None)], fmt="none",
                    ecolor=TINTE2, elinewidth=0.8)
    ax.set_xticks(x, [f"{t}\n(n = {kb[t]['abgedeckt']})" for t in typen])
    ax.set_ylabel("erkannt (%)")
    ax.set_ylim(0, 105)
    ax.legend(loc="upper right", ncol=2)
    ax.set_title("Erkannte Bord-Trigger je Einordnung im Katalog (alle Tage)", loc="left")
    fig.tight_layout()
    fig.savefig(ziel / "abb5_katalog.png")
    plt.close(fig)


def abb6(plt, ziel, det):
    tr = [x for d in det for x in d["trigger"] if x["typ"] == "GRB" and "modell" in x and "muster" in x]
    fig, (l, r_) = plt.subplots(1, 2, figsize=(6.3, 2.6), gridspec_kw={"width_ratios": [1.4, 1]})
    xs, ys, kor = [], [], []
    for x in tr:
        m = np.array(x["muster"], float)
        if m.max() <= 0:
            continue
        xs.append(np.array(x["modell"]))
        ys.append(m / m.max())
        if x.get("korrelation") is not None:
            kor.append(x["korrelation"])
    if xs:
        l.scatter(np.concatenate(xs), np.concatenate(ys), s=8, color=BLAU, alpha=0.45, edgecolors="none")
    l.set_xlabel("Modell max(cos θ, 0)")
    l.set_ylabel("Überschuss / hellster Detektor")
    l.set_title("Muster über die zwölf Detektoren", loc="left")
    if kor:
        r_.hist(kor, bins=np.linspace(-1, 1, 21), color=BLAU, edgecolor="white", linewidth=1.0)
    r_.set_xlabel("Korrelation Messung–Modell")
    r_.set_ylabel("Gammablitze")
    r_.set_title("je Ausbruch", loc="left")
    _komma(l)
    _komma(r_, "x")
    fig.tight_layout()
    fig.savefig(ziel / "abb6_antwort.png")
    plt.close(fig)
