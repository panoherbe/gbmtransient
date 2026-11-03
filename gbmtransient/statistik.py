"""Statistische Verfahren: Signifikanz je Detektor, Koinzidenz über Detektoren, Ereignisbildung, Intervalle.

Signifikanz: Likelihood-Quotienten-Test für Poisson-Zählungen mit Untergrund aus Nebenfenstern nach Li & Ma
(1983, ApJ 272, 317, Gl. 17); alpha = Messzeit im Signalfenster / Messzeit in den Untergrundfenstern. Unter der
Nullhypothese (kein Signal) ist S asymptotisch standardnormalverteilt.

Untergrund: zweiseitig (Fenster der Breite w vor und nach dem Signalfenster, Abstand g; ein linearer Verlauf des
Untergrunds hebt sich heraus) oder kausal (nur das Fenster davor; für Echtzeit-Alarme).

Koinzidenz: die k-größte Signifikanz über die zwölf Detektoren (Ordnungsstatistik). Ein Trigger mit "mindestens
k Detektoren über s" ist genau C_k >= s.
"""

from __future__ import annotations

import math

import numpy as np

from .einlesen import DT


def li_ma(n_on, n_off, alpha):
    """Signifikanz nach Li & Ma (1983), Gl. 17, mit Vorzeichen des Überschusses."""
    n_on = np.asarray(n_on, dtype=float)
    n_off = np.asarray(n_off, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    tot = n_on + n_off
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(n_on > 0, n_on * np.log((1 + alpha) / alpha * n_on / tot), 0.0)
        b = np.where(n_off > 0, n_off * np.log((1 + alpha) * n_off / tot), 0.0)
    s = np.sqrt(np.maximum(2.0 * (a + b), 0.0))
    return np.sign(n_on - alpha * n_off) * s


def kumulieren(x: np.ndarray, rand: int) -> np.ndarray:
    """Kumulierte Summe über die letzte Achse mit 'rand' Nullzellen auf jeder Seite (für Fenstersummen)."""
    p = np.zeros(x.shape[:-1] + (x.shape[-1] + 2 * rand + 1,), dtype=np.float64)
    p[..., rand + 1:rand + 1 + x.shape[-1]] = x
    return np.cumsum(p, axis=-1)


def _fenster(c: np.ndarray, m: int, start: int, anzahl: int) -> np.ndarray:
    """Summen über m Zellen, die an den Zellgrenzen start-rand .. start-rand+anzahl-1 enden."""
    return c[..., start:start + anzahl] - c[..., start - m:start - m + anzahl]


def signifikanz(cz: np.ndarray, ce: np.ndarray, n: int, rand: int, m: int, w: int, g: int, modus: str,
                voll: float, min_on: float = 0.75, min_off: float = 0.5) -> np.ndarray:
    """S(j) für ein Signalfenster von m Zellen, das an Zellgrenze j endet; j = 0..n, NaN wo ungültig.

    cz, ce: kumulierte Zählungen und Lebenszeiten (kumulieren(x, rand), rand >= m + g + w + 1). w, g in Zellen.
    Gültig, wenn das Signalfenster mindestens min_on und jedes Untergrundfenster mindestens min_off seiner
    Sollmesszeit hat (voll = Sollmesszeit einer Zelle).
    """
    if rand < m + g + w + 1:
        raise ValueError("rand zu klein")
    k = n + 1
    n_on, e_on = _fenster(cz, m, rand, k), _fenster(ce, m, rand, k)
    n_off, e_links = _fenster(cz, w, rand - m - g, k), _fenster(ce, w, rand - m - g, k)
    e_off = e_links.copy()
    ok = (e_on >= min_on * m * voll) & (e_links >= min_off * w * voll)
    if modus == "zweiseitig":
        e_rechts = _fenster(ce, w, rand + g + w, k)
        n_off = n_off + _fenster(cz, w, rand + g + w, k)
        e_off += e_rechts
        ok &= e_rechts >= min_off * w * voll
    elif modus != "kausal":
        raise ValueError(modus)
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = np.where(ok, e_on / np.where(e_off > 0, e_off, 1.0), 1.0)
        s = li_ma(n_on, n_off, alpha)
    s[~ok] = np.nan
    return s.astype(np.float32)


def ueberschuss(cz: np.ndarray, ce: np.ndarray, n: int, rand: int, m: int, w: int, g: int,
                voll: float) -> tuple[np.ndarray, np.ndarray]:
    """Überschuss N_on - alpha*N_off (zweiseitiger Untergrund) und Messzeit E_on je Fensterende j = 0..n."""
    k = n + 1
    n_on, e_on = _fenster(cz, m, rand, k), _fenster(ce, m, rand, k)
    n_off = _fenster(cz, w, rand - m - g, k) + _fenster(cz, w, rand + g + w, k)
    e_l, e_r = _fenster(ce, w, rand - m - g, k), _fenster(ce, w, rand + g + w, k)
    ok = (e_on >= 0.75 * m * voll) & (e_l >= 0.5 * w * voll) & (e_r >= 0.5 * w * voll)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(ok, n_on - e_on / (e_l + e_r) * n_off, np.nan)
    return x, np.where(ok, e_on, np.nan)


def ordnungsstatistik(s: np.ndarray, ks=(1, 2, 3)) -> dict[int, np.ndarray]:
    """k-größter Wert über die Detektoren (Achse 0); NaN zählt als -inf."""
    x = np.where(np.isnan(s), -np.inf, s)
    x = np.sort(x, axis=0)
    return {k: x[-k].copy() for k in ks}


def ereignisse(c: np.ndarray, schwelle: float, verbinden: int):
    """Zusammenhängende Bereiche mit c >= schwelle (Lücken bis 'verbinden' Indizes werden überbrückt).

    Liefert Arrays (anfang, ende, spitze_index, spitzenwert); Indizes in c.
    """
    idx = np.flatnonzero(c >= schwelle)
    if idx.size == 0:
        leer = np.zeros(0, dtype=np.int64)
        return leer, leer, leer, np.zeros(0, dtype=np.float32)
    neu = np.flatnonzero(np.diff(idx) > verbinden) + 1
    starts = np.concatenate([[0], neu])
    enden = np.concatenate([neu, [idx.size]]) - 1
    werte = c[idx]
    spitzen = np.maximum.reduceat(werte, starts)
    # Index der Spitze: erster Index je Gruppe mit Wert == Spitze
    gruppe = np.repeat(np.arange(starts.size), np.diff(np.concatenate([starts, [idx.size]])))
    ist_spitze = werte == spitzen[gruppe]
    erste = np.full(starts.size, -1, dtype=np.int64)
    pos = np.flatnonzero(ist_spitze)
    g = gruppe[pos]
    _, u = np.unique(g, return_index=True)
    erste[g[u]] = idx[pos[u]]
    return idx[starts], idx[enden], erste, spitzen.astype(np.float32)


# ---------------------------------------------------------------------------------------------- Intervalle

def wilson(k, n, z: float = 1.959964):
    """Wilson-Intervall (95 %) für einen Anteil k/n."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, k / n, np.nan)
        nenner = 1 + z * z / n
        mitte = (p + z * z / (2 * n)) / nenner
        halb = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / nenner
    return mitte - halb, mitte + halb


def _poisson_cdf(k: int, mu: float) -> float:
    if mu <= 0:
        return 1.0
    log_term = -mu
    s = math.exp(log_term)
    for i in range(1, k + 1):
        log_term += math.log(mu) - math.log(i)
        s += math.exp(log_term)
    return min(s, 1.0)


def poisson_intervall(k: int, cl: float = 0.95) -> tuple[float, float]:
    """Exaktes (Garwood-)Intervall für den Erwartungswert einer Poisson-Zählung k."""
    a = (1 - cl) / 2

    def suche(f, lo, hi):
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid):
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    hi_max = k + 20 + 10 * math.sqrt(k + 1)
    unten = 0.0 if k == 0 else suche(lambda mu: 1 - _poisson_cdf(k - 1, mu) >= a, 0.0, hi_max)
    oben = suche(lambda mu: _poisson_cdf(k, mu) <= a, 0.0, hi_max)
    return unten, oben


def logistik_anpassen(x, y, iterationen: int = 100, ridge: float = 1e-6):
    """Maximum-Likelihood-Anpassung P(y=1) = 1/(1+exp(-(a+b x))) per Newton-Raphson.

    Liefert (a, b, Kovarianzmatrix aus der Fisher-Information); NaN, wenn alle y gleich sind (nicht bestimmbar).
    Ein winziger Ridge-Term hält die Iteration bei vollständiger Trennung endlich.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if y.size == 0 or y.min() == y.max():
        return float("nan"), float("nan"), np.full((2, 2), np.nan)
    X = np.stack([np.ones_like(x), x], axis=1)
    beta = np.zeros(2)
    for _ in range(iterationen):
        p = 1 / (1 + np.exp(-np.clip(X @ beta, -500, 500)))
        H = X.T @ (X * (p * (1 - p))[:, None]) + ridge * np.eye(2)
        schritt = np.linalg.solve(H, X.T @ (y - p) - ridge * beta)
        beta = beta + schritt
        if np.max(np.abs(schritt)) < 1e-10:
            break
    p = 1 / (1 + np.exp(-np.clip(X @ beta, -500, 500)))
    H = X.T @ (X * (p * (1 - p))[:, None]) + ridge * np.eye(2)
    return float(beta[0]), float(beta[1]), np.linalg.inv(H)


def logistik_quantil(a: float, b: float, cov: np.ndarray, q: float) -> tuple[float, float]:
    """x mit P = q und Standardfehler (Delta-Methode); NaN, wenn die Kurve nicht bestimmt ist."""
    if not (math.isfinite(a) and math.isfinite(b)) or b == 0:
        return float("nan"), float("nan")
    L = math.log(q / (1 - q))
    x = (L - a) / b
    grad = np.array([-1 / b, -(L - a) / b ** 2])
    return x, float(math.sqrt(max(grad @ cov @ grad, 0.0)))
