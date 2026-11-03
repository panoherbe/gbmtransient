"""Künstliche Ausbrüche für die Gütemessung: Pulsform, Richtung, Detektorantwort, Poisson-Zählungen.

Pulsform nach Norris et al. (2005, ApJ 627, 324): I(t) = lambda * exp(-tau1/t - t/tau2) für t > 0, Spitze 1 bei
t = sqrt(tau1*tau2). Detektorantwort vereinfacht: Zählrate proportional zu max(cos theta_d, 0), theta_d = Winkel
zwischen Quellrichtung und Achse des NaI-Detektors d (Ausrichtungen nach Meegan et al. 2009, ApJ 702, 791, Tab. 1).
Streuung an Raumfahrzeug und Erdatmosphäre fehlen; im Bericht wird das Modell gegen echte Ausbrüche geprüft.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .einlesen import DT

# (Azimut, Zenitwinkel) in Grad im Koordinatensystem des Raumfahrzeugs, n0..nb
NAI_RICHTUNG = [(45.89, 20.58), (45.11, 45.31), (58.44, 90.21), (314.87, 45.24), (303.15, 90.27), (3.35, 89.79),
                (224.93, 20.43), (224.62, 46.18), (236.61, 89.97), (135.19, 45.55), (123.73, 90.42),
                (183.74, 90.32)]


def einheitsvektor(az_grad, zen_grad) -> np.ndarray:
    az, zen = np.radians(az_grad), np.radians(zen_grad)
    return np.stack([np.sin(zen) * np.cos(az), np.sin(zen) * np.sin(az), np.cos(zen)], axis=-1)


ACHSEN = einheitsvektor(np.array([a for a, _ in NAI_RICHTUNG]), np.array([z for _, z in NAI_RICHTUNG]))


def antwort(theta_grad, phi_grad) -> np.ndarray:
    """Relative Zählrate je Detektor (12,) für eine Quelle bei Zenitwinkel theta, Azimut phi."""
    q = einheitsvektor(phi_grad, theta_grad)
    return np.clip(ACHSEN @ q, 0.0, None)


@dataclass(frozen=True)
class Pulsklasse:
    name: str
    tau1: float  # s
    tau2: float  # s
    laenge: float  # s, ab Pulsbeginn simuliert
    verhaeltnis: float  # Überschuss 27-50 keV / 50-300 keV (Annahme, im Bericht gegen echte GRB geprüft)
    a_min: float  # cts/s, Spitzenrate 50-300 keV bei senkrechtem Einfall
    a_max: float


KLASSEN = {
    "kurz": Pulsklasse("kurz", 0.02, 0.25, 4.0, 0.3, 50.0, 5000.0),
    "lang": Pulsklasse("lang", 1.5, 6.0, 60.0, 0.6, 10.0, 1000.0),
}


def puls(t, tau1: float, tau2: float) -> np.ndarray:
    t = np.asarray(t, float)
    lam = np.exp(2 * np.sqrt(tau1 / tau2))
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        y = lam * np.exp(-tau1 / t - t / tau2)
    return np.where(t > 0, y, 0.0)


def zellmittel(k: Pulsklasse, versatz: float, n_zellen: int, unterteilung: int = 32) -> np.ndarray:
    """Mittelwert des Pulses je Rasterzelle; der Puls beginnt 'versatz' Sekunden nach dem Beginn der ersten Zelle."""
    t = (np.arange(n_zellen * unterteilung) + 0.5) * (DT / unterteilung) - versatz
    return puls(t, k.tau1, k.tau2).reshape(n_zellen, unterteilung).mean(axis=1)


def spitze_1024(k: Pulsklasse) -> float:
    """Höchstes Mittel des Pulses über 1,024 s (Spitze = 1): Umrechnung Spitzenrate -> Rate im 1,024-s-Fenster."""
    t = np.arange(0, k.laenge, 0.001) + 0.0005
    y = puls(t, k.tau1, k.tau2)
    c = np.concatenate([[0], np.cumsum(y) * 0.001])
    return float(np.max(c[1024:] - c[:-1024]) / 1.024)


def zufallsrichtung(rng: np.random.Generator) -> tuple[float, float]:
    """Isotrop verteilt: (Zenitwinkel, Azimut) in Grad."""
    return float(np.degrees(np.arccos(rng.uniform(-1, 1)))), float(rng.uniform(0, 360))


def einbauen(zaehl: dict[str, np.ndarray], dauer: np.ndarray, j0: int, versatz: float, k: Pulsklasse,
             amplitude: float, faktor: np.ndarray, rng: np.random.Generator) -> None:
    """Poisson-Zählungen eines Pulses in die Bänder 27-50 und 50-300 keV ab Zelle j0 addieren (in place)."""
    n_z = int(np.ceil((k.laenge + versatz) / DT))
    j1 = min(j0 + n_z, dauer.shape[1])
    mittel = zellmittel(k, versatz, j1 - j0)
    lebend = dauer[:, j0:j1]  # Lebenszeit je Zelle (s)
    erwartung = amplitude * faktor[:, None] * mittel[None, :] * lebend
    zaehl["50-300"][:, j0:j1] += rng.poisson(erwartung).astype(np.int32)
    zaehl["27-50"][:, j0:j1] += rng.poisson(erwartung * k.verhaeltnis).astype(np.int32)
