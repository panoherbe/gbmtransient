"""Zeitrechnung: Fermi-Missionszeit (MET) und UTC.

MET zählt SI-Sekunden seit 2001-01-01 00:00:00 UTC (FITS-Kopf: TIMESYS = TT, MJDREF = 51910 + 7,4287e-4). Seit
2001 wurden fünf Schaltsekunden eingefügt; die letzte am 31.12.2016. Die Kataloge des HEASARC führen die
Triggerzeit als MJD in UTC.
"""

from __future__ import annotations

from datetime import date, datetime, timedelta, timezone

import numpy as np

MET0 = datetime(2001, 1, 1, tzinfo=timezone.utc)
MJD_MET0 = 51910.0
# UTC-Tage, ab deren Beginn eine weitere Schaltsekunde gilt (eingefügt am Ende des Vortags)
SCHALTSEKUNDEN = [date(2006, 1, 1), date(2009, 1, 1), date(2012, 7, 1), date(2015, 7, 1), date(2017, 1, 1)]
_SCHALT_UTC_S = np.array([(datetime(d.year, d.month, d.day, tzinfo=timezone.utc) - MET0).total_seconds()
                          for d in SCHALTSEKUNDEN])
_SCHALT_MET = _SCHALT_UTC_S + np.arange(1, len(SCHALTSEKUNDEN) + 1)


def mjd_aus_datum(d: date) -> int:
    return (d - date(1858, 11, 17)).days


def met_aus_mjd(mjd_utc):
    """MJD (UTC) -> MET in Sekunden."""
    s = (np.asarray(mjd_utc, dtype=float) - MJD_MET0) * 86400.0
    return s + np.searchsorted(_SCHALT_UTC_S, s, side="right")


def met_aus_utc(t: datetime) -> float:
    s = (t - MET0).total_seconds()
    return float(s + np.searchsorted(_SCHALT_UTC_S, s, side="right"))


def utc_aus_met(met: float) -> datetime:
    n = int(np.searchsorted(_SCHALT_MET, met, side="right"))
    return MET0 + timedelta(seconds=float(met) - n)


def utc_text(met: float) -> str:
    return utc_aus_met(met).strftime("%Y-%m-%d %H:%M:%S")
