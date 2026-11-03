"""Tests ohne Netz: Zeitrechnung, Statistik, Raster, Simulation, Optimierung und ein Ende-zu-Ende-Lauf auf
synthetischen Dateien im Format der GBM-Daten und der HEASARC-Kataloge."""

import hashlib
import io
import json
import math
from datetime import date, datetime, timezone

import numpy as np
import pytest

from gbmtransient import einlesen, optimierung, quellen, simulation, statistik, suche, zeit


# --------------------------------------------------------------------------------------------- Zeitrechnung

def test_met_grb170817a():
    # GBM-Trigger zu GW170817: 17.08.2017 12:41:06 UTC = MET 524666471 (fünf Schaltsekunden seit 2001)
    t = datetime(2017, 8, 17, 12, 41, 6, tzinfo=timezone.utc)
    assert zeit.met_aus_utc(t) == 524666471.0
    assert zeit.utc_aus_met(524666471.0) == t


def test_met_tagesbeginn_2025():
    # erste Zelle der CTIME-Datei vom 01.06.2025 beginnt 0,5 s vor Mitternacht (MET 770428804,49)
    assert zeit.met_aus_utc(datetime(2025, 6, 1, tzinfo=timezone.utc)) == 770428805.0
    mjd = zeit.mjd_aus_datum(date(2025, 6, 1))
    assert float(zeit.met_aus_mjd(mjd)) == 770428805.0


# -------------------------------------------------------------------------------------------------- Statistik

def _lr_direkt(n_on, n_off, alpha):
    """-2 ln Lambda für Poisson an/aus, direkt aus den Likelihoods (Fakultäten kürzen sich)."""
    def ll(n, mu):
        return n * math.log(mu) - mu if n > 0 else -mu
    b = (n_on + n_off) / (1 + alpha)
    l0 = ll(n_on, alpha * b) + ll(n_off, b)
    l1 = ll(n_on, max(n_on, 1e-300)) + ll(n_off, max(n_off, 1e-300))
    return 2 * (l1 - l0)


@pytest.mark.parametrize("n_on,n_off,alpha", [(120, 400, 0.25), (80, 400, 0.25), (15, 1000, 0.01),
                                              (0, 50, 0.2), (300, 2600, 0.125)])
def test_li_ma_gleich_likelihood_quotient(n_on, n_off, alpha):
    s = float(statistik.li_ma(n_on, n_off, alpha))
    assert s * s == pytest.approx(_lr_direkt(n_on, n_off, alpha), rel=1e-9, abs=1e-9)
    assert np.sign(s) == np.sign(n_on - alpha * n_off)


def test_li_ma_nullverteilung():
    rng = np.random.default_rng(1)
    alpha = 0.125
    b = 400.0
    s = statistik.li_ma(rng.poisson(alpha * b, 40000), rng.poisson(b, 40000), alpha)
    assert abs(s.mean()) < 0.03
    assert abs(s.std() - 1) < 0.03


def _reihe(n=6000, rate=300.0, steigung=0.0, seed=2):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * einlesen.DT
    mu = (rate + steigung * t) * einlesen.DT
    z = rng.poisson(np.tile(mu, (12, 1))).astype(np.int32)
    d = np.full((12, n), einlesen.DT, np.float32)
    return z, d


def _s(z, d, m, w, g, modus):
    rand = m + g + w + 1
    return statistik.signifikanz(statistik.kumulieren(z, rand), statistik.kumulieren(d, rand), z.shape[1], rand,
                                 m, w, g, modus, float(np.median(d)))


def test_signifikanz_ruhig_und_luecke():
    z, d = _reihe()
    d[:, 3000:3400] = 0  # Datenlücke
    z[:, 3000:3400] = 0
    s = _s(z, d, 4, 64, 8, "zweiseitig")
    assert s.shape == (12, z.shape[1] + 1)
    assert np.all(np.isnan(s[:, 3050:3350]))
    ok = np.isfinite(s)
    assert abs(np.nanmean(s)) < 0.05 and abs(np.nanstd(s) - 1) < 0.1
    assert ok.mean() > 0.8


def test_linearer_untergrund_zweiseitig_ohne_verzerrung_kausal_mit():
    z, d = _reihe(n=20000, rate=200.0, steigung=0.2)  # 0,2 cts/s pro s: Anstieg wie im Orbit
    m, w, g = 32, 128, 8
    zw = np.nanmean(_s(z, d, m, w, g, "zweiseitig"))
    ka = np.nanmean(_s(z, d, m, w, g, "kausal"))
    assert abs(zw) < 0.05
    assert ka > 0.3  # kausal hinkt dem Anstieg hinterher und meldet Überschuss


def test_ordnungsstatistik():
    s = np.array([[1.0, np.nan], [5.0, 2.0], [3.0, np.nan]])
    c = statistik.ordnungsstatistik(s, (1, 2, 3))
    assert list(c[1]) == [5.0, 2.0]
    assert list(c[2]) == [3.0, -np.inf]
    assert c[3][0] == 1.0


def test_ereignisse_verbinden_und_spitze():
    c = np.zeros(100, np.float32)
    c[10:13] = [5, 7, 6]
    c[15] = 5.5  # Lücke 2 -> gleiches Ereignis bei verbinden=3
    c[60] = 9
    a, e, p, sp = statistik.ereignisse(c, 5.0, 3)
    assert list(a) == [10, 60] and list(e) == [15, 60] and list(p) == [11, 60]
    assert list(sp) == [7, 9]


def test_intervalle():
    lo, hi = statistik.wilson(0, 10)
    assert lo == pytest.approx(0, abs=1e-12) and hi == pytest.approx(0.2775, abs=1e-3)
    assert statistik.poisson_intervall(0)[1] == pytest.approx(3.689, abs=1e-3)
    lo, hi = statistik.poisson_intervall(10)
    assert lo == pytest.approx(4.795, abs=1e-3) and hi == pytest.approx(18.390, abs=1e-3)


def test_logistik_findet_parameter():
    rng = np.random.default_rng(3)
    x = rng.uniform(1, 4, 20000)
    p = 1 / (1 + np.exp(-(-10 + 4 * x)))
    y = rng.random(x.size) < p
    a, b, cov = statistik.logistik_anpassen(x, y)
    assert a == pytest.approx(-10, abs=0.5) and b == pytest.approx(4, abs=0.2)
    x50, se = statistik.logistik_quantil(a, b, cov, 0.5)
    assert x50 == pytest.approx(2.5, abs=0.02) and 0 < se < 0.02


# ------------------------------------------------------------------------------------------- Raster, Modell

def test_auf_raster_erhaelt_zaehlungen_und_messzeit():
    t0 = 1000.0
    start = list(t0 + 0.256 * np.arange(10))  # zehn Standard-Intervalle
    fein = t0 + 2.56 + 0.064 * np.arange(16)  # danach 1,024 s in 0,064-s-Intervallen
    start += list(fein)
    start = np.array(start)
    ende = start + np.where(np.arange(start.size) < 10, 0.256, 0.064)
    start = np.append(start, ende[-1] + 5.0)  # nach einer Lücke von 5 s
    ende = np.append(ende, start[-1] + 0.256)
    zaehl = np.ones((start.size, 8), np.int32) * np.arange(1, 9)
    exposure = (ende - start) * 0.99
    quality = np.zeros(start.size, int)
    quality[3] = 1  # verworfen
    z, e = einlesen.auf_raster(start, ende, exposure, quality, zaehl, t0, 40)
    assert z.sum(axis=0)[0] == start.size - 1
    assert e.sum() == pytest.approx(exposure.sum() - exposure[3])
    assert z[3].sum() == 0 and e[3] == 0  # QUALITY != 0
    assert z[10, 0] == 4 and e[10] == pytest.approx(4 * 0.064 * 0.99)  # vier feine Intervalle in einer Zelle
    assert e[14:33].sum() == 0  # Lücke


def test_antwort_und_puls():
    assert np.allclose(np.linalg.norm(simulation.ACHSEN, axis=1), 1)
    az, zen = simulation.NAI_RICHTUNG[0]
    f = simulation.antwort(zen, az)
    assert np.argmax(f) == 0 and f[0] == pytest.approx(1)
    k = simulation.KLASSEN["lang"]
    t_spitze = math.sqrt(k.tau1 * k.tau2)
    assert float(simulation.puls(t_spitze, k.tau1, k.tau2)) == pytest.approx(1)
    assert 0.9 < simulation.spitze_1024(k) <= 1
    assert simulation.spitze_1024(simulation.KLASSEN["kurz"]) < 0.5


def test_einbau_erwartungswert():
    z = {"27-50": np.zeros((12, 400), np.int32), "50-300": np.zeros((12, 400), np.int32)}
    d = np.full((12, 400), 0.256, np.float32)
    k = simulation.KLASSEN["lang"]
    f = np.zeros(12)
    f[2] = 1.0
    rng = np.random.default_rng(4)
    for _ in range(50):
        simulation.einbauen(z, d, 10, 0.1, k, 100.0, f, rng)
    erwartet = 50 * 100.0 * simulation.zellmittel(k, 0.1, 235).sum() * 0.256
    assert z["50-300"][2].sum() == pytest.approx(erwartet, rel=0.03)
    assert z["27-50"][2].sum() == pytest.approx(erwartet * k.verhaeltnis, rel=0.04)
    assert z["50-300"][0].sum() == 0


def test_verschiebungen_trennen_jedes_paar():
    n = suche.PHI.size
    assert n == suche.N_VERSCHIEBUNGEN and np.unique(suche.PHI).size == n
    for p in suche.PHI:
        for q in range(1, 12):
            assert abs(q * p - round(q * p)) >= suche.MIN_ABSTAND - 1e-12


def test_konfigurationen_vollstaendig():
    assert len(suche.KONFIGS) == 162 and len(set(suche.KONFIGS)) == 162


def test_schwelle_aus_histogramm():
    nb = suche.KANTEN.size - 1
    agg = {"hist_verschoben": np.zeros((2, nb)), "gueltig_verschoben_s": np.full(2, 10 * 86400.0),
           "hist_unbestaetigt": np.zeros((2, nb)), "gueltig_s": np.full(2, 86400.0)}
    agg["hist_verschoben"][0, 50] = 5  # 5 Ereignisse bei 4,5 in 10 Tagen -> 0,5/Tag
    agg["hist_unbestaetigt"][0, 100] = 3  # 3 bei 5,5 an einem Tag
    agg["hist_unbestaetigt"][1, -1] = 10  # nie zulässig
    s, bindend = optimierung.schwellen(agg, far_max=0.1, ru_max=2.0)
    assert s[0] == pytest.approx(suche.KANTEN[101]) and bindend[0] == "pruefaufwand"
    assert not np.isfinite(s[1])


def test_verzeichnis_hoechste_version():
    html = ('<a href="glg_ctime_n0_250901_v00.pha">x</a><a href="glg_ctime_n0_250901_v01.pha">x</a>'
            '<a href="glg_ctime_nb_250901_v00.pha">x</a><a href="glg_poshist_all_250901_v00.fit">x</a>'
            '<a href="glg_ctime_n1_250902_v00.pha">x</a>')
    d = quellen.dateien_im_verzeichnis(html, date(2025, 9, 1))
    assert d == {"n0": "glg_ctime_n0_250901_v01.pha", "nb": "glg_ctime_nb_250901_v00.pha",
                 "poshist": "glg_poshist_all_250901_v00.fit"}


# ------------------------------------------------------------------------------------------- Ende zu Ende

def _ctime_fits(t_start, n, rate_kanal, extra=None, seed=0):
    from astropy.io import fits

    rng = np.random.default_rng(seed)
    start = t_start + 0.256 * np.arange(n)
    mu = np.tile(np.array(rate_kanal) * 0.256, (n, 1))
    if extra is not None:
        mu = mu + extra
    counts = rng.poisson(mu).astype(np.int16)
    cols = [fits.Column("COUNTS", "8I", array=counts), fits.Column("EXPOSURE", "E", array=np.full(n, 0.2547)),
            fits.Column("QUALITY", "I", array=np.zeros(n, np.int16)), fits.Column("TIME", "D", array=start),
            fits.Column("ENDTIME", "D", array=start + 0.256)]
    sp = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    e = [4.5, 11.8, 26.9, 50.4, 101.6, 293.8, 537.8, 983.3, 2000.0]
    eb = fits.BinTableHDU.from_columns([fits.Column("CHANNEL", "I", array=np.arange(8)),
                                        fits.Column("E_MIN", "E", array=e[:-1]),
                                        fits.Column("E_MAX", "E", array=e[1:])], name="EBOUNDS")
    b = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), eb, sp]).writeto(b)
    return b.getvalue()


def _poshist(t_start, n_s):
    from astropy.io import fits

    t = t_start + np.arange(n_s, dtype=float)
    cols = [fits.Column("SCLK_UTC", "D", array=t)] + [fits.Column(f"QSJ_{i}", "D", array=np.zeros(n_s) + (i == 4))
                                                      for i in range(1, 5)]
    cols += [fits.Column("SC_LAT", "E", array=np.zeros(n_s)), fits.Column("SC_LON", "E", array=np.zeros(n_s)),
             fits.Column("FLAGS", "I", array=np.ones(n_s, np.int16))]
    b = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols, name="GLAST POS HIST")]).writeto(b)
    return b.getvalue()


def _votable(spalten: dict):
    from astropy.io.votable import from_table
    from astropy.table import Table

    b = io.BytesIO()
    from_table(Table(spalten)).to_xml(b)
    return b.getvalue()


def test_ende_zu_ende(tmp_path):
    """Zwei kurze synthetische Tage (je 1 h), am ersten ein Ausbruch in vier Detektoren mit Katalog-Trigger."""
    raw = tmp_path / "raw" / "2026-01-01"
    man = []

    def ablegen(rel, b):
        p = raw / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b)
        man.append({"pfad": rel, "url": "https://example.invalid/" + rel, "bytes": len(b),
                    "sha256": hashlib.sha256(b).hexdigest(), "abgerufen": "2026-01-01T00:00:00Z"})

    rate = [80, 340, 200, 150, 130, 30, 20, 60]
    n = 14063  # 1 h
    tage = [date(2025, 9, 1), date(2025, 9, 2)]
    t_burst = zeit.met_aus_utc(datetime(2025, 9, 1, 0, 30, tzinfo=timezone.utc))
    for i, d in enumerate(tage):
        t0 = zeit.met_aus_utc(datetime(d.year, d.month, d.day, tzinfo=timezone.utc)) - 0.5
        for k, det in enumerate(einlesen.NAI):
            extra = None
            if i == 0 and k in (0, 1, 3, 5):
                tt = t0 + 0.256 * np.arange(n) + 0.128 - t_burst
                puls = simulation.puls(tt, 0.5, 3.0) * 1500 * 0.256
                extra = np.outer(puls, [0, 0.2, 0.5, 0.6, 0.4, 0.1, 0, 0])
            rel = f"heasarc/FTP/fermi/data/gbm/daily/{d:%Y/%m/%d}/current/glg_ctime_{det}_{d:%y%m%d}_v00.pha"
            ablegen(rel, _ctime_fits(t0, n, rate, extra, seed=10 * i + k))
        ablegen(f"heasarc/FTP/fermi/data/gbm/daily/{d:%Y/%m/%d}/current/glg_poshist_all_{d:%y%m%d}_v00.fit",
                _poshist(t0, 3600))
    mjd = zeit.mjd_aus_datum(tage[0]) + (t_burst - zeit.met_aus_utc(datetime(2025, 9, 1, tzinfo=timezone.utc))) / 86400
    ablegen("heasarc/tap/fermigtrig_2025-09-01_2025-09-02.vot",
            _votable({"trigger_name": ["bn250901021"], "trigger_time": [mjd], "trigger_type": ["GRB"],
                      "theta": [30.0], "phi": [45.0]}))
    ablegen("heasarc/tap/fermigbrst_2025-09-01_2025-09-02.vot",
            _votable({"trigger_name": ["bn250901021"], "trigger_time": [mjd], "t90": [12.0], "flux_1024": [5.0]}))
    (raw / "manifest.jsonl").write_text("".join(json.dumps(e) + "\n" for e in man), encoding="utf-8")

    from gbmtransient.bericht import auswerten

    aus = tmp_path / "nachweis"
    try:
        import matplotlib  # noqa: F401
        mit_abb = True
    except ImportError:
        mit_abb = False
    r = auswerten(raw, aus, prozesse=1, abbildungen=mit_abb)
    assert (aus / "ergebnisse.json").exists() and (aus / "BERICHT.md").exists()
    assert r["daten"]["tage"] == 2 and r["katalog"]["trigger_im_zeitraum"] == 1
    grb = r["bewertung"]["alle"]["beste"]["katalog"]["GRB"]
    assert grb["abgedeckt"] == 1 and grb["erkannt"] == 1
    assert r["nullverteilung"]["std"] == pytest.approx(1, abs=0.1)
    if mit_abb:
        assert len(list((aus / "abbildungen").glob("*.png"))) == 6
    zeilen = (aus / "katalog_abgleich.csv").read_text(encoding="utf-8").splitlines()
    assert zeilen[1].startswith("bn250901021,GRB")


def test_pruefsumme_wird_geprueft(tmp_path):
    p = tmp_path / "a.bin"
    p.write_bytes(b"abc")
    e = {"pfad": "a.bin", "sha256": hashlib.sha256(b"abd").hexdigest()}
    with pytest.raises(einlesen.PruefsummenFehler):
        einlesen.lies_geprueft(tmp_path, e)


def test_holen_ohne_netz(tmp_path, monkeypatch):
    listing = '<a href="glg_ctime_n0_250901_v00.pha">'.encode() + b"".join(
        f'<a href="glg_ctime_{d}_250901_v00.pha">'.encode() for d in quellen.NAI[1:]) + \
        b'<a href="glg_poshist_all_250901_v00.fit">'
    aufrufe = []

    def hole(url, timeout=0, versuche=0):
        aufrufe.append(url)
        return listing if url.endswith("/") else url.encode()

    monkeypatch.setattr(quellen, "hole", hole)
    r = quellen.holen(tmp_path, date(2025, 9, 1), date(2025, 9, 1), "test", arbeiter=2)
    assert r["fehler"] == [] and r["dateien"] == 2 + 1 + 13
    n = len(aufrufe)
    r = quellen.holen(tmp_path, date(2025, 9, 1), date(2025, 9, 1), "test", arbeiter=2)
    assert len(aufrufe) == n  # zweiter Lauf: nichts neu geholt
