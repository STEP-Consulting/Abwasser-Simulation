# -*- coding: utf-8 -*-
"""
Rechenkern der Abwasser-Simulation (ohne Streamlit).

Konzepte:
  "ozon": Biologie + Ozonung
  "nf"  : Biologie + Nanofiltration (Bypass, Konzentrat-Rückführung, optional
          Produktwasser-GAK und Konzentrat-GAK mit Wiedervermischung)

Die Funktion `simuliere(p)` nimmt ein Parameter-Dictionary (siehe
`standard_parameter()`) und liefert Ergebnistabelle, Mengenströme und Kennzahlen.
Dadurch lassen sich Konzepte vergleichen, Sensitivitäten rechnen und Tests schreiben.
"""

import copy

import numpy as np
import pandas as pd
from scipy.stats import lognorm

# =========================================================
# KONSTANTEN
# =========================================================

ETA_PUMPE = 0.70
BETRIEBSSTUNDEN_JAHR = 8760
PH_SYS_FIX = 7.5
GAK_SCHUETTDICHTE_KG_M3 = 450.0          # Annahme
SZENARIO_ZIELWERT_PFOA_PFOS = 300.0      # ng/L, eigener Szenariowert (kein gesetzlicher Grenzwert)
ENERGIE_OZONUNG_AACHEN = 0.04824         # kWh/m³, nur Ozonung, KA Aachen-Sörs 2018
Z90 = 1.2816                             # 90 %-Quantil der Normalverteilung

GRENZWERTE_ABWV = {
    "Größenklasse 1": {"CSB": 150.0, "BSB5": 40.0},
    "Größenklasse 2": {"CSB": 110.0, "BSB5": 25.0},
    "Größenklasse 3": {"CSB": 90.0, "BSB5": 20.0},
    "Größenklasse 4": {"CSB": 90.0, "BSB5": 20.0},
    "Größenklasse 5": {"CSB": 75.0, "BSB5": 15.0},
}

VORKLAERUNG_ETA = {"0,75 - 1 h": 0.30, "1,5 - 2 h": 0.35, "> 2,5 h": 0.40}

# =========================================================
# STOFFDATENBANK
# =========================================================

_NAN = np.nan

BASE_DATA = {
    "Substanzname": [
        "CSB", "BSB5",
        "PFOA", "PFOS", "TFA", "PFHxA", "PFBS", "PFNA",
        "Diclofenac", "Ibuprofen", "Carbamazepin", "Paracetamol", "Naproxen",
    ],
    "Gruppe": [
        "Summenparameter", "Summenparameter",
        "PFAS", "PFAS", "PFAS", "PFAS", "PFAS", "PFAS",
        "PhAC", "PhAC", "PhAC", "PhAC", "PhAC",
    ],
    "Einheit": [
        "mg/L", "mg/L",
        "ng/L", "ng/L", "ng/L", "ng/L", "ng/L", "ng/L",
        "ng/L", "ng/L", "ng/L", "ng/L", "ng/L",
    ],
    # Summenparameter haben kein Molekulargewicht, keinen pKa, kein log K_ow (NaN = ohne Bedeutung)
    "MW": [
        _NAN, _NAN,
        414.0, 500.0, 114.0, 314.1, 300.1, 564.1,
        318.14, 206.28, 236.27, 151.17, 230.26,
    ],
    "pKa": [
        _NAN, _NAN,
        0.5, -3.27, 0.23, -0.16, -3.31, -0.21,
        4.08, 4.4, 13.9, 9.86, 4.15,
    ],
    "Ladung_z": [
        0, 0,
        -1, -1, -1, -1, -1, -1,
        -1, -1, 0, 0, -1,
    ],
    # Log_Kow wird im Modell derzeit nicht verwendet (nur Information)
    "Log_Kow": [
        _NAN, _NAN,
        4.67, 4.49, -0.2, 3.99, 1.05, 5.02,
        4.64, 3.68, 2.45, 0.49, 2.82,
    ],
    "Kd_m3_kg": [
        0.05, 0.02,
        0.329, 2.226, 0.0, 0.15, 0.007, 2.5,
        0.016, 0.007, 0.0, 0.0004, 0.013,
    ],
    # k_biol in L/(g SS*d) nach Joss et al. (2006); Paracetamol: van Bergen et al. (2021), umgerechnet
    # CSB nutzt f_S statt k_biol, BSB5 nutzt k als 1/d (siehe _biologie)
    "k_biol_L_gSS_d": [
        0.0, 3.5,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.1, 21.0, 0.01, 58, 1.0,
    ],
    "KF": [
        0.0, 0.0,
        0.29122326, 0.07315197, 0.5, 0.50608802, 1.15938066, 0.5,
        0.8865, 45.8, 24.4341655, 37.0566988, 12.1398782,
    ],
}

DEFAULT_CONCENTRATIONS = {
    "CSB": 450.0, "BSB5": 220.0,
    "PFOA": 300.0, "PFOS": 250.0, "TFA": 1000.0, "PFHxA": 400.0, "PFBS": 350.0, "PFNA": 150.0,
    "Diclofenac": 2000.0, "Ibuprofen": 5000.0, "Carbamazepin": 1500.0,
    "Paracetamol": 3000.0, "Naproxen": 1200.0,
}

OZON_PARAMETER = {
    "Diclofenac": {"k_o3": 5.26e5, "k_oh": 7.4e9},
    "Ibuprofen": {"k_o3": 10.0, "k_oh": 7.5e9},
    "Carbamazepin": {"k_o3": 3.00e5, "k_oh": 8.8e9},
    "Paracetamol": {"k_o3": 9.44e4, "k_oh": 5.2e9},
    "Naproxen": {"k_o3": 4.42e3, "k_oh": 9.3e9},
}


def get_stoffdaten():
    return pd.DataFrame(BASE_DATA)


def standard_parameter():
    """Alle Eingaben der Simulation mit Standardwerten."""
    return {
        "konzept": "nf",                # "ozon" oder "nf"
        "q_feed": 100.0,                # m³/h
        "bypass": 0.0,                  # Anteil des Bioablaufs ohne NF
        # Ozonung
        "ozon_dosis": 4.0,              # mg O3/L
        "opex_ozon": 0.018,             # €/m³ (pauschaler Richtwert)
        # NF
        "recovery": 0.80,
        "recycle": 0.20,
        "mwco": 200.0,                  # Da, Datenblattwert bei 90 % Rückhalt
        "zeta": -15.0,                  # mV
        "sigma": 0.3,
        "delta_p": 10.0,                # bar
        # GAK
        "prod_gak": False,
        "konz_gak": False,
        "prod_ebct": 15.0, "prod_tage": 180.0, "prod_doc": 0.35, "prod_dp": 0.6,
        "konz_ebct": 15.0, "konz_tage": 180.0, "konz_doc": 0.35, "konz_dp": 0.6,
        "r_doc": 0.90,                  # DOC-Rückhalt der NF (Aufkonzentrierung im Konzentrat)
        "mittelwert": True,             # mittlere Elimination über die Standzeit
        "lastabhaengig": False,         # BV50 abhängig von der Zulaufkonzentration (Freundlich)
        "n_inv": 0.5,                   # 1/n der Freundlich-Isotherme (Platzhalter)
        # Biologie
        "vkl_eta": 0.35,
        "hrt_h": 12.0, "srt": 15.0, "x_ts": 3.5, "f_s": 0.05,
        # Unsicherheitsfaktoren (für Sensitivität)
        "kd_faktor": 1.0, "kbiol_faktor": 1.0,
        "c_roh": copy.deepcopy(DEFAULT_CONCENTRATIONS),
    }


# =========================================================
# KOSTENANNAHMEN
# =========================================================

KOSTEN = {
    "strompreis": 0.20,                 # €/kWh

    # Ozonung
    "sauerstoffpreis": 0.15,            # €/kg O2
    "ozonkonzentration": 0.10,          # kg O3/kg O2
    "ozon_wartung": 0.005,              # €/m³ Zulauf
    "ozon_gesamt_opex": 0.018,          # €/m³ Zulauf

    # Nanofiltration
    "membranpreis": 40.0,               # €/m²
    "membranlebensdauer": 6.0,          # Jahre
    "membranfluss": 15.0,               # L/(m²*h)

    "antiscalant_dosis": 3.0,           # g/m³ NF-Zulauf
    "antiscalant_preis": 4.0,           # €/kg

    "cip_kosten": 0.005,                # €/m³ NF-Zulauf
    "nf_wartung": 0.010,                # €/m³ Permeat

    # Personal
    "personal_nf": 0.30,                # FTE
    "personal_gak": 0.10,               # FTE je GAK-Stufe
    "personal_ozon": 0.10,              # FTE
    "personalkosten": 75000.0,          # €/FTE*a

    # Aktivkohle
    "aktivkohlepreis": 1800.0,          # €/t
    "gak_betrieb": 0.003,               # €/m³ behandeltes Wasser

    # Konzentrat
    "konzentratentsorgung": 0.0,        # €/m³ Konzentrat
}


# =========================================================
# HYDRAULIK
# =========================================================

def mwco_zu_median(mwco, sigma):
    """MWCO (Datenblattwert bei 90 % Rückhalt) -> Median der lognormalen Porenverteilung."""
    return mwco / np.exp(Z90 * sigma)


def kreislauf_fluesse(q_feed, bypass, Y, r):
    Y = min(Y, 0.99)
    a = 1.0 - bypass
    q_bio = q_feed / (1.0 - r * (1.0 - Y) * a)
    q_nf = a * q_bio
    q_bypass = q_bio - q_nf
    q_perm = Y * q_nf
    q_konz_ges = q_nf - q_perm
    q_rez = r * q_konz_ges
    q_konz_aus = q_konz_ges - q_rez
    q_prod = q_perm + q_bypass
    return dict(a=a, q_bio=q_bio, q_nf=q_nf, q_bypass=q_bypass, q_perm=q_perm,
                q_konz_ges=q_konz_ges, q_rez=q_rez, q_konz_aus=q_konz_aus, q_prod=q_prod)


def fluesse_ozon(q_feed):
    return dict(a=1.0, q_bio=q_feed, q_nf=0.0, q_bypass=0.0, q_perm=0.0,
                q_konz_ges=0.0, q_rez=0.0, q_konz_aus=0.0, q_prod=q_feed)


# =========================================================
# VERFAHRENSSCHRITTE
# =========================================================

def _vorklaerung(c_roh, gruppe, eta):
    return c_roh * (1.0 - eta) if gruppe == "Summenparameter" else c_roh


def _biologie(c_in, kd, k_bio, name, gruppe, x_ts, srt, hrt_d, f_s):
    """Gibt (c_out, E_total, E_sorp, E_bio) als Anteile zurück."""
    ress = x_ts * hrt_d / srt
    e_sorp = (kd * ress) / (1.0 + kd * ress)

    if name == "CSB":
        e_sorp = 0.0
        e_bio = 1.0 - f_s
    elif gruppe == "Summenparameter":
        e_bio = (k_bio * srt) / (1.0 + k_bio * srt) if k_bio > 0 else 0.0
    else:
        kxh = k_bio * x_ts * hrt_d
        e_bio = kxh / (1.0 + kxh) if k_bio > 0 else 0.0

    e_total = float(np.clip(e_sorp + e_bio - e_sorp * e_bio, 0.0, 0.98))
    return c_in * (1.0 - e_total), e_total, e_sorp, e_bio


def _ozon(c_in, name, gruppe, dosis):
    """Gibt (c_out, Elimination als Anteil) zurück."""
    if gruppe != "PhAC":
        return c_in, 0.0
    par = OZON_PARAMETER.get(name)
    ct_o3 = (dosis / 4.0) * 1.5e-3
    ct_oh = (dosis / 4.0) * 1.0e-10
    if par:
        rem = 1.0 - np.exp(-(par["k_o3"] * ct_o3 + par["k_oh"] * ct_oh))
    else:
        rem = 0.85 * min(1.0, dosis / 4.0)
    rem = float(np.clip(rem, 0.0, 0.99))
    return c_in * (1.0 - rem), rem


def _nf_rueckhalt(row, median, sigma, zeta):
    """Gesamtrückhalt R und Größenanteil R_size."""
    if row["Gruppe"] == "Summenparameter":
        return 0.0, 0.0
    z = int(row["Ladung_z"])
    alpha = 1.0 / (1.0 + 10.0 ** (float(row["pKa"]) - PH_SYS_FIX)) if z < 0 else 0.0
    r_size = float(lognorm.cdf(float(row["MW"]), s=sigma, scale=median))
    charge = float(np.clip(abs(zeta / 30.0) * alpha * abs(z), 0.0, 0.95))
    r = float(np.clip(1.0 - (1.0 - r_size) * (1.0 - charge), 0.0, 0.99))
    return r, r_size


def _kreislauf_stoffstrom(c_vkl, E, R, Y, r, q_feed, fl, E_rez):
    Y = min(Y, 0.99)
    cpf = 1.0 - float(np.clip(R, 0.0, 1.0))
    q_ratio = q_feed / fl["q_bio"]
    nenner = 1.0 - (1.0 - E_rez) * r * fl["a"] * (1.0 - Y * cpf)
    c_bio = (1.0 - E) * c_vkl * q_ratio / nenner
    c_perm = c_bio * cpf
    c_konz = c_bio * (1.0 - Y * cpf) / (1.0 - Y)
    c_prod = (fl["q_perm"] * c_perm + fl["q_bypass"] * c_bio) / fl["q_prod"]
    return c_bio, c_perm, c_konz, c_prod


def gak_entfernung(c_in, row, ebct, tage, matrix_faktor, mittelwert=True,
                   lastabhaengig=False, n_inv=0.5, c_ref=0.0):
    """
    Elimination (Anteil 0..0,99) einer GAK-Stufe über Durchbruchskurve.

    BV50 gilt bei der Konzentration c_ref (Konzentration des NF-Zulaufs). Ist
    `lastabhaengig` aktiv, wird BV50 mit (c_in / c_ref)^(1/n - 1) skaliert
    (Freundlich: Kapazität q ~ C^(1/n), BV50 ~ q / C), begrenzt auf Faktor 0,1 bis 10.
    """
    if row["Gruppe"] == "Summenparameter" or (not np.isfinite(c_in)) or c_in <= 0.0:
        return 0.0
    kf = float(row["KF"])
    if kf <= 0.0 or ebct <= 0.0:
        return 0.0

    bed_volumes = 1440.0 * tage / ebct

    if row["Substanzname"] in ("TFA", "PFBS"):
        bv50_base = 4000.0
    elif row["Gruppe"] == "PFAS":
        bv50_base = 15000.0
    else:
        bv50_base = max(3000.0, kf * 600.0)

    bv50 = bv50_base * matrix_faktor
    if lastabhaengig and c_ref > 0.0:
        bv50 *= float(np.clip((c_in / c_ref) ** (n_inv - 1.0), 0.1, 10.0))

    steilheit = 3.5
    max_rem = 0.98 if row["Gruppe"] == "PFAS" else 0.99

    def rem_bei_bv(bv):
        durchbruch = 1.0 / (1.0 + np.exp(-steilheit * ((bv / bv50) - 1.0)))
        return max_rem * (1.0 - durchbruch)

    if mittelwert:
        rem = float(np.mean(rem_bei_bv(np.linspace(0.0, bed_volumes, 400))))
    else:
        rem = float(rem_bei_bv(bed_volumes))
    return float(np.clip(rem, 0.0, 0.99))


# =========================================================
# HAUPTFUNKTION
# =========================================================

def simuliere(p):
    """
    Führt die komplette Simulation aus.

    Rückgabe: dict mit
      df    - Ergebnistabelle je Substanz
      fl    - Mengenströme (m³/h)
      kenn  - Kennzahlen (Konzentrationen, Energie, Kohle, Bilanzen)
    """
    ist_ozon = (p["konzept"] == "ozon")
    q_feed = float(p["q_feed"])

    if ist_ozon:
        fl = fluesse_ozon(q_feed)
        prod_gak = konz_gak = False
    else:
        fl = kreislauf_fluesse(q_feed, p["bypass"], p["recovery"], p["recycle"])
        prod_gak = bool(p["prod_gak"])
        konz_gak = bool(p["konz_gak"])

    q_prod = fl["q_prod"]
    q_k = fl["q_konz_aus"]
    hrt_eff_h = p["hrt_h"] * q_feed / fl["q_bio"]
    median = mwco_zu_median(p["mwco"], p["sigma"])
    Y = min(p["recovery"], 0.99)

    # Matrix-Faktor der Konzentrat-GAK: DOC im Konzentrat ist aufkonzentriert
    aufkonz_doc = (1.0 - Y * (1.0 - p["r_doc"])) / (1.0 - Y)
    konz_doc_eff = max(0.02, p["konz_doc"] / aufkonz_doc)

    df = get_stoffdaten()
    rows = []

    for _, row in df.iterrows():
        name = row["Substanzname"]
        gruppe = row["Gruppe"]
        c_roh = float(p["c_roh"][name])
        c_vkl = _vorklaerung(c_roh, gruppe, p["vkl_eta"])

        kd = float(row["Kd_m3_kg"])
        k_bio = float(row["k_biol_L_gSS_d"])
        if gruppe != "Summenparameter":
            kd *= p["kd_faktor"]
            k_bio *= p["kbiol_faktor"]

        _, E, e_sorp, e_bio = _biologie(c_vkl, kd, k_bio, name, gruppe, p["x_ts"], p["srt"],
                                        hrt_eff_h / 24.0, p["f_s"])
        c_bio_ref, _, _, _ = _biologie(c_vkl, kd, k_bio, name, gruppe, p["x_ts"], p["srt"],
                                       p["hrt_h"] / 24.0, p["f_s"])

        fracht_in = q_feed * c_vkl
        r_nf = r_size = 0.0
        ozon_pct = prod_gak_pct = konz_gak_pct = 0.0
        akkum = 1.0
        f_bio = f_pg = f_kg = f_abzug = f_end = np.nan

        if ist_ozon:
            c_bio = c_vkl * (1.0 - E)
            c_end, ozon_rem = _ozon(c_bio, name, gruppe, p["ozon_dosis"])
            ozon_pct = ozon_rem * 100.0
            c_perm = c_prod_roh = c_prod_gak = c_end
            c_konz = c_konz_gak = np.nan
            bilanz = np.nan
        else:
            r_nf, r_size = _nf_rueckhalt(row, median, p["sigma"], p["zeta"])
            E_rez = 0.0 if name == "CSB" else E
            c_bio, c_perm, c_konz, c_prod_roh = _kreislauf_stoffstrom(
                c_vkl, E, r_nf, p["recovery"], p["recycle"], q_feed, fl, E_rez)
            akkum = (c_bio / c_bio_ref) if c_bio_ref > 0 else 1.0

            rem_p = 0.0
            if prod_gak:
                rem_p = gak_entfernung(c_prod_roh, row, p["prod_ebct"], p["prod_tage"], p["prod_doc"],
                                       p["mittelwert"], p["lastabhaengig"], p["n_inv"], c_bio)
            c_prod_gak = c_prod_roh * (1.0 - rem_p)
            prod_gak_pct = rem_p * 100.0

            if konz_gak:
                rem_k = gak_entfernung(c_konz, row, p["konz_ebct"], p["konz_tage"], konz_doc_eff,
                                       p["mittelwert"], p["lastabhaengig"], p["n_inv"], c_bio)
                c_konz_gak = c_konz * (1.0 - rem_k)
                konz_gak_pct = rem_k * 100.0
                c_end = (q_prod * c_prod_gak + q_k * c_konz_gak) / (q_prod + q_k)
                q_end = q_prod + q_k
            else:
                c_konz_gak = c_konz
                c_end = c_prod_gak
                q_end = q_prod

            if fracht_in > 0:
                f_bio = (E * q_feed * c_vkl + E_rez * fl["q_rez"] * c_konz) / fracht_in * 100.0
                f_pg = q_prod * (c_prod_roh - c_prod_gak) / fracht_in * 100.0
                if konz_gak:
                    f_kg = q_k * (c_konz - c_konz_gak) / fracht_in * 100.0
                    f_abzug = 0.0
                else:
                    f_kg = 0.0
                    f_abzug = q_k * c_konz / fracht_in * 100.0
                f_end = q_end * c_end / fracht_in * 100.0
                bilanz = f_bio + f_pg + f_kg + f_abzug + f_end
            else:
                bilanz = 100.0

        rows.append({
            "Substanzname": name, "Gruppe": gruppe, "Einheit": row["Einheit"],
            "C_Rohwasser": c_roh, "C_Vorklaerung": c_vkl,
            "Bio_Elimination_%": E * 100.0,
            "Bio_Sorptionsanteil_%": e_sorp * 100.0, "Bio_Abbauanteil_%": e_bio * 100.0,
            "C_Bio_Ablauf": c_bio, "Akkumulation_Bio_x": akkum,
            "R_size_%": r_size * 100.0, "R_NF_%": r_nf * 100.0,
            "Ozon_Elimination_%": ozon_pct,
            "C_Permeat": c_perm,
            "C_Produktwasser_roh": c_prod_roh,
            "GAK_Elimination_%": prod_gak_pct,
            "C_Produktwasser": c_prod_gak,
            "C_Konzentrat": c_konz,
            "Konz_GAK_Elimination_%": konz_gak_pct,
            "C_Konzentrat_nach_GAK": c_konz_gak,
            "C_Endwasser": c_end,
            "Fracht_Biologie_%": f_bio, "Fracht_Prod_GAK_%": f_pg, "Fracht_Konz_GAK_%": f_kg,
            "Fracht_Konz_Abzug_%": f_abzug, "Fracht_Endwasser_%": f_end,
            "Bilanzsumme_%": bilanz,
            "Konzentrationsreduktion_%": (1.0 - c_end / c_roh) * 100.0 if c_roh > 0 else np.nan,
            "Reduktion_Produktwasser_%": (1.0 - c_prod_gak / c_roh) * 100.0 if c_roh > 0 else np.nan,
        })

    df_res = pd.DataFrame(rows)

    # ---------------- Kennzahlen ----------------
    def conc(stoff, spalte="C_Endwasser"):
        return float(df_res.loc[df_res["Substanzname"] == stoff, spalte].iloc[0])

    phac = df_res["Gruppe"] == "PhAC"
    pfas = df_res["Gruppe"] == "PFAS"

    kenn = {
        "csb_end": conc("CSB"), "bsb_end": conc("BSB5"),
        "pfas_sum_end": float(df_res.loc[pfas, "C_Endwasser"].sum()),
        "pfas_sum_prod": float(df_res.loc[pfas, "C_Produktwasser"].sum()),
        "pfoa_pfos_end": conc("PFOA") + conc("PFOS"),
        "pfoa_pfos_prod": conc("PFOA", "C_Produktwasser") + conc("PFOS", "C_Produktwasser"),
        "avg_phac_red": float(df_res.loc[phac, "Konzentrationsreduktion_%"].mean()),
        "avg_phac_red_prod": float(df_res.loc[phac, "Reduktion_Produktwasser_%"].mean()),
        "hrt_eff_h": hrt_eff_h,
        "ress_eff_g_m3": p["x_ts"] * (hrt_eff_h / 24.0) / p["srt"] * 1000.0,
        "v_bb_m3": p["hrt_h"] * q_feed,
        "q_prod": q_prod, "q_perm": fl["q_perm"], "q_konz_aus": q_k, "q_rez": fl["q_rez"],
        "q_bio": fl["q_bio"], "q_nf": fl["q_nf"], "q_bypass": fl["q_bypass"],
        "wasserbilanz_abw": abs(q_prod + q_k - q_feed),
    }

    # Energie
    if ist_ozon:
        p_ges = ENERGIE_OZONUNG_AACHEN * q_feed
        kenn.update(sec_nf=np.nan, p_nf_kw=0.0, p_gak_prod_kw=0.0, p_gak_konz_kw=0.0, p_ges_kw=p_ges)
    else:
        sec_nf = p["delta_p"] / (36.0 * Y * ETA_PUMPE)
        p_nf = fl["q_perm"] * sec_nf
        p_gp = q_prod * p["prod_dp"] / (36.0 * ETA_PUMPE) if prod_gak else 0.0
        p_gk = q_k * p["konz_dp"] / (36.0 * ETA_PUMPE) if konz_gak else 0.0
        p_ges = p_nf + p_gp + p_gk
        kenn.update(sec_nf=sec_nf, p_nf_kw=p_nf, p_gak_prod_kw=p_gp, p_gak_konz_kw=p_gk, p_ges_kw=p_ges)
    kenn["e_jahr_kwh"] = p_ges * BETRIEBSSTUNDEN_JAHR
    kenn["spez_energie_kwh_m3"] = p_ges / q_feed          # je m³ behandeltem Zulauf

    # Aktivkohle
    def kohle(aktiv, ebct, tage, q):
        if not aktiv or ebct <= 0 or tage <= 0:
            return 0.0, 0.0
        bv_wechsel = 1440.0 * tage / ebct
        g_m3 = GAK_SCHUETTDICHTE_KG_M3 * 1000.0 / bv_wechsel
        return g_m3, g_m3 * q * BETRIEBSSTUNDEN_JAHR / 1.0e6

    kenn["kohle_prod_g_m3"], kenn["kohle_prod_t_a"] = kohle(prod_gak, p["prod_ebct"], p["prod_tage"], q_prod)
    kenn["kohle_konz_g_m3"], kenn["kohle_konz_t_a"] = kohle(konz_gak, p["konz_ebct"], p["konz_tage"], q_k)
    kenn["kohle_gesamt_t_a"] = kenn["kohle_prod_t_a"] + kenn["kohle_konz_t_a"]
    kenn["konz_doc_eff"] = konz_doc_eff
    kenn["aufkonz_doc"] = aufkonz_doc

   
    # =====================================================
    # KOSTENBERECHNUNG
    # Bezugsgröße: €/m³ behandelter Zulauf
    # =====================================================

    jahresmenge = q_feed * BETRIEBSSTUNDEN_JAHR

    if ist_ozon:

        # Pauschaler Gesamt-OPEX der Ozonung
        # Sauerstoff, Strom, Wartung und Personal werden
        # nicht zusätzlich addiert (Vermeidung Doppelzählung).

        kosten_ozon = float(p["opex_ozon"])

       kenn["opex_ozon_eur_a"] = (
    float(p["opex_ozon"]) * jahresmenge
    if ist_ozon else np.nan
)

        kenn["kosten_nf_gesamt_m3"] = np.nan
        kenn["kosten_nf_gesamt_eur_a"] = np.nan

    else:

        # ---------------------------------------------
        # 1. Stromkosten
        # ---------------------------------------------

        kosten_strom = (
            kenn["spez_energie_kwh_m3"]
            * KOSTEN["strompreis"]
        )

        # ---------------------------------------------
        # 2. Membrankosten (annualisierte CAPEX)
        # ---------------------------------------------

        membranflaeche = (
            fl["q_perm"] * 1000.0
            / KOSTEN["membranfluss"]
        )

        membrankosten_jahr = (
            membranflaeche
            * KOSTEN["membranpreis"]
            / KOSTEN["membranlebensdauer"]
        )

        kosten_membran = membrankosten_jahr / jahresmenge

        # ---------------------------------------------
        # 3. Antiscalant
        # ---------------------------------------------

        kosten_antiscalant = (
            KOSTEN["antiscalant_dosis"] / 1000.0
            * KOSTEN["antiscalant_preis"]
            * fl["q_nf"] / q_feed
        )

        # ---------------------------------------------
        # 4. Reinigungschemikalien (CIP)
        # ---------------------------------------------

        kosten_cip = (
            KOSTEN["cip_kosten"]
            * fl["q_nf"] / q_feed
        )

        # ---------------------------------------------
        # 5. Wartung NF
        # ---------------------------------------------

        kosten_wartung = (
            KOSTEN["nf_wartung"]
            * fl["q_perm"] / q_feed
        )

        # ---------------------------------------------
        # 6. Personalkosten NF
        # ---------------------------------------------

        kosten_personal_nf = (
            KOSTEN["personal_nf"]
            * KOSTEN["personalkosten"]
            / jahresmenge
        )

        # ---------------------------------------------
        # 7. Aktivkohle
        # ---------------------------------------------

        kosten_aktivkohle = (
            kenn["kohle_gesamt_t_a"]
            * KOSTEN["aktivkohlepreis"]
            / jahresmenge
        )

        # ---------------------------------------------
        # 8. GAK-Betrieb
        # ---------------------------------------------

        kosten_gak_betrieb = (
            (fl["q_prod"] if prod_gak else 0.0)
            + (q_k if konz_gak else 0.0)
        ) / q_feed * KOSTEN["gak_betrieb"]

        # ---------------------------------------------
        # 9. Personal GAK
        # ---------------------------------------------

        anzahl_gak = int(prod_gak) + int(konz_gak)

        kosten_personal_gak = (
            anzahl_gak
            * KOSTEN["personal_gak"]
            * KOSTEN["personalkosten"]
            / jahresmenge
        )

        # ---------------------------------------------
        # 10. Konzentratentsorgung
        # ---------------------------------------------

        kosten_konzentrat = (
            q_k / q_feed
            * KOSTEN["konzentratentsorgung"]
        )

        # ---------------------------------------------
        # 11. Gesamtkosten NF
        # ---------------------------------------------

        kosten_nf = (
            kosten_strom
            + kosten_membran
            + kosten_antiscalant
            + kosten_cip
            + kosten_wartung
            + kosten_personal_nf
            + kosten_aktivkohle
            + kosten_gak_betrieb
            + kosten_personal_gak
            + kosten_konzentrat
        )

        # Einzelkosten speichern

        kenn.update({
            "kosten_strom": kosten_strom,
            "kosten_membran": kosten_membran,
            "kosten_antiscalant": kosten_antiscalant,
            "kosten_cip": kosten_cip,
            "kosten_wartung": kosten_wartung,
            "kosten_personal_nf": kosten_personal_nf,
            "kosten_aktivkohle": kosten_aktivkohle,
            "kosten_gak_betrieb": kosten_gak_betrieb,
            "kosten_personal_gak": kosten_personal_gak,
            "kosten_konzentrat": kosten_konzentrat,

            "kosten_nf_gesamt_m3": kosten_nf,
            "kosten_nf_gesamt_eur_a": kosten_nf * jahresmenge,

            "membranflaeche_m2": membranflaeche,
            "membrankosten_eur_a": membrankosten_jahr,
        })

        kenn["kosten_ozon_gesamt_m3"] = np.nan
        kenn["kosten_ozon_gesamt_eur_a"] = np.nan

    # Bestehende Ozon-OPEX-Kennzahl beibehalten
    kenn["opex_ozon_eur_a"] = (
        KOSTEN["ozon_gesamt_opex"] * jahresmenge
        if ist_ozon else np.nan
    )

# ---------------------------------------------------------
    # FALLBACK / SICHERSTELLUNG ALLER KENNZAHLEN (Ozon vs. NF)
    # ---------------------------------------------------------
    if "kosten_nf_gesamt_m3" not in kenn:
        kenn["kosten_nf_gesamt_m3"] = np.nan
        kenn["kosten_nf_gesamt_eur_a"] = np.nan
    if "kosten_ozon_gesamt_m3" not in kenn:
        kenn["kosten_ozon_gesamt_m3"] = np.nan
        kenn["kosten_ozon_gesamt_eur_a"] = np.nan
        
    return {"df": df_res, "fl": fl, "kenn": kenn}

# =========================================================
# SENSITIVITÄT (one-at-a-time)
# =========================================================

ZIELGROESSEN = {
    "Ø Arzneimittel-Reduktion Endwasser (%)": lambda k: k["avg_phac_red"],
    "PFAS-Summe Endwasser (ng/L)": lambda k: k["pfas_sum_end"],
    "PFAS-Summe Produktwasser (ng/L)": lambda k: k["pfas_sum_prod"],
    "Spez. Energiebedarf (kWh/m³ Zulauf)": lambda k: k["spez_energie_kwh_m3"],
}

# (Label, Schlüssel, niedrig, hoch, Modus "rel" = Faktor, "abs" = Absolutwert, nur_wenn)
SENS_PARAMETER = [
    ("MWCO (±20 %)", ["mwco"], 0.8, 1.2, "rel", None),
    ("Porenverteilung σ (0,2 / 0,5)", ["sigma"], 0.2, 0.5, "abs", None),
    ("Zeta-Potential (−25 / −5 mV)", ["zeta"], -25.0, -5.0, "abs", None),
    ("Ausbeute Y (0,7 / 0,9)", ["recovery"], 0.7, 0.9, "abs", None),
    ("Membrandruck ΔP (±30 %)", ["delta_p"], 0.7, 1.3, "rel", None),
    ("Bypass-Anteil (0 / 0,2)", ["bypass"], 0.0, 0.2, "abs", None),
    ("Konzentrat-Rückführung (0 / 0,4)", ["recycle"], 0.0, 0.4, "abs", None),
    ("K_d (×0,5 / ×2)", ["kd_faktor"], 0.5, 2.0, "abs", None),
    ("k_biol (×0,5 / ×2)", ["kbiol_faktor"], 0.5, 2.0, "abs", None),
    ("Matrix-Faktor GAK (±40 %)", ["prod_doc", "konz_doc"], 0.6, 1.4, "rel", "gak"),
    ("GAK-Standzeit (±50 %)", ["prod_tage", "konz_tage"], 0.5, 1.5, "rel", "gak"),
]


def sensitivitaet(p_basis, zielgroesse):
    """One-at-a-time-Sensitivität für das NF-Konzept. Gibt DataFrame zurück."""
    ziel = ZIELGROESSEN[zielgroesse]
    p0 = copy.deepcopy(p_basis)
    p0["konzept"] = "nf"
    basis = ziel(simuliere(p0)["kenn"])
    gak_aktiv = bool(p0["prod_gak"] or p0["konz_gak"])

    zeilen = []
    for label, keys, lo, hi, modus, nur_wenn in SENS_PARAMETER:
        if nur_wenn == "gak" and not gak_aktiv:
            continue
        ergebnisse = []
        for wert in (lo, hi):
            pv = copy.deepcopy(p0)
            for k in keys:
                pv[k] = pv[k] * wert if modus == "rel" else wert
            ergebnisse.append(ziel(simuliere(pv)["kenn"]))
        zeilen.append({"Parameter": label, "Ergebnis_niedrig": ergebnisse[0],
                       "Ergebnis_hoch": ergebnisse[1], "Basis": basis,
                       "Spanne": abs(ergebnisse[1] - ergebnisse[0])})
    out = pd.DataFrame(zeilen)
    return out.sort_values("Spanne", ascending=True).reset_index(drop=True)
