# -*- coding: utf-8 -*-
"""
Integrierte Abwasser-Simulation: Verfahrensvergleich & Pilotanlagen-Vorbereitung
1. Bestandssystem: Biologie + Ozonung
2. Advanced System: Nanofiltration mit Konzentrat-GAK & Produkt-GAK 
   - Exaktes Schwarz-Weiß-Design
   - Saubere, übersichtliche Leitungsführung ohne Linien-Überschneidungen
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import altair as alt
from scipy.stats import lognorm
from html import escape

# =========================================================
# 1. STREAMLIT-KONFIGURATION & PASSWORDSCHUTZ
# =========================================================

st.set_page_config(
    page_title="Abwasser-Simulation: Verfahrensvergleich",
    page_icon="💧",
    layout="wide"
)

PASSWORD = "Passwort"  # Hier dein Wunschpasswort eintragen

def check_password():
    """Gibt True zurück, wenn das Passwort korrekt eingegeben wurde."""
    if "password_correct" not in st.session_state:
        st.session_state["password_correct"] = False

    if st.session_state["password_correct"]:
        return True

    st.subheader("🔒 Passwortgeschützter Bereich")
    pwd = st.text_input("Bitte gib das Passwort ein, um die Simulation zu starten:", type="password")
    
    if st.button("Anmelden"):
        if pwd == PASSWORD:
            st.session_state["password_correct"] = True
            st.rerun()
        else:
            st.error("Falsches Passwort!")
    return False

if not check_password():
    st.stop()

st.title("💧 Abwasser-Simulation: Verfahrens- & Konzeptvergleich")

with st.expander("💡 Bedienung & Modellübersicht", expanded=False):
    st.markdown("""
    ## So funktioniert das Konzept (Konzentrat-GAK & Vermischung)

    ### 1. Konzentrat-Management
    - Ein Teil des NF-Konzentrats wird in die Biologie zurückgeführt.
    - Der Reststrom des Konzentrats wird über eine **Konzentrat-GAK** gereinigt und hinter der Haupt-/Polishing-Stufe mit dem gereinigten Endwasser **vermischt**.
    - Das Fließbild zeigt nun die exakte Einmündung hinter den jeweiligen GAK-Stufen ohne Überlappungen in klarem Schwarz-Weiß-Design.
    """)

# =========================================================
# 2. STOFFDATENBANK & PARAMETER
# =========================================================

base_data = {
    "Substanzname": [
        "CSB", "BSB5",
        "PFOA", "PFOS", "TFA", "PFHxA", "PFBS", "PFNA",
        "Diclofenac", "Ibuprofen", "Carbamazepin",
        "Paracetamol", "Naproxen"
    ],
    "Gruppe": [
        "Summenparameter", "Summenparameter",
        "PFAS", "PFAS", "PFAS", "PFAS", "PFAS", "PFAS",
        "PhAC", "PhAC", "PhAC", "PhAC", "PhAC"
    ],
    "Einheit": [
        "mg/L", "mg/L",
        "ng/L", "ng/L", "ng/L", "ng/L", "ng/L", "ng/L",
        "ng/L", "ng/L", "ng/L", "ng/L", "ng/L"
    ],
    "MW": [
        0, 0,
        414.0, 500.0, 114.0, 314.1, 300.1, 564.1,
        318.14, 206.28, 236.27, 151.17, 230.26
    ],
    "pKa": [
        0.0, 0.0,
        0.5, -3.27, 0.23, -0.16, -3.31, -0.21,
        4.08, 4.4, 13.9, 9.86, 4.15
    ],
    "Ladung_z": [
        0, 0,
        -1, -1, -1, -1, -1, -1,
        -1, -1, 0, 0, -1
    ],
    "Log_Kow": [
        0.0, 0.0,
        4.67, 4.49, -0.2, 3.99, 1.05, 5.02,
        4.64, 3.68, 2.45, 0.49, 2.82
    ],
    "Kd_m3_kg": [
        0.05, 0.02,
        0.329, 2.226, 0.0, 0.15, 0.007, 2.5,
        0.016, 0.007, 0.0, 0.0, 0.013
    ],
    "k_biol_L_gSS_d": [
        0, 3.5,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.1, 21.0, 0.01, 3.5, 1.0
    ],
    "KF": [
        0.0, 0.0,
        0.29122326, 0.07315197, 0.5, 0.50608802, 1.15938066, 0.5,
        0.8865, 45.8, 24.4341655, 37.0566988, 12.1398782
    ]
}

df_stoff = pd.DataFrame(base_data)

default_concentrations = {
    "CSB": 450.0, "BSB5": 220.0,
    "PFOA": 300.0, "PFOS": 250.0, "TFA": 1000.0, "PFHxA": 400.0, "PFBS": 350.0, "PFNA": 150.0,
    "Diclofenac": 2000.0, "Ibuprofen": 5000.0, "Carbamazepin": 1500.0, "Paracetamol": 3000.0, "Naproxen": 1200.0
}

OZON_PARAMETER = {
    "Diclofenac": {"k_o3": 5.26e5, "k_oh": 7.4e9},
    "Ibuprofen": {"k_o3": 10.0, "k_oh": 7.5e9},
    "Carbamazepin": {"k_o3": 3.00e5, "k_oh": 8.8e9},
    "Paracetamol": {"k_o3": 9.44e4, "k_oh": 5.2e9},
    "Naproxen": {"k_o3": 4.42e3, "k_oh": 9.3e9},
}

# =========================================================
# 3. SEITENLEISTE
# =========================================================

st.sidebar.header("⚙️ Verfahrensauswahl")

verfahrens_wahl = st.sidebar.radio(
    "Konzept",
    ["Bestandssystem (Ozonung)", "Nanofiltration (NF)"]
)

is_ozon = (verfahrens_wahl == "Bestandssystem (Ozonung)")

with st.sidebar.expander("🔀 Hydraulik & Mengenströme", expanded=True):
    q_feed = st.slider("Rohwasserzulauf Q (m³/h)", 10.0, 500.0, 100.0, 5.0)
    if not is_ozon:
        bypass_fraction = st.slider("Bypass-Anteil Bioablauf", 0.0, 0.80, 0.0, 0.05)
    else:
        bypass_fraction = 0.0

has_bypass = (not is_ozon) and (bypass_fraction > 0.0)

if is_ozon:
    with st.sidebar.expander("🧪 Ozonung", expanded=True):
        ozon_dosis = st.slider("Ozon-Dosis (mg O₃ / L)", 1.0, 10.0, 4.0, 0.5)
    aktivieren_gak = False
    aktivieren_konz_gak = False
    gak_ebct, gak_delta_p, gak_standzeit_tage, doc_faktor = 0.0, 0.0, 0, 0.35
else:
    ozon_dosis = 0.0

if not is_ozon:
    with st.sidebar.expander("🧬 Nanofiltration (Membran)", expanded=True):
        recovery = st.slider("Permeatausbeute Y", 0.50, 0.95, 0.80, 0.05)
        recycle_fraction = st.slider("Konzentrat-Rückführung in die Biologie (Anteil)", 0.0, 0.90, 0.20, 0.05)
        mwco = st.slider("Trenngrenze MWCO bei 90 % Rückhalt (Da)", 50, 500, 200, 10)
        zeta_pot = st.slider("Zeta-Potential (mV)", -30.0, 0.0, -15.0, 1.0)
        sigma_poren = st.slider("Porenweitenverteilung (σ)", 0.1, 0.6, 0.3, 0.05)
        delta_p = st.slider("Membrandruck ΔP_NF (bar)", 4.0, 25.0, 10.0, 0.5)

    with st.sidebar.expander("🪨 GAK-Polishing & Konzentrat-GAK", expanded=False):
        aktivieren_gak = st.checkbox("Produktwasser-GAK aktivieren", value=False)
        aktivieren_konz_gak = st.checkbox("Konzentrat-GAK (Restkonzentrat-Behandlung & Vermischung)", value=True)
        
        if aktivieren_gak or aktivieren_konz_gak:
            gak_ebct = st.slider("Kontaktzeit EBCT (min)", 5.0, 30.0, 15.0, 1.0)
            gak_delta_p = st.slider("GAK-Druckverlust ΔP (bar)", 0.2, 2.0, 0.6, 0.1)
            gak_standzeit_tage = st.slider("Filter-Betriebstage (Erschöpfung)", 30, 365, 180, 30)
            doc_faktor = st.slider("Matrix-Faktor (DOC-Konkurrenz)", 0.1, 0.8, 0.35, 0.05)
        else:
            gak_ebct, gak_delta_p, gak_standzeit_tage, doc_faktor = 0.0, 0.0, 0, 0.35
else:
    recovery, recycle_fraction, mwco, zeta_pot, sigma_poren, delta_p = 0.8, 0.2, 200, -15.0, 0.3, 0.0

mwco_median = mwco / np.exp(1.2816 * sigma_poren)

with st.sidebar.expander("🧫 Vorklärung & Biologie", expanded=False):
    vkl_option = st.selectbox("Vorklärung Verweilzeit", ["0,75 - 1 h", "1,5 - 2 h", "> 2,5 h"], index=1)
    eta_vorklaerung = {"0,75 - 1 h": 0.30, "1,5 - 2 h": 0.35, "> 2,5 h": 0.40}[vkl_option]
    hrt_h = st.number_input("Aufenthaltszeit Belebung HRT (h)", 2.0, 48.0, 12.0, 1.0)
    srt = st.number_input("Schlammalter SRT (Tage)", 0.1, 100.0, 15.0, 1.0)
    x_ts = st.number_input("Trockensubstanz X_TS (kg/m³)", 0.5, 8.0, 3.5, 0.1)
    f_s_csb = st.number_input("Inerter gelöster CSB-Anteil f_S (DWA-A 131)", 0.0, 0.20, 0.05, 0.01)

with st.sidebar.expander("⚖️ Abwasserverordnung (AbwV)", expanded=False):
    groessenklasse = st.selectbox("Anlagengrößenklasse",
                                  ["Größenklasse 1", "Größenklasse 2", "Größenklasse 3", "Größenklasse 4",
                                   "Größenklasse 5"], index=3)

grenzwerte = {
    "Größenklasse 1": {"CSB": 150.0, "BSB5": 40.0},
    "Größenklasse 2": {"CSB": 110.0, "BSB5": 25.0},
    "Größenklasse 3": {"CSB": 90.0, "BSB5": 20.0},
    "Größenklasse 4": {"CSB": 90.0, "BSB5": 20.0},
    "Größenklasse 5": {"CSB": 75.0, "BSB5": 15.0}
}
grenzwert_csb = grenzwerte[groessenklasse]["CSB"]
grenzwert_bsb5 = grenzwerte[groessenklasse]["BSB5"]

ETA_PUMPE = 0.70
BETRIEBSSTUNDEN_JAHR = 8760
PH_SYS_FIX = 7.5
GAK_SCHUETTDICHTE_KG_M3 = 450.0
SZENARIO_ZIELWERT_PFOA_PFOS = 300.0

st.sidebar.markdown("---")
st.sidebar.subheader("📄 Dokumentation")

try:
    with open("Begleitdokument.pdf", "rb") as pdf_file:
        PDF_BYTE = pdf_file.read()
    
    st.sidebar.download_button(
        label="📥 Begleitdokument öffnen / laden",
        data=PDF_BYTE,
        file_name="Begleitdokument_Abwasser_Simulation.pdf",
        mime="application/pdf"
    )
except FileNotFoundError:
    st.sidebar.info("ℹ️ Begleitdokument derzeit nicht verfügbar.")


def build_process_svg(
    is_ozon,
    has_bypass,
    aktivieren_gak,
    aktivieren_konz_gak,
    recycle_fraction=0.0,
    bypass_fraction=0.0
):
    import html

    # Prozentwerte direkt aus den Sidebar-Einstellungen
    recycle_percent = recycle_fraction * 100
    konzentrat_percent = (1 - recycle_fraction) * 100
    bypass_percent = bypass_fraction * 100

    # Schwarz-Weiß-Farbschema
    stroke = "#222222"
    fill = "#ffffff"
    text_color = "#222222"
    muted = "#666666"

    parts = []

    parts.append("""
    <svg xmlns="http://www.w3.org/2000/svg"
         viewBox="0 0 1420 470"
         width="100%"
         height="100%"
         style="background:#ffffff;font-family:Arial,sans-serif">

    <defs>
      <marker id="arrow"
              markerWidth="10"
              markerHeight="10"
              refX="8"
              refY="3"
              orient="auto">
        <path d="M0,0 L0,6 L9,3 z" fill="#222222"/>
      </marker>
    </defs>
    """)

    def line(x1, y1, x2, y2, dashed=False, arrow=True):
        dash = 'stroke-dasharray="7,5"' if dashed else ""
        marker = 'marker-end="url(#arrow)"' if arrow else ""
        parts.append(
            f'<path d="M{x1},{y1} L{x2},{y2}" '
            f'fill="none" stroke="{stroke}" stroke-width="2" '
            f'{dash} {marker}/>'
        )

    def box(x, y, w, h, title, subtitle=None):
        parts.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" '
            f'rx="3" fill="{fill}" stroke="{stroke}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{x+w/2}" y="{y+30}" text-anchor="middle" '
            f'font-size="15" font-weight="bold" fill="{text_color}">'
            f'{html.escape(title)}</text>'
        )
        if subtitle:
            parts.append(
                f'<text x="{x+w/2}" y="{y+52}" text-anchor="middle" '
                f'font-size="12" fill="{muted}">'
                f'{html.escape(subtitle)}</text>'
            )

    def label(x, y, value, anchor="middle", size=12):
        parts.append(
            f'<text x="{x}" y="{y}" text-anchor="{anchor}" '
            f'font-size="{size}" fill="{text_color}">'
            f'{html.escape(str(value))}</text>'
        )

    def junction(x, y):
        parts.append(
            f'<circle cx="{x}" cy="{y}" r="4" fill="{stroke}"/>'
        )

    # ---------------------------------------------------------
    # Ozon-Variante
    # ---------------------------------------------------------
    if is_ozon:
        box(40, 120, 160, 70, "Vorbehandlung")
        box(300, 120, 160, 70, "Biologie")
        box(580, 120, 160, 70, "Ozonung")
        box(860, 120, 160, 70, "Nachbehandlung")

        line(200, 155, 300, 155)
        line(460, 155, 580, 155)
        line(740, 155, 860, 155)
        line(1020, 155, 1150, 155)

        label(40, 100, "Zulauf", "start")
        label(1150, 140, "Ablauf", "start")

    # ---------------------------------------------------------
    # Nanofiltrations-Variante
    # ---------------------------------------------------------
    else:
        # Hauptaggregate
        box(35, 115, 145, 70, "Vorbehandlung")
        box(300, 115, 155, 70, "Biologie")
        box(620, 115, 160, 70, "Nanofiltration")

        # Zulauf und Biologie
        line(180, 150, 300, 150)
        line(455, 150, 530, 150)
        junction(530, 150)

        label(35, 95, "Zulauf", "start")

        # Bypass: biologischer Ablauf umgeht die NF
        if has_bypass:
            line(530, 150, 530, 55, arrow=False)
            line(530, 55, 930, 55, arrow=False)
            line(930, 55, 930, 150)

            label(
                720, 42,
                f"Bypass: {bypass_percent:.1f} %"
            )

        # NF-Zulauf und Permeat
        line(530, 150, 620, 150)
        line(780, 150, 930, 150)
        junction(930, 150)

        label(850, 132, "NF-Permeat")

        # Produktwasser-GAK nach der Zusammenführung
        if aktivieren_gak:
            box(990, 115, 150, 70, "Produktwasser-GAK")
            line(930, 150, 990, 150)
            line(1140, 150, 1220, 150)
        else:
            line(930, 150, 1220, 150)

        # Konzentrat aus der NF
        line(700, 185, 700, 260)
        junction(700, 260)

        # Konzentrat-Rückführung zur Biologie
        if recycle_fraction > 0:
            line(700, 260, 700, 350, arrow=False)
            line(700, 350, 380, 350, arrow=False)
            line(380, 350, 380, 185)

            label(
                540, 340,
                f"Rückführung: {recycle_percent:.1f} %"
            )

        # Verbleibendes Konzentrat
        if konzentrat_percent > 0:
            line(700, 260, 900, 260)
            label(
                800, 245,
                f"Restkonzentrat: {konzentrat_percent:.1f} %"
            )

            if aktivieren_konz_gak:
                box(900, 230, 170, 70, "Konzentrat-GAK")
                line(1070, 265, 1220, 265)
                line(1220, 265, 1220, 150)
            else:
                line(900, 260, 900, 390)
                label(915, 385, "Konzentratabzug", "start")

        # Ablauf
        label(1240, 140, "Produktwasser", "start")
        line(1220, 150, 1360, 150)
        label(1360, 135, "Ablauf", "start")

    # Abschluss
    parts.append("</svg>")

    return "".join(parts)

# =========================================================
# 5. ZULAUFKONZENTRATIONEN
# =========================================================

st.header("🧪 Rohwasser-Zulaufkonzentrationen")
if "c_roh" not in st.session_state:
    st.session_state.c_roh = default_concentrations.copy()

if st.button("↺ Standardwerte wiederherstellen"):
    st.session_state.c_roh = default_concentrations.copy()
    for key in [k for k in st.session_state.keys() if str(k).startswith("inp_")]:
        del st.session_state[key]
    st.rerun()

sum_cols = st.columns(2)
for i, stoff in enumerate(["CSB", "BSB5"]):
    with sum_cols[i]:
        st.session_state.c_roh[stoff] = st.number_input(
            f"{stoff} (mg/L)", min_value=0.0,
            value=float(st.session_state.c_roh[stoff]), step=10.0, key=f"inp_{stoff}")

trace_stoffe = [s for s in default_concentrations if s not in ["CSB", "BSB5"]]
for start in range(0, len(trace_stoffe), 4):
    cols = st.columns(4)
    for col, stoff in zip(cols, trace_stoffe[start:start + 4]):
        with col:
            st.session_state.c_roh[stoff] = st.number_input(
                f"{stoff} (ng/L)", min_value=0.0,
                value=float(st.session_state.c_roh[stoff]), step=50.0, key=f"inp_{stoff}")

df = df_stoff.copy()
df["C_Rohwasser"] = df["Substanzname"].map(st.session_state.c_roh)

# =========================================================
# 6. SIMULATIONSLOGIK & VERMISCHUNGSLOGIK
# =========================================================

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

if is_ozon:
    fl, q_produkt = None, q_feed
    hrt_eff_h = hrt_h
else:
    fl = kreislauf_fluesse(q_feed, bypass_fraction, recovery, recycle_fraction)
    q_produkt = fl["q_prod"]
    hrt_eff_h = hrt_h * q_feed / fl["q_bio"]

def berechne_vorklaerung(row, eta):
    c = float(row["C_Rohwasser"])
    return c * (1.0 - eta) if row["Gruppe"] == "Summenparameter" else c

df["C_Vorklaerung"] = df.apply(lambda row: berechne_vorklaerung(row, eta_vorklaerung), axis=1)

def berechne_biologie(row, x_ts, srt, hrt_d):
    c_in = float(row["C_Vorklaerung"])
    kd = float(row["Kd_m3_kg"])
    k_bio = float(row["k_biol_L_gSS_d"])

    ress = x_ts * hrt_d / srt
    e_sorp = (kd * ress) / (1.0 + kd * ress)

    if row["Substanzname"] == "CSB":
        e_sorp = 0.0
        e_bio = 1.0 - f_s_csb
    elif row["Gruppe"] == "Summenparameter":
        e_bio = (k_bio * srt) / (1.0 + k_bio * srt) if k_bio > 0 else 0.0
    else:
        kxh = k_bio * x_ts * hrt_d
        e_bio = kxh / (1.0 + kxh) if k_bio > 0 else 0.0

    e_total = e_sorp + e_bio - e_sorp * e_bio
    e_total = float(np.clip(e_total, 0.0, 0.98))
    c_out = c_in * (1.0 - e_total)

    return pd.Series({
        "C_Bio_Ablauf": c_out,
        "Bio_Elimination_%": e_total * 100.0,
        "Bio_Sorptionsanteil_%": e_sorp * 100.0,
        "Bio_Abbauanteil_%": e_bio * 100.0
    })

bio_eff = df.apply(lambda row: berechne_biologie(row, x_ts, srt, hrt_eff_h / 24.0), axis=1)
df[["C_Bio_Ablauf", "Bio_Elimination_%", "Bio_Sorptionsanteil_%", "Bio_Abbauanteil_%"]] = bio_eff

bio_ref = df.apply(lambda row: berechne_biologie(row, x_ts, srt, hrt_h / 24.0), axis=1)
df["C_Bio_Ablauf_ref"] = bio_ref["C_Bio_Ablauf"]

def wende_ozon_an(row, dosis):
    substanz = row["Substanzname"]
    gruppe = row["Gruppe"]
    c_in = row["C_Bio_Ablauf"]

    if gruppe == "Summenparameter" or gruppe == "PFAS":
        return c_in, 0.0

    if gruppe == "PhAC":
        p = OZON_PARAMETER.get(substanz)
        ct_o3 = (dosis / 4.0) * 1.5e-3
        ct_oh = (dosis / 4.0) * 1.0e-10

        if p and p["k_o3"] is not None:
            restanteil = np.exp(-(p["k_o3"] * ct_o3 + p["k_oh"] * ct_oh))
            ozon_rem = 1.0 - restanteil
        else:
            ozon_rem = 0.85 * min(1.0, dosis / 4.0)

        ozon_rem = float(np.clip(ozon_rem, 0.0, 0.99))
        c_out = max(0.0, float(c_in * (1.0 - ozon_rem)))
        return c_out, ozon_rem * 100.0

    return c_in, 0.0

df["alpha"] = np.where(df["Ladung_z"] < 0, 1.0 / (1.0 + 10.0 ** (df["pKa"] - PH_SYS_FIX)), 0.0)
df["R_size"] = lognorm.cdf(df["MW"], s=sigma_poren, scale=mwco_median)
df.loc[df["Gruppe"] == "Summenparameter", "R_size"] = 0.0

def get_nf_retention(row):
    if row["Gruppe"] == "Summenparameter":
        return 0.0
    charge_factor = float(np.clip(abs(zeta_pot / 30.0) * row["alpha"] * abs(row["Ladung_z"]), 0.0, 0.95))
    p_size = float(np.clip(1.0 - row["R_size"], 0.0, 1.0))
    p_charge = float(np.clip(1.0 - charge_factor, 0.0, 1.0))
    return float(np.clip(1.0 - (p_size * p_charge), 0.0, 0.99))

df["R_NF"] = df.apply(get_nf_retention, axis=1)

def wende_gak_an(c_in, row, ebct, q_ref_m3h, betriebstage, doc_faktor):
    if row["Gruppe"] == "Summenparameter" or pd.isna(c_in) or c_in <= 0:
        return (c_in if not pd.isna(c_in) else 0.0), 0.0

    substanz = row["Substanzname"]
    gruppe = row["Gruppe"]
    kf = float(row["KF"])

    if kf <= 0.0:
        return c_in, 0.0

    q_ref = q_ref_m3h if q_ref_m3h > 0 else 1.0
    v_filter_m3 = q_ref * (ebct / 60.0)
    bed_volumes = (q_ref * 24.0 * betriebstage) / v_filter_m3

    if substanz in ["TFA", "PFBS"]:
        bv_50_base = 4000.0
    elif gruppe == "PFAS":
        bv_50_base = 15000.0
    else:
        bv_50_base = max(3000.0, kf * 600.0)

    bv_50_eff = bv_50_base * doc_faktor
    steepness = 3.5

    breakthrough_ratio = 1.0 / (1.0 + np.exp(-steepness * ((bed_volumes / bv_50_eff) - 1.0)))
    max_rem = 0.98 if gruppe == "PFAS" else 0.99
    rem = max_rem * (1.0 - breakthrough_ratio)
    rem = float(np.clip(rem, 0.0, 0.99))
    c_out = max(0.0, float(c_in * (1.0 - rem)))
    return c_out, rem * 100.0

result_rows = []
for _, row in df.iterrows():
    c_roh = float(row["C_Rohwasser"])
    c_vkl = float(row["C_Vorklaerung"])
    E = float(row["Bio_Elimination_%"]) / 100.0
    c_bio_ref = float(row["C_Bio_Ablauf_ref"])

    nf_elim_pct = gak_elim_pct = ozon_elim_pct = konz_gak_elim_pct = 0.0
    c_konz = np.nan
    c_konz_gac = np.nan
    f_bio = f_konz = f_gak = f_konz_gak = 0.0
    akkum = 1.0

    fracht_in = q_feed * c_vkl

    if is_ozon:
        c_bio = float(row["C_Bio_Ablauf"])
        c_final, ozon_elim_pct = wende_ozon_an(row, ozon_dosis)
    else:
        r_nf = float(row["R_NF"])
        nf_elim_pct = r_nf * 100.0
        E_rez = 0.0 if row["Substanzname"] == "CSB" else E
        
        c_bio = c_vkl * (1.0 - E)
        for _ in range(15):
            c_perm = c_bio * (1.0 - r_nf)
            c_konz = c_bio * (1.0 - recovery * (1.0 - r_nf)) / (1.0 - recovery)
            c_zu_bio = (q_feed * c_vkl + fl["q_rez"] * c_konz) / fl["q_bio"]
            c_bio = c_zu_bio * (1.0 - E)
            
        akkum = (c_bio / c_bio_ref) if c_bio_ref > 0 else 1.0
        c_nach_nf = (fl["q_perm"] * c_perm + fl["q_bypass"] * c_bio) / fl["q_prod"]

        if aktivieren_gak:
            c_prod_gac, gak_elim_pct = wende_gak_an(c_nach_nf, row, gak_ebct, q_produkt, gak_standzeit_tage, doc_faktor)
        else:
            c_prod_gac = c_nach_nf

        if aktivieren_konz_gak:
            c_konz_gac, konz_gak_elim_pct = wende_gak_an(c_konz, row, gak_ebct, fl["q_konz_aus"], gak_standzeit_tage, doc_faktor)
        else:
            c_konz_gac = c_konz

        q_ges_end = q_produkt + (fl["q_konz_aus"] if aktivieren_konz_gak else 0.0)
        if q_ges_end > 0:
            c_final = (q_produkt * c_prod_gac + (fl["q_konz_aus"] if aktivieren_konz_gak else 0.0) * c_konz_gac) / q_ges_end
        else:
            c_final = c_prod_gac

        if fracht_in > 0:
            f_bio = (E * (q_feed * c_vkl + fl["q_rez"] * c_konz)) / fracht_in * 100.0
            f_konz = (fl["q_konz_aus"] * c_konz_gac) / fracht_in * 100.0
            if aktivieren_konz_gak and c_konz > 0:
                f_konz_gak = (fl["q_konz_aus"] * (c_konz - c_konz_gac)) / fracht_in * 100.0
            f_gak = (q_produkt * (c_nach_nf - c_prod_gac)) / fracht_in * 100.0

    f_prod = (q_ges_end * c_final / fracht_in * 100.0) if (not is_ozon and fracht_in > 0) else ((q_feed * c_final / fracht_in * 100.0) if fracht_in > 0 else 0.0)

    if is_ozon:
        f_bilanz = np.nan
    else:
        f_bilanz = (f_bio + f_konz + f_konz_gak + f_gak + (q_ges_end * c_final / fracht_in * 100.0)) if fracht_in > 0 else 100.0

    konz_red = (1.0 - (c_final / c_roh)) * 100.0 if c_roh > 0 else np.nan

    result_rows.append({
        "Substanzname": row["Substanzname"], "Gruppe": row["Gruppe"], "Einheit": row["Einheit"],
        "C_Rohwasser": c_roh, "Bio_Elimination_%": row["Bio_Elimination_%"],
        "Akkumulation_Bio_x": akkum,
        "R_NF_%": nf_elim_pct, "Ozon_Elimination_%": ozon_elim_pct, "GAK_Elimination_%": gak_elim_pct,
        "C_Endwasser": c_final,
        "Fracht_Biologie_%": f_bio, "Fracht_Produkt_%": f_prod,
        "Fracht_Konzentrat_%": f_konz, "Fracht_Konz_GAK_%": f_konz_gak, "Fracht_GAK_%": f_gak,
        "Bilanzsumme_%": f_bilanz,
        "Konzentrationsreduktion_%": konz_red
    })

df_res = pd.DataFrame(result_rows)

# =========================================================
# 7. DASHBOARD: KERN-METRIKEN
# =========================================================

verfahren_label = "Bestandssystem (Ozonung)" if is_ozon else f"Nanofiltration (Bypass: {int(bypass_fraction*100)}%)"
st.header(f"🎯 Zielwertüberwachung: {verfahren_label}")

def get_concentration(stoff):
    values = df_res.loc[df_res["Substanzname"] == stoff, "C_Endwasser"].values
    return float(values[0]) if len(values) else np.nan

csb_end = get_concentration("CSB")
bsb_end = get_concentration("BSB5")
pfas_sum_end = df_res.loc[df_res["Gruppe"] == "PFAS", "C_Endwasser"].sum()
pfoa_pfos_end = float(np.nan_to_num(get_concentration("PFOA")) + np.nan_to_num(get_concentration("PFOS")))
avg_phac_elim = df_res.loc[df_res["Gruppe"] == "PhAC", "Konzentrationsreduktion_%"].mean()

pfas_ziel_ok = pfoa_pfos_end <= SZENARIO_ZIELWERT_PFOA_PFOS

z1, z2, z3, z4, z5 = st.columns(5)
z1.metric(f"CSB (max. {grenzwert_csb})", f"{csb_end:.1f} mg/L",
          delta="Eingehalten" if csb_end <= grenzwert_csb else "Überschritten",
          delta_color="normal" if csb_end <= grenzwert_csb else "inverse")
z2.metric(f"BSB5 (max. {grenzwert_bsb5})", f"{bsb_end:.1f} mg/L",
          delta="Eingehalten" if bsb_end <= grenzwert_bsb5 else "Überschritten",
          delta_color="normal" if bsb_end <= grenzwert_bsb5 else "inverse")
z3.metric(f"PFOA + PFOS (Szenario-Ziel {SZENARIO_ZIELWERT_PFOA_PFOS:.0f})", f"{pfoa_pfos_end:.1f} ng/L",
          delta="Unter Szenario-Zielwert" if pfas_ziel_ok else "Über Szenario-Zielwert",
          delta_color="normal" if pfas_ziel_ok else "inverse",
          help="Eigener Szenario-Zielwert des Projekts, kein gesetzlicher Grenzwert.")
z4.metric("PFAS-Summe", f"{pfas_sum_end:.1f} ng/L")
z5.metric("Ø Arzneimittel-Red.", f"{avg_phac_elim:.1f}%")

if not is_ozon:
    abw = (df_res["Bilanzsumme_%"] - 100.0).abs().max()
    if abw > 1e-4:
        st.warning(f"Massenbilanz-Abweichung: {abw:.3g} %")
    else:
        st.caption("✅ Frachtbilanz perfekt geschlossen (100.0 %).")

st.divider()

# =========================================================
# 8. KOSTENBILANZ & ENERGIEBEDARFS-BILANZ
# =========================================================

if is_ozon:
    st.header("💶 Betriebskostenbilanz (OPEX) – Ozonung")
    spez_opex_ozon = 0.018
    wassermenge_jahr = q_feed * BETRIEBSSTUNDEN_JAHR
    opex_gesamt_jahr = wassermenge_jahr * spez_opex_ozon

    k1, k2 = st.columns(2)
    k1.metric("Spezifische OPEX (Ozonung)", f"{spez_opex_ozon:.3f} €/m³", help="Basierend auf Daten Aachen-Sörs")
    k2.metric("Jährliche Betriebskosten", f"{opex_gesamt_jahr:,.0f} €/a".replace(",", "."))
    st.divider()

st.header("⚡ Energiebedarfsbilanz & Vergleich")
ENERGIE_BENCHMARK_AACHEN = 0.04824  # kWh/m³ (nur Ozonung, KA Aachen-Sörs 2018)

if not is_ozon:
    sec_permeat = delta_p / (36.0 * recovery * ETA_PUMPE)
    leistungsbedarf_nf_kw = fl["q_perm"] * sec_permeat
    jahresenergie_nf_kwh = leistungsbedarf_nf_kw * BETRIEBSSTUNDEN_JAHR

    spez_energie_gesamt_m3 = (leistungsbedarf_nf_kw / q_produkt) if q_produkt > 0 else 0.0

    e1, e2, e3 = st.columns(3)
    e1.metric("Spez. Energiebedarf (NF)", f"{sec_permeat:.2f} kWh/m³ Permeat")
    e2.metric("Elektrische Leistung (NF)", f"{leistungsbedarf_nf_kw:.1f} kW")
    e3.metric("Jahresenergie (NF)", f"{jahresenergie_nf_kwh:,.0f} kWh/a".replace(",", "."))

    if aktivieren_gak or aktivieren_konz_gak:
        sec_gak = gak_delta_p / (36.0 * ETA_PUMPE)
        q_gak_ref = q_produkt + (fl["q_konz_aus"] if aktivieren_konz_gak else 0.0)
        leistungsbedarf_gak_kw = q_gak_ref * sec_gak
        jahresenergie_gak_kwh = leistungsbedarf_gak_kw * BETRIEBSSTUNDEN_JAHR
        spez_energie_gesamt_m3 += (leistungsbedarf_gak_kw / q_produkt)
        st.metric("Gesamtenergie NF + GAK", f"{jahresenergie_nf_kwh + jahresenergie_gak_kwh:,.0f} kWh/a".replace(",", "."))

        bv_wechsel = 1440.0 * gak_standzeit_tage / gak_ebct
        gak_spez_verbrauch = GAK_SCHUETTDICHTE_KG_M3 * 1000.0 / bv_wechsel
        gak_jahresbedarf_t = gak_spez_verbrauch * q_gak_ref * BETRIEBSSTUNDEN_JAHR / 1.0e6
        g1, g2 = st.columns(2)
        g1.metric("Spez. Kohleverbrauch (berechnet)", f"{gak_spez_verbrauch:.1f} g/m³")
        g2.metric("Kohlebedarf pro Jahr", f"{gak_jahresbedarf_t:.1f} t/a")

    st.markdown("---")
    st.subheader("📊 Energetischer Vergleich zur Ozonung (Großmaßstab)")
    c_ben1, c_ben2 = st.columns(2)
    c_ben1.metric("Spez. Energiebedarf Advanced System", f"{spez_energie_gesamt_m3:.3f} kWh/m³")
    c_ben2.metric("Vergleichswert Ozonung KA Aachen-Sörs (2018)", f"{ENERGIE_BENCHMARK_AACHEN:.3f} kWh/m³")

else:
    jahresenergie_bestand_kwh = q_feed * ENERGIE_BENCHMARK_AACHEN * BETRIEBSSTUNDEN_JAHR
    e1, e2 = st.columns(2)
    e1.metric("Energiebedarf Ozonung (KA Aachen-Sörs 2018)", f"{ENERGIE_BENCHMARK_AACHEN:.3f} kWh/m³")
    e2.metric("Jahresenergiebedarf Ozonung", f"{jahresenergie_bestand_kwh:,.0f} kWh/a".replace(",", "."))

st.divider()

# =========================================================
# 9. TABELLE & DIAGRAMME
# =========================================================

st.header("📊 Detaillierte Stoffbilanz")

display_cols = ["Substanzname", "Gruppe", "Einheit", "C_Rohwasser", "Bio_Elimination_%"]
if is_ozon:
    display_cols.append("Ozon_Elimination_%")
else:
    display_cols.extend(["R_NF_%"])
    if aktivieren_gak:
        display_cols.append("GAK_Elimination_%")
display_cols.extend(["C_Endwasser", "Konzentrationsreduktion_%"])

display_df = df_res[display_cols].copy()
for col in display_df.columns:
    if col not in ["Substanzname", "Gruppe", "Einheit"]:
        display_df[col] = display_df[col].round(2)

st.dataframe(display_df, use_container_width=True, hide_index=True)

st.subheader("📈 Verursacher- & Stufen-Analyse: Fracht- und Eliminationsaufteilung")

if not is_ozon:
    stack_cols = ["Substanzname", "Fracht_Biologie_%", "Fracht_Konzentrat_%", "Fracht_Produkt_%"]
    if aktivieren_konz_gak:
        stack_cols.append("Fracht_Konz_GAK_%")
    if aktivieren_gak:
        stack_cols.append("Fracht_GAK_%")

    stack_df = df_res.loc[df_res["Gruppe"] != "Summenparameter", stack_cols].melt(
        id_vars=["Substanzname"],
        var_name="Prozessstufe",
        value_name="Anteil_Prozent"
    )

    stage_mapping = {
        "Fracht_Biologie_%": "1. Biologie (Abbau / Sorption)",
        "Fracht_Konzentrat_%": "2. Konzentrat-Restfracht (Endwasser)",
        "Fracht_Konz_GAK_%": "3. Konzentrat-GAK Adsorption",
        "Fracht_GAK_%": "4. Produktwasser-GAK Adsorption",
        "Fracht_Produkt_%": "5. Verbleib im Permeat/Endwasser"
    }
    stack_df["Prozessstufe"] = stack_df["Prozessstufe"].map(stage_mapping)

    chart = alt.Chart(stack_df).mark_bar().encode(
        x=alt.X("Substanzname:N", title="Substanz", sort=None),
        y=alt.Y("Anteil_Prozent:Q", title="Anteil am Rohwasser-Eintrag (%)", stack="zero"),
        color=alt.Color("Prozessstufe:N", title="Verbleib / Elimination über Stufe", scale=alt.Scale(scheme="greys")),
        tooltip=["Substanzname", "Prozessstufe", alt.Tooltip("Anteil_Prozent:Q", format=".1f")]
    ).properties(height=400)
    st.altair_chart(chart, use_container_width=True)

else:
    ozon_stack_rows = []
    for _, row in df_res.loc[df_res["Gruppe"] != "Summenparameter"].iterrows():
        c_roh = float(row["C_Rohwasser"])
        c_bio = float(row.get("C_Bio_Ablauf", row["C_Endwasser"]))
        c_end = float(row["C_Endwasser"])
        
        if c_roh > 0:
            bio_pct = max(0.0, (1.0 - (c_bio / c_roh))) * 100.0
            ozon_pct = max(0.0, ((c_bio - c_end) / c_roh)) * 100.0
            end_pct = max(0.0, (c_end / c_roh)) * 100.0
        else:
            bio_pct, ozon_pct, end_pct = 0.0, 0.0, 0.0

        ozon_stack_rows.append({"Substanzname": row["Substanzname"], "Prozessstufe": "1. Biologischer Abbau", "Anteil_Prozent": bio_pct})
        ozon_stack_rows.append({"Substanzname": row["Substanzname"], "Prozessstufe": "2. Ozonung", "Anteil_Prozent": ozon_pct})
        ozon_stack_rows.append({"Substanzname": row["Substanzname"], "Prozessstufe": "3. Verbleib im Endwasser", "Anteil_Prozent": end_pct})
        
    stack_df = pd.DataFrame(ozon_stack_rows)
    chart = alt.Chart(stack_df).mark_bar().encode(
        x=alt.X("Substanzname:N", title="Substanz", sort=None),
        y=alt.Y("Anteil_Prozent:Q", title="Anteil am Rohwasser (%)", stack="zero"),
        color=alt.Color("Prozessstufe:N", title="Prozessstufe", scale=alt.Scale(scheme="greys")),
        tooltip=["Substanzname", "Prozessstufe", alt.Tooltip("Anteil_Prozent:Q", format=".1f")]
    ).properties(height=400)
    st.altair_chart(chart, use_container_width=True)
