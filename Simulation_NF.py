# -*- coding: utf-8 -*-
"""
Integrierte Abwasser-Simulation: Verfahrensvergleich & Pilotanlagen-Vorbereitung
1. Bestandssystem: Biologie + Ozonung
2. Advanced System: Nanofiltration mit Konzentrat-Rückführung, Produkt-GAK und Konzentrat-GAK

Version 5 (Korrekturen gegenüber v4/Eigenversion):
- Kreislauf wieder geschlossen gelöst (statt 15 Iterationen), inkl. E_rez für CSB
- Frachtbilanz ohne Doppelzählung des Konzentrats bei aktiver Konzentrat-GAK
- Konzentrat-GAK mit eigenem Matrix-Faktor (DOC im Konzentrat aufkonzentriert)
- GAK-Elimination wahlweise als Mittelwert über die Standzeit statt Momentaufnahme
- Energie/Kohle: nur die tatsächlich durchströmten Ströme, Bezug auf Q_Zulauf
- Fließbild: saubere Verzweigung, Mischpunkt hinter Konz-GAK, Bypass vor der Produkt-GAK
- Ozon-Diagramm: C_Bio_Ablauf wird wieder übergeben
- Passwort direkt im Skript hinterlegt
"""

import hmac
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

# Festes Passwort direkt im Skript hinterlegen
PASSWORD = "Passwort"  # Hier Wunschpasswort eintragen

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
    Willkommen bei der **Abwasser-Simulation**! Mit diesem Werkzeug kannst du verschiedene Verfahren zur Abwasserreinigung (Spurenstoff- und PFAS-Elimination) direkt miteinander vergleichen.

    ---

    ### So bedienst du die Simulation

    1. **Verfahren auswählen (Seitenleiste links):**
       * Wähle oben links zwischen dem **Bestandssystem** (Klassische Ozonung) und dem **Advanced System** (Nanofiltration/NF).
    
    2. **Parameter anpassen:**
       * **Hydraulik:** Stelle die Wassermenge ($Q$) und evtl. einen Bypass ein.
       * **Reinigungsstufen:** Schalte nach Bedarf die *Produktwasser-GAK* oder die *Konzentrat-GAK* an oder aus.
       * **Stoffkonzentrationen:** Weiter unten kannst du die Zulaufwerte für Schadstoffe (z. B. PFAS, Medikamente) individuell anpassen.

    3. **Ergebnisse ablesen:**
       * Die Simulation rechnet alle Frachten und Konzentrationen **automatisch in Echtzeit** um.
       * Achte ganz oben im Dashboard auf die **Ampel-Indikatoren**: Sie zeigen dir sofort, ob gesetzliche Grenzwerte oder Projekt-Szenarioziele eingehalten werden.

    ---

    ### Wie funktioniert das Advanced-Konzept? (Kurzerklärt)

    * **Das Membran-Prinzip (NF):** Die Nanofiltration hält Schadstoffe wie PFAS extrem effizient zurück. Dabei entstehen zwei Ströme: gereinigtes *Permeat* (Hauptstrom) und aufkonzentriertes *Konzentrat*.
    * **Kreislauf-Rückführung:** Ein Teil des Konzentrats wird zurück in die Biologie geleitet, um den Abbau zu maximieren.
    * **Konzentrat-GAK & Vermischung:** Das restliche Konzentrat läuft durch einen Aktivkohlefilter (GAK) und wird anschließend wieder sauber mit dem Hauptstrom vermischt.

    > ℹ️ **Tipp:** Für weitere Informationen, lade das Begleitdokument ganz unten in der Seitenleiste herunter.
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

gak_ebct, gak_delta_p, gak_standzeit_tage, doc_faktor = 0.0, 0.0, 0, 0.35
gak_mittelwert, r_doc = True, 0.90

if is_ozon:
    with st.sidebar.expander("🧪 Ozonung", expanded=True):
        ozon_dosis = st.slider("Ozon-Dosis (mg O₃ / L)", 1.0, 10.0, 4.0, 0.5)
    aktivieren_gak = False
    aktivieren_konz_gak = False
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
            gak_standzeit_tage = st.slider("Filter-Standzeit bis zum Kohletausch (Tage)", 30, 365, 180, 30)
            gak_mittelwert = st.checkbox(
                "Mittlere Elimination über die Standzeit (statt Momentaufnahme am Standzeitende)", value=True,
                help="Realistisch ist der Mittelwert über den Zyklus bis zum Kohletausch.")
            doc_faktor = st.slider("Matrix-Faktor Produktwasser (DOC-Konkurrenz)", 0.1, 0.8, 0.35, 0.05)
            if aktivieren_konz_gak:
                r_doc = st.slider("DOC-Rückhalt der NF (für Konzentrat-GAK)", 0.50, 0.99, 0.90, 0.01,
                                  help="Bestimmt, wie stark der DOC im Konzentrat aufkonzentriert wird.")
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
    if srt * 24.0 <= hrt_h:
        st.warning("Schlammalter muss größer als die Aufenthaltszeit sein.")

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

# =========================================================
# 4. SIMULATIONSLOGIK: HYDRAULIK
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

# =========================================================
# 5. DYNAMISCHES PROZESSFLIEßBILD (SCHWARZ-WEISS)
# =========================================================

st.subheader("🗺️ Anlagen-Prozessfließbild")

def build_process_svg(is_ozon, has_bypass, aktivieren_gak, aktivieren_konz_gak, fl, q_feed):
    parts = []

    def text(x, y, label, cls="label", anchor="middle"):
        return f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{escape(label)}</text>'

    def line(points, arrow=False, cls="pipe"):
        coords = " ".join(f"{x},{y}" for x, y in points)
        marker = ' marker-end="url(#arrow)"' if arrow else ""
        return f'<polyline points="{coords}" class="{cls}"{marker}/>'

    def box(x, y, w, h, label):
        return f"""
        <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5"/>
        <text x="{x + w / 2}" y="{y + h / 2 + 5}" class="box-label" text-anchor="middle">{escape(label)}</text>
        """

    def junction(x, y):
        return f'<circle cx="{x}" cy="{y}" r="3.5" class="junction"/>'

    parts.append("""
    <svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 1280 310" preserveAspectRatio="xMidYMid meet" role="img">
      <defs>
        <marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="userSpaceOnUse">
          <path d="M0,0 L8,4 L0,8 Z" class="arrowhead"/>
        </marker>
      </defs>
      <style>
        text { font-family: Arial, sans-serif; font-size: 13px; fill: #000; }
        .subtitle { font-size: 11px; fill: #444; }
        .box-label { font-size: 14px; font-weight: bold; }
        .label { font-size: 11px; fill: #222; }
        .pipe { fill: none; stroke: #000; stroke-width: 1.8; stroke-linecap: round; stroke-linejoin: round; }
        .arrowhead { fill: #000; stroke: none; }
        rect { fill: #fff; stroke: #000; stroke-width: 1.8; }
        .junction { fill: #000; stroke: none; }
      </style>
    """)

    parts.append(text(150, 15, "Vorklärung", "subtitle"))
    parts.append(text(12, 91, "Kommunales", anchor="start"))
    parts.append(text(12, 108, "Abwasser", anchor="start"))
    parts.append(text(12, 126, f"Q = {q_feed:.0f} m³/h", "label", "start"))
    parts.append(line([(100, 100), (128, 100)], arrow=True))
    parts.append('<rect x="128" y="72" width="45" height="56" rx="3"/>')
    for i in range(5):
        x = 133 + i * 7
        parts.append(f'<line x1="{x}" y1="77" x2="{x + 12}" y2="123" class="pipe"/>')
    parts.append(line([(173, 100), (225, 100)], arrow=True))

    parts.append(text(300, 15, "Biologische Reinigung", "subtitle"))
    parts.append(box(225, 65, 150, 70, "Biologie"))
    parts.append(line([(300, 135), (300, 175)], arrow=True))
    parts.append(text(310, 168, "Schlammabzug", "label", "start"))

    if is_ozon:
        parts.append(text(625, 15, "Ozonung", "subtitle"))
        parts.append(line([(375, 100), (550, 100)], arrow=True))
        parts.append(box(550, 65, 150, 70, "Ozonung"))
        parts.append(line([(700, 100), (770, 100)], arrow=True))
        parts.append(text(810, 94, "Gereinigtes", anchor="start"))
        parts.append(text(810, 111, "Endwasser", anchor="start"))
    else:
        parts.append(text(600, 15, "Nanofiltration", "subtitle"))

        parts.append(line([(375, 100), (535, 100)], arrow=True))
        if has_bypass:
            parts.append(junction(430, 100))
            parts.append(line([(430, 100), (430, 45), (700, 45), (700, 100)], arrow=True))
            parts.append(text(565, 36, f"Bypass ({fl['q_bypass']:.1f} m³/h)", "subtitle"))
            parts.append(junction(700, 100))

        parts.append("""
        <rect x="535" y="60" width="120" height="80" rx="5"/>
        <text x="560" y="105" class="box-label" text-anchor="middle">NF</text>
        <line x1="535" y1="61" x2="654" y2="139" stroke="#000" stroke-width="1.8"/>
        """)

        current_x = 655
        if aktivieren_gak:
            parts.append(text(795, 40, "Produktwasser-GAK", "subtitle"))
            parts.append(line([(current_x, 100), (745, 100)], arrow=True))
            parts.append(box(745, 70, 100, 60, "GAK"))
            current_x = 845

        mix_x = 900 if aktivieren_konz_gak else current_x + 50
        parts.append(line([(current_x, 100), (mix_x, 100)], arrow=False))

        parts.append(line([(595, 140), (595, 215)], arrow=False))
        parts.append(junction(595, 215))
        parts.append(line([(595, 215), (185, 215), (185, 125), (225, 125)], arrow=True))
        parts.append(text(390, 206, f"Rückführung in Biologie ({fl['q_rez']:.1f} m³/h)", "label"))

        if aktivieren_konz_gak:
            parts.append(text(785, 222, "Konzentrat-GAK", "subtitle"))
            parts.append(line([(595, 215), (595, 260), (730, 260)], arrow=True))
            parts.append(box(730, 230, 110, 60, "Konz-GAK"))
            parts.append(line([(840, 260), (mix_x, 260), (mix_x, 100)], arrow=True))
            parts.append(junction(mix_x, 100))
            parts.append(text(mix_x + 8, 200, "Gereinigtes Konzentrat", "label", "start"))
            parts.append(text(610, 252, f"{fl['q_konz_aus']:.1f} m³/h", "label", "start"))
        else:
            parts.append(line([(595, 215), (595, 260), (760, 260)], arrow=True))
            parts.append(text(770, 264, f"Konzentrat-Abzug ({fl['q_konz_aus']:.1f} m³/h)", "label", "start"))

        parts.append(line([(mix_x, 100), (mix_x + 70, 100)], arrow=True))
        parts.append(text(mix_x + 85, 94, "Gereinigtes", anchor="start"))
        parts.append(text(mix_x + 85, 111, "Endwasser", anchor="start"))

    parts.append("</svg>")
    return "".join(parts)

svg_schema = build_process_svg(
    is_ozon=is_ozon,
    has_bypass=has_bypass,
    aktivieren_gak=aktivieren_gak,
    aktivieren_konz_gak=aktivieren_konz_gak,
    fl=fl,
    q_feed=q_feed
)

html_schema = f"""
<div style="background: #ffffff; padding: 16px; border-radius: 8px; border: 1px solid #000000; width: 100%; box-sizing: border-box; overflow-x: auto;">
    <div style="min-width: 820px; max-width: 1280px; margin: 0 auto;">{svg_schema}</div>
</div>
"""
components.html(html_schema, height=370, scrolling=False)
st.divider()

# =========================================================
# 6. ZULAUFKONZENTRATIONEN
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
# 7. SIMULATIONSLOGIK: STOFFE
# =========================================================

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

def kreislauf_stoffstrom(c_vkl, E, R, Y, r, q_feed, fl, E_rez=None):
    if E_rez is None:
        E_rez = E
    Y = min(Y, 0.99)
    cpf = 1.0 - float(np.clip(R, 0.0, 1.0))
    q_ratio = q_feed / fl["q_bio"]
    nenner = 1.0 - (1.0 - E_rez) * r * fl["a"] * (1.0 - Y * cpf)

    c_bio = (1.0 - E) * c_vkl * q_ratio / nenner
    c_perm = c_bio * cpf
    c_konz = c_bio * (1.0 - Y * cpf) / (1.0 - Y)
    c_prod = (fl["q_perm"] * c_perm + fl["q_bypass"] * c_bio) / fl["q_prod"]
    c_zu_bio = (q_feed * c_vkl + fl["q_rez"] * c_konz) / fl["q_bio"]
    return c_zu_bio, c_bio, c_perm, c_konz, c_prod

def wende_gak_an(c_in, row, ebct, betriebstage, matrix_faktor, mittelwert=True):
    if row["Gruppe"] == "Summenparameter" or pd.isna(c_in) or c_in <= 0:
        return (c_in if not pd.isna(c_in) else 0.0), 0.0

    substanz = row["Substanzname"]
    gruppe = row["Gruppe"]
    kf = float(row["KF"])

    if kf <= 0.0:
        return c_in, 0.0

    bed_volumes = 1440.0 * betriebstage / ebct

    if substanz in ["TFA", "PFBS"]:
        bv_50_base = 4000.0
    elif gruppe == "PFAS":
        bv_50_base = 15000.0
    else:
        bv_50_base = max(3000.0, kf * 600.0)

    bv_50_eff = bv_50_base * matrix_faktor
    steepness = 3.5
    max_rem = 0.98 if gruppe == "PFAS" else 0.99

    def rem_bei_bv(bv):
        durchbruch = 1.0 / (1.0 + np.exp(-steepness * ((bv / bv_50_eff) - 1.0)))
        return max_rem * (1.0 - durchbruch)

    if mittelwert:
        rem = float(np.mean(rem_bei_bv(np.linspace(0.0, bed_volumes, 400))))
    else:
        rem = float(rem_bei_bv(bed_volumes))

    rem = float(np.clip(rem, 0.0, 0.99))
    c_out = max(0.0, float(c_in * (1.0 - rem)))
    return c_out, rem * 100.0

if (not is_ozon) and aktivieren_konz_gak:
    y_c = min(recovery, 0.99)
    aufkonz_doc = (1.0 - y_c * (1.0 - r_doc)) / (1.0 - y_c)
    doc_faktor_konz = max(0.02, doc_faktor / aufkonz_doc)
else:
    aufkonz_doc = 1.0
    doc_faktor_konz = doc_faktor

# =========================================================
# 8. HAUPTSCHLEIFE
# =========================================================

result_rows = []
for _, row in df.iterrows():
    c_roh = float(row["C_Rohwasser"])
    c_vkl = float(row["C_Vorklaerung"])
    E = float(row["Bio_Elimination_%"]) / 100.0
    c_bio_ref = float(row["C_Bio_Ablauf_ref"])

    nf_elim_pct = gak_elim_pct = ozon_elim_pct = konz_gak_elim_pct = 0.0
    f_bio = f_gak = f_konz_gak = f_abzug = f_end = 0.0
    akkum = 1.0
    fracht_in = q_feed * c_vkl

    if is_ozon:
        c_bio = float(row["C_Bio_Ablauf"])
        c_final, ozon_elim_pct = wende_ozon_an(row, ozon_dosis)
        c_nach_nf = c_bio  # Für einheitliche Variable im Profil
    else:
        r_nf = float(row["R_NF"])
        nf_elim_pct = r_nf * 100.0
        E_rez = 0.0 if row["Substanzname"] == "CSB" else E

        c_zu_bio, c_bio, c_perm, c_konz, c_nach_nf = kreislauf_stoffstrom(
            c_vkl, E, r_nf, recovery, recycle_fraction, q_feed, fl, E_rez=E_rez)
        akkum = (c_bio / c_bio_ref) if c_bio_ref > 0 else 1.0

        if aktivieren_gak:
            c_prod_gac, gak_elim_pct = wende_gak_an(c_nach_nf, row, gak_ebct, gak_standzeit_tage,
                                                    doc_faktor, gak_mittelwert)
        else:
            c_prod_gac = c_nach_nf

        q_k = fl["q_konz_aus"]
        if aktivieren_konz_gak:
            c_konz_gac, konz_gak_elim_pct = wende_gak_an(c_konz, row, gak_ebct, gak_standzeit_tage,
                                                         doc_faktor_konz, gak_mittelwert)
            q_ges_end = q_produkt + q_k
            c_final = (q_produkt * c_prod_gac + q_k * c_konz_gac) / q_ges_end
        else:
            c_konz_gac = c_konz
            q_ges_end = q_produkt
            c_final = c_prod_gac

        if fracht_in > 0:
            f_bio = (E * q_feed * c_vkl + E_rez * fl["q_rez"] * c_konz) / fracht_in * 100.0
            f_gak = q_produkt * (c_nach_nf - c_prod_gac) / fracht_in * 100.0
            if aktivieren_konz_gak:
                f_konz_gak = q_k * (c_konz - c_konz_gac) / fracht_in * 100.0
            else:
                f_abzug = q_k * c_konz / fracht_in * 100.0
            f_end = q_ges_end * c_final / fracht_in * 100.0

    if is_ozon:
        f_bilanz = np.nan
    else:
        f_bilanz = (f_bio + f_gak + f_konz_gak + f_abzug + f_end) if fracht_in > 0 else 100.0

    konz_red = (1.0 - (c_final / c_roh)) * 100.0 if c_roh > 0 else np.nan

    result_rows.append({
        "Substanzname": row["Substanzname"], "Gruppe": row["Gruppe"], "Einheit": row["Einheit"],
        "C_Rohwasser": c_roh, "Bio_Elimination_%": row["Bio_Elimination_%"],
        "C_Bio_Ablauf": c_bio, "C_Nach_NF": c_nach_nf, "Akkumulation_Bio_x": akkum,
        "R_NF_%": nf_elim_pct, "Ozon_Elimination_%": ozon_elim_pct,
        "GAK_Elimination_%": gak_elim_pct, "Konz_GAK_Elimination_%": konz_gak_elim_pct,
        "C_Endwasser": c_final,
        "Fracht_Biologie_%": f_bio, "Fracht_Prod_GAK_%": f_gak, "Fracht_Konz_GAK_%": f_konz_gak,
        "Fracht_Konz_Abzug_%": f_abzug, "Fracht_Endwasser_%": f_end,
        "Bilanzsumme_%": f_bilanz,
        "Konzentrationsreduktion_%": konz_red
    })

df_res = pd.DataFrame(result_rows)

# =========================================================
# 9. DASHBOARD: KERN-METRIKEN
# =========================================================

verfahren_label = "Bestandssystem (Ozonung)" if is_ozon else f"Nanofiltration (Bypass: {int(bypass_fraction * 100)} %)"
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
    if abw > 1e-6:
        st.error(f"Massenbilanz nicht geschlossen (max. Abweichung {abw:.3g} %)")
    else:
        st.caption("✅ Frachtbilanz geschlossen (Biologie + GAK-Stufen + Konzentrat-Abzug + Endwasser = 100 %).")
    if not aktivieren_konz_gak:
        st.caption("ℹ️ Ohne Konzentrat-GAK verlässt das Konzentrat die Anlage als Abzug.")

st.divider()

# =========================================================
# 10. MENGENSTRÖME (NF)
# =========================================================

if not is_ozon:
    st.header("🔁 Mengenströme")
    h1, h2, h3, h4 = st.columns(4)
    h1.metric("Biologie-Zulauf Q_bio", f"{fl['q_bio']:.1f} m³/h")
    h2.metric("NF-Zulauf Q_NF", f"{fl['q_nf']:.1f} m³/h")
    h3.metric("Rücklaufstrom Q_rez", f"{fl['q_rez']:.1f} m³/h")
    h4.metric("Ausgeschl. Konzentrat", f"{fl['q_konz_aus']:.1f} m³/h")
    st.divider()


    # =========================================================
# 11. TABELLE & LINIEN-KONZENTRATIONSPROFIL
# =========================================================

st.header("📊 Detaillierte Stoffbilanz & Verlauf")

display_cols = ["Substanzname", "Gruppe", "Einheit", "C_Rohwasser", "Bio_Elimination_%"]
if is_ozon:
    display_cols.append("Ozon_Elimination_%")
else:
    display_cols.extend(["R_NF_%"])
    if aktivieren_gak:
        display_cols.append("GAK_Elimination_%")
    if aktivieren_konz_gak:
        display_cols.append("Konz_GAK_Elimination_%")
display_cols.extend(["C_Endwasser", "Konzentrationsreduktion_%"])

display_df = df_res[display_cols].copy()
for col in display_df.columns:
    if col not in ["Substanzname", "Gruppe", "Einheit"]:
        display_df[col] = display_df[col].round(2)

st.dataframe(display_df, use_container_width=True, hide_index=True)

st.subheader("📈 Linien-Konzentrationsprofil über die Reinigungsstufen")

verfügbare_substanzen = df_res.loc[df_res["Gruppe"] != "Summenparameter", "Substanzname"].tolist()
ausgewählte_substanzen = st.multiselect(
    "Substanzen für das Profil auswählen:",
    options=verfügbare_substanzen,
    default=["PFOA", "Diclofenac", "Carbamazepin", "Ibuprofen"] if "Diclofenac" in verfügbare_substanzen else verfügbare_substanzen[:3]
)

if ausgewählte_substanzen:
    profile_rows = []
    
    for _, row in df_res[df_res["Substanzname"].isin(ausgewählte_substanzen)].iterrows():
        substanz = row["Substanzname"]
        c_roh = float(row["C_Rohwasser"])
        
        if c_roh <= 0:
            continue
            
        c_vkl_rel = (float(row.get("C_Vorklaerung", c_roh)) / c_roh) * 100.0
        c_bio_rel = (float(row.get("C_Bio_Ablauf", c_roh)) / c_roh) * 100.0
        c_end_rel = (float(row["C_Endwasser"]) / c_roh) * 100.0
        
        if is_ozon:
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "1. Rohwasser", "Rel_Konzentration": 100.0})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "2. Vorklärung", "Rel_Konzentration": c_vkl_rel})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "3. Biologie", "Rel_Konzentration": c_bio_rel})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "4. Endwasser (Ozon)", "Rel_Konzentration": c_end_rel})
        else:
            c_nach_nf_rel = (float(row.get("C_Nach_NF", row["C_Endwasser"])) / c_roh) * 100.0
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "1. Rohwasser", "Rel_Konzentration": 100.0})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "2. Vorklärung", "Rel_Konzentration": c_vkl_rel})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "3. Biologie", "Rel_Konzentration": c_bio_rel})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "4. Nach Nanofiltration", "Rel_Konzentration": c_nach_nf_rel})
            profile_rows.append({"Substanz": substanz, "Prozessstufe": "5. Gereinigtes Endwasser", "Rel_Konzentration": c_end_rel})

    profile_df = pd.DataFrame(profile_rows)

    if not profile_df.empty:
        # HINWEIS: Hier von ":N" auf ":O" (Ordinal) geändert, da Linien eine Reihenfolge benötigen!
        base_chart = alt.Chart(profile_df).encode(
            x=alt.X("Prozessstufe:O", title="Verfahrensstufe"),
            y=alt.Y("Rel_Konzentration:Q", title="Verbleibende Konzentration (% vom Rohwasser)", scale=alt.Scale(domain=[0, 105])),
            color=alt.Color("Substanz:N", title="Substanz", scale=alt.Scale(scheme="category10")),
            tooltip=[
                alt.Tooltip("Substanz:N", title="Substanz"),
                alt.Tooltip("Prozessstufe:O", title="Stufe"),
                alt.Tooltip("Rel_Konzentration:Q", title="Verbleibend", format=".1f", suffix="%")
            ]
        )

        lines = base_chart.mark_line(strokeWidth=2.5)
        points = base_chart.mark_circle(size=70)

        chart = (lines + points).properties(height=420)
        st.altair_chart(chart, use_container_width=True)
    else:
        st.info("Bitte wähle mindestens eine Substanz aus.")
else:
    st.info("Keine Substanzen ausgewählt.")

