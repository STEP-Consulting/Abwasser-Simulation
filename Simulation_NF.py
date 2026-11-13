# -*- coding: utf-8 -*-
"""
Abwasser-Simulation: Verfahrensvergleich & Pilotanlagen-Vorbereitung
Oberfläche (Streamlit). Der Rechenkern liegt in simulation.py.

Start:  streamlit run app.py
Passwort: in .streamlit/secrets.toml  ->  password = "..."    oder Umgebungsvariable APP_PASSWORD
"""

import hmac
import json
import os
from html import escape

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

import simulation as sim

# =========================================================
# 1. KONFIGURATION & PASSWORTSCHUTZ
# =========================================================

st.set_page_config(
    page_title="Abwasser-Simulation: Verfahrensvergleich",
    page_icon="💧",
    layout="wide",
)


def hole_passwort():
    try:
        return st.secrets["password"]
    except Exception:
        return os.environ.get("APP_PASSWORD")


def check_password():
    """True, wenn das Passwort korrekt eingegeben wurde."""
    if st.session_state.get("password_correct", False):
        return True

    soll = hole_passwort()
    st.subheader("🔒 Passwortgeschützter Bereich")
    if not soll:
        st.error("Es ist kein Passwort konfiguriert. Bitte `password` in `.streamlit/secrets.toml` "
                 "(bzw. in den Secrets der Hosting-Plattform) oder die Umgebungsvariable "
                 "`APP_PASSWORD` setzen.")
        return False

    pwd = st.text_input("Bitte gib das Passwort ein, um die Simulation zu starten:", type="password")
    if st.button("Anmelden"):
        if hmac.compare_digest(pwd.encode("utf-8"), str(soll).encode("utf-8")):
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
    Mit diesem Werkzeug kannst du Verfahren zur Spurenstoff- und PFAS-Entfernung vergleichen
    (Bestandssystem **Ozonung** gegenüber **Nanofiltration**, optional mit Aktivkohle).

    **So bedienst du die Simulation**

    1. **Konzept wählen** (Seitenleiste): bestimmt die *Einzelansicht* und das Fließbild.
    2. **Parameter anpassen:** Alle Einstellungen in der Seitenleiste gelten auch für die
       Registerkarten *Vergleich* und *Sensitivität*.
    3. **Zulaufkonzentrationen** unter dem Fließbild anpassen.
    4. **Ergebnisse** in den drei Registerkarten ansehen:
       * *Einzelkonzept:* Kennzahlen, Mengenströme, Energie, Stoffbilanz, Diagramme.
       * *Vergleich:* Ozonung und NF nebeneinander mit denselben Eingaben.
       * *Sensitivität:* Welche Parameter bewegen das Ergebnis am stärksten?

    **Wichtige Hinweise**

    * **Produktwasser** = Permeat + Bypass (ggf. nach Produktwasser-GAK).
      **Endwasser** = Produktwasser, bei aktiver Konzentrat-GAK zusätzlich das behandelte, wieder
      zugemischte Konzentrat. Ohne Konzentrat-GAK wird das Konzentrat abgezogen.
    * Viele Parameter sind Platzhalter, die mit Pilotdaten kalibriert werden sollen.
    * Der **PFOA+PFOS-Wert** ist ein eigener Szenario-Zielwert, kein gesetzlicher Grenzwert.
    """)

# =========================================================
# 2. SEITENLEISTE
# =========================================================

st.sidebar.header("⚙️ Verfahrensauswahl")
verfahrens_wahl = st.sidebar.radio(
    "Konzept (Einzelansicht & Fließbild)",
    ["Bestandssystem (Ozonung)", "Nanofiltration (NF)"],
)
ist_ozon = (verfahrens_wahl == "Bestandssystem (Ozonung)")
st.sidebar.caption("Die Registerkarten *Vergleich* und *Sensitivität* nutzen immer alle Einstellungen unten.")

with st.sidebar.expander("🔀 Hydraulik & Mengenströme", expanded=True):
    q_feed = st.slider("Rohwasserzulauf Q (m³/h)", 10.0, 500.0, 100.0, 5.0)
    bypass = st.slider("Bypass-Anteil Bioablauf (nur NF)", 0.0, 0.80, 0.0, 0.05)

with st.sidebar.expander("🧪 Ozonung", expanded=ist_ozon):
    ozon_dosis = st.slider("Ozon-Dosis (mg O₃ / L)", 1.0, 10.0, 4.0, 0.5)
    opex_ozon = st.slider("Spezifische Betriebskosten Ozonung (€/m³, Richtwert)", 0.010, 0.050, 0.018, 0.001)

with st.sidebar.expander("🧬 Nanofiltration (Membran)", expanded=not ist_ozon):
    recovery = st.slider("Permeatausbeute Y", 0.50, 0.95, 0.80, 0.05)
    recycle = st.slider("Konzentrat-Rückführung in die Biologie (Anteil)", 0.0, 0.90, 0.20, 0.05)
    mwco = st.slider("Trenngrenze MWCO bei 90 % Rückhalt (Da)", 50, 500, 200, 10)
    zeta = st.slider("Zeta-Potential (mV)", -30.0, 0.0, -15.0, 1.0)
    sigma = st.slider("Porenweitenverteilung (σ)", 0.1, 0.6, 0.3, 0.05)
    delta_p = st.slider("Membrandruck ΔP_NF (bar)", 4.0, 25.0, 10.0, 0.5)

with st.sidebar.expander("🪨 Aktivkohle (GAK)", expanded=False):
    prod_gak = st.checkbox("Produktwasser-GAK aktivieren", value=False)
    konz_gak = st.checkbox("Konzentrat-GAK aktivieren (behandeltes Konzentrat wird wieder zugemischt)",
                           value=False)

    prod_ebct, prod_tage, prod_doc, prod_dp = 15.0, 30, 0.35, 0.6
    konz_ebct, konz_tage, konz_doc, konz_dp = 15.0, 30, 0.35, 0.6
    r_doc = 0.90

    if prod_gak:
        st.markdown("**Produktwasser-GAK**")
        prod_ebct = st.slider("EBCT Produktwasser (min)", 5.0, 30.0, 15.0, 1.0)
        prod_tage = st.slider("Standzeit bis Kohletausch Produktwasser (Tage)", 30, 360, 30, 30)
        prod_doc = st.slider("Matrix-Faktor Produktwasser (DOC-Konkurrenz)", 0.1, 0.8, 0.35, 0.05)
        prod_dp = st.slider("Druckverlust Produktwasser-GAK (bar)", 0.2, 2.0, 0.6, 0.1)
    if konz_gak:
        st.markdown("**Konzentrat-GAK**")
        konz_ebct = st.slider("EBCT Konzentrat (min)", 5.0, 30.0, 15.0, 1.0)
        konz_tage = st.slider("Standzeit bis Kohletausch Konzentrat (Tage)", 30, 360, 30, 30)
        konz_doc = st.slider("Matrix-Faktor Konzentrat (vor DOC-Aufkonzentrierung)", 0.1, 0.8, 0.35, 0.05)
        konz_dp = st.slider("Druckverlust Konzentrat-GAK (bar)", 0.2, 2.0, 0.6, 0.1)
        r_doc = st.slider("DOC-Rückhalt der NF", 0.50, 0.99, 0.90, 0.01,
                          help="Bestimmt, wie stark der DOC im Konzentrat aufkonzentriert wird.")
    if prod_gak or konz_gak:
        st.markdown("**Gemeinsame GAK-Einstellungen**")
        mittelwert = st.checkbox("Mittlere Elimination über die Standzeit (statt Momentaufnahme am Ende)",
                                 value=True)
        lastabhaengig = st.checkbox("Lastabhängige Kapazität (Freundlich, experimentell)", value=False,
                                    help="BV50 wird mit (C_Zulauf / C_NF-Zulauf)^(1/n − 1) skaliert.")
        n_inv = st.slider("Freundlich 1/n (Platzhalter)", 0.2, 0.9, 0.5, 0.05) if lastabhaengig else 0.5
    else:
        mittelwert, lastabhaengig, n_inv = True, False, 0.5

with st.sidebar.expander("🧫 Vorklärung & Biologie", expanded=False):
    vkl_option = st.selectbox("Vorklärung Verweilzeit", list(sim.VORKLAERUNG_ETA.keys()), index=1)
    hrt_h = st.number_input("Aufenthaltszeit Belebung HRT (h)", 2.0, 48.0, 12.0, 1.0)
    srt = st.number_input("Schlammalter SRT (Tage)", 0.1, 100.0, 15.0, 1.0)
    x_ts = st.number_input("Trockensubstanz X_TS (kg/m³)", 0.5, 8.0, 3.5, 0.1)
    f_s = st.number_input("Inerter gelöster CSB-Anteil f_S (DWA-A 131)", 0.0, 0.20, 0.05, 0.01)
    if srt * 24.0 <= hrt_h:
        st.warning("Schlammalter muss größer als die Aufenthaltszeit sein.")

with st.sidebar.expander("⚖️ Abwasserverordnung (AbwV)", expanded=False):
    groessenklasse = st.selectbox("Anlagengrößenklasse", list(sim.GRENZWERTE_ABWV.keys()), index=3)
grenzwert_csb = sim.GRENZWERTE_ABWV[groessenklasse]["CSB"]
grenzwert_bsb5 = sim.GRENZWERTE_ABWV[groessenklasse]["BSB5"]

st.sidebar.markdown("---")
st.sidebar.subheader("📄 Dokumentation")
try:
    with open("Begleitdokument.pdf", "rb") as pdf_file:
        st.sidebar.download_button(
            label="📥 Begleitdokument öffnen / laden",
            data=pdf_file.read(),
            file_name="Begleitdokument_Abwasser_Simulation.pdf",
            mime="application/pdf",
        )
except FileNotFoundError:
    st.sidebar.info("ℹ️ Begleitdokument derzeit nicht verfügbar.")

# =========================================================
# 3. PROZESSFLIEßBILD (SVG, schwarz-weiß)
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


fl_svg = None if ist_ozon else sim.kreislauf_fluesse(q_feed, bypass, recovery, recycle)
svg_schema = build_process_svg(
    is_ozon=ist_ozon,
    has_bypass=(not ist_ozon) and bypass > 0.0,
    aktivieren_gak=(prod_gak and not ist_ozon),
    aktivieren_konz_gak=(konz_gak and not ist_ozon),
    fl=fl_svg,
    q_feed=q_feed,
)
html_schema = f"""
<div style="background: #ffffff; padding: 16px; border-radius: 8px; border: 1px solid #000000; width: 100%; box-sizing: border-box; overflow-x: auto;">
    <div style="min-width: 820px; max-width: 1280px; margin: 0 auto;">{svg_schema}</div>
</div>
"""
components.html(html_schema, height=370, scrolling=False)
st.divider()

# =========================================================
# 4. ZULAUFKONZENTRATIONEN
# =========================================================

st.header("🧪 Rohwasser-Zulaufkonzentrationen")
if "c_roh" not in st.session_state:
    st.session_state.c_roh = dict(sim.DEFAULT_CONCENTRATIONS)

if st.button("↺ Standardwerte wiederherstellen"):
    st.session_state.c_roh = dict(sim.DEFAULT_CONCENTRATIONS)
    for key in [k for k in st.session_state.keys() if str(k).startswith("inp_")]:
        del st.session_state[key]
    st.rerun()

sum_cols = st.columns(2)
for i, stoff in enumerate(["CSB", "BSB5"]):
    with sum_cols[i]:
        st.session_state.c_roh[stoff] = st.number_input(
            f"{stoff} (mg/L)", min_value=0.0,
            value=float(st.session_state.c_roh[stoff]), step=10.0, key=f"inp_{stoff}")

trace_stoffe = [s for s in sim.DEFAULT_CONCENTRATIONS if s not in ("CSB", "BSB5")]
for start in range(0, len(trace_stoffe), 4):
    cols = st.columns(4)
    for col, stoff in zip(cols, trace_stoffe[start:start + 4]):
        with col:
            st.session_state.c_roh[stoff] = st.number_input(
                f"{stoff} (ng/L)", min_value=0.0,
                value=float(st.session_state.c_roh[stoff]), step=50.0, key=f"inp_{stoff}")

# Parameter-Dictionary für den Rechenkern
p_basis = sim.standard_parameter()
p_basis.update({
    "q_feed": q_feed, "bypass": bypass,
    "ozon_dosis": ozon_dosis, "opex_ozon": opex_ozon,
    "recovery": recovery, "recycle": recycle, "mwco": float(mwco), "zeta": zeta,
    "sigma": sigma, "delta_p": delta_p,
    "prod_gak": prod_gak, "konz_gak": konz_gak,
    "prod_ebct": prod_ebct, "prod_tage": float(prod_tage), "prod_doc": prod_doc, "prod_dp": prod_dp,
    "konz_ebct": konz_ebct, "konz_tage": float(konz_tage), "konz_doc": konz_doc, "konz_dp": konz_dp,
    "r_doc": r_doc, "mittelwert": mittelwert, "lastabhaengig": lastabhaengig, "n_inv": n_inv,
    "vkl_eta": sim.VORKLAERUNG_ETA[vkl_option],
    "hrt_h": hrt_h, "srt": srt, "x_ts": x_ts, "f_s": f_s,
    "c_roh": dict(st.session_state.c_roh),
})

st.divider()


# =========================================================
# 5. HILFSFUNKTIONEN FÜR DIE ANZEIGE
# =========================================================

def fmt_de(wert, stellen=0):
    """Zahl im deutschen Format (Tausenderpunkt, Dezimalkomma)."""
    if wert is None or (isinstance(wert, float) and not np.isfinite(wert)):
        return "–"
    text = f"{wert:,.{stellen}f}"
    return text.replace(",", "X").replace(".", ",").replace("X", ".")


def zeige_zielwerte(kenn, ist_ozon_ansicht, konz_gak_aktiv):
    pfas_ziel_ok = kenn["pfoa_pfos_end"] <= sim.SZENARIO_ZIELWERT_PFOA_PFOS
    z1, z2, z3, z4, z5 = st.columns(5)
    z1.metric(f"CSB (max. {grenzwert_csb:.0f})", f"{kenn['csb_end']:.1f} mg/L",
              delta="Eingehalten" if kenn["csb_end"] <= grenzwert_csb else "Überschritten",
              delta_color="normal" if kenn["csb_end"] <= grenzwert_csb else "inverse")
    z2.metric(f"BSB5 (max. {grenzwert_bsb5:.0f})", f"{kenn['bsb_end']:.1f} mg/L",
              delta="Eingehalten" if kenn["bsb_end"] <= grenzwert_bsb5 else "Überschritten",
              delta_color="normal" if kenn["bsb_end"] <= grenzwert_bsb5 else "inverse")
    z3.metric(f"PFOA + PFOS (Szenario-Ziel {sim.SZENARIO_ZIELWERT_PFOA_PFOS:.0f})",
              f"{kenn['pfoa_pfos_end']:.1f} ng/L",
              delta="Unter Szenario-Zielwert" if pfas_ziel_ok else "Über Szenario-Zielwert",
              delta_color="normal" if pfas_ziel_ok else "inverse",
              help="Eigener Szenario-Zielwert des Projekts, kein gesetzlicher Grenzwert.")
    z4.metric("PFAS-Summe (Endwasser)", f"{kenn['pfas_sum_end']:.1f} ng/L")
    z5.metric("Ø Arzneimittel-Red. (Endwasser)", f"{kenn['avg_phac_red']:.1f} %")

    if (not ist_ozon_ansicht) and konz_gak_aktiv:
        st.markdown("**Nur Produktwasser (ohne wieder zugemischtes Konzentrat):**")
        y1, y2, y3 = st.columns(3)
        y1.metric("PFOA + PFOS (Produktwasser)", f"{kenn['pfoa_pfos_prod']:.1f} ng/L")
        y2.metric("PFAS-Summe (Produktwasser)", f"{kenn['pfas_sum_prod']:.1f} ng/L")
        y3.metric("Ø Arzneimittel-Red. (Produktwasser)", f"{kenn['avg_phac_red_prod']:.1f} %")


def plot_konzentrationsverlauf(df_res, substanz, ist_ozon_ansicht):
    row = df_res[df_res["Substanzname"] == substanz].iloc[0]
    if ist_ozon_ansicht:
        labels = ["Rohwasser", "Bioablauf", "Endwasser"]
        werte = [row["C_Rohwasser"], row["C_Bio_Ablauf"], row["C_Endwasser"]]
    else:
        labels = ["Rohwasser", "Bioablauf", "Produktwasser (roh)", "Produktwasser (nach GAK)", "Endwasser"]
        werte = [row["C_Rohwasser"], row["C_Bio_Ablauf"], row["C_Produktwasser_roh"],
                 row["C_Produktwasser"], row["C_Endwasser"]]
    werte = [max(float(w), 1e-3) for w in werte]
    fig = go.Figure(go.Scatter(
        x=labels, y=werte, mode="lines+markers+text",
        text=[f"{w:,.1f}".replace(",", ".") for w in werte], textposition="top center",
        line=dict(color="black", width=2), marker=dict(size=9, color="white", line=dict(color="black", width=2)),
    ))
    fig.update_layout(
        title=f"Konzentrationsverlauf <b>{substanz}</b> ({row['Einheit']})",
        template="plotly_white", height=400, showlegend=False,
        yaxis=dict(type="log", title=f"Konzentration ({row['Einheit']}, log)", showline=True, linecolor="black"),
        xaxis=dict(showline=True, linecolor="black"),
        margin=dict(l=50, r=25, t=60, b=50), font=dict(family="Arial", color="black"),
    )
    return fig


def plot_frachtverteilung(df_res):
    d = df_res[df_res["Gruppe"] != "Summenparameter"].copy()
    
    spalten = [
        ("Fracht_Biologie_%", "Biologie (Abbau + Sorption)", "#2a7f7a"),  
        ("Fracht_Prod_GAK_%", "Produktwasser-GAK", "#5bb5ad"),           
        ("Fracht_Konz_GAK_%", "Konzentrat-GAK", "#a5d9d3"),           
        ("Fracht_Konz_Abzug_%", "Konzentrat-Abzug", "#b0b0b0"),       
        ("Fracht_Endwasser_%", "Verbleib im Endwasser", "#e4572e"),   
    ]
    
    if "Fracht_Endwasser_%" in d.columns:
        d = d.sort_values(by="Fracht_Endwasser_%", ascending=True)

    fig = go.Figure()
    
    for col, name, farbe in spalten:
        if col not in d.columns or d[col].fillna(0).abs().sum() <= 1e-9:
            continue
        
        # Textwerte nur anzeigen, wenn der Anteil > 1% ist
        text_werte = [f"{val:.1f}%" if val > 1.0 else "" for val in d[col]]
        
        fig.add_trace(go.Bar(
            y=d["Substanzname"], 
            x=d[col], 
            name=name, 
            orientation="h",
            text=text_werte,
            textposition="inside",
            insidetextanchor="middle",
            marker=dict(
                color=farbe, 
                line=dict(color="black", width=0.8)
            ),
            hovertemplate="<b>%{y}</b><br>%{name}: %{x:.1f}%<extra></extra>"
        ))

    fig.update_layout(
        barmode="stack", 
        title="<b>Frachtverteilung</b> (Anteil an der Fracht nach Vorklärung)",
        template="plotly_white", 
        height=450,
        xaxis=dict(
            title="Anteil an der Fracht [%]", 
            range=[0, 100], 
            showline=True, 
            linecolor="black",
            gridcolor="#eeeeee"
        ),
        yaxis=dict(
            showline=True, 
            linecolor="black"
        ), 
        margin=dict(l=80, r=25, t=60, b=80),
        legend=dict(
            orientation="h", 
            y=-0.25, 
            x=0.5, 
            xanchor="center"
        ), 
        font=dict(family="Arial", color="black"),
    )
    return fig

# =========================================================
# 6. REGISTERKARTEN
# =========================================================

tab_ein, tab_vgl = st.tabs(["📌 Einzelkonzept", "⚖️ Vergleich Ozonung ↔ NF"])

# ---------------------------------------------------------
# 6.1 EINZELKONZEPT
# ---------------------------------------------------------
with tab_ein:
    p_ein = dict(p_basis)
    p_ein["konzept"] = "ozon" if ist_ozon else "nf"
    res = sim.simuliere(p_ein)
    df_res, fl, kenn = res["df"], res["fl"], res["kenn"]

    label = ("Bestandssystem (Ozonung)" if ist_ozon
             else f"Nanofiltration (Bypass: {int(round(bypass * 100))} %)")
    st.header(f"🎯 Zielwertüberwachung: {label}")
    zeige_zielwerte(kenn, ist_ozon, konz_gak)

    if not ist_ozon:
        if kenn["bilanz_abw"] > 1e-6:
            st.error(f"Massenbilanz nicht geschlossen (max. Abweichung {kenn['bilanz_abw']:.3g} %)")
        else:
            st.caption("✅ Frachtbilanz geschlossen (Biologie + GAK-Stufen + Konzentrat-Abzug + Endwasser = 100 %).")
        if not konz_gak:
            st.caption("ℹ️ Ohne Konzentrat-GAK verlässt das Konzentrat die Anlage als Abzug; "
                       "Endwasser = Produktwasser.")
        if konz_gak:
            st.caption(f"ℹ️ DOC im Konzentrat {kenn['aufkonz_doc']:.1f}-fach aufkonzentriert → "
                       f"wirksamer Matrix-Faktor Konzentrat-GAK: {kenn['konz_doc_eff']:.3f}.")

    st.divider()

    # Mengenströme
    if not ist_ozon:
        st.header("🔁 Mengenströme")
        h1, h2, h3, h4, h5 = st.columns(5)
        h1.metric("Biologie-Zulauf", f"{fl['q_bio']:.1f} m³/h",
                  delta=f"+{(fl['q_bio'] / q_feed - 1.0) * 100.0:.1f} % hydraulisch", delta_color="off")
        h2.metric("NF-Zulauf", f"{fl['q_nf']:.1f} m³/h")
        h3.metric("Produktwasser", f"{fl['q_prod']:.1f} m³/h")
        h4.metric("Rückführung", f"{fl['q_rez']:.1f} m³/h")
        h5.metric("Konzentrat (Abzug bzw. zur GAK)", f"{fl['q_konz_aus']:.1f} m³/h")
        st.divider()

    # Energie & Betriebsmittel
    st.header("⚡ Energie, Betriebsmittel & Kosten")
    ref = sim.ENERGIE_OZONUNG_AACHEN
    if ist_ozon:
        e1, e2, e3 = st.columns(3)
        e1.metric("Energiebedarf Ozonung (KA Aachen-Sörs 2018)", f"{ref:.3f} kWh/m³",
                  help="Spezifischer Energiebedarf der Ozonung der KA Aachen-Sörs (2018), je m³ Zulauf.")
        e2.metric("Jahresenergiebedarf Ozonung", f"{fmt_de(kenn['e_jahr_kwh'])} kWh/a")
        e3.metric("Jährliche Betriebskosten (pauschal)", f"{fmt_de(kenn['opex_ozon_eur_a'])} €/a")
    else:
        e1, e2, e3, e4 = st.columns(4)
        e1.metric("Spez. Energie NF", f"{kenn['sec_nf']:.2f} kWh/m³ Permeat")
        e2.metric("Elektrische Leistung gesamt", f"{kenn['p_ges_kw']:.1f} kW")
        e3.metric("Jahresenergie gesamt", f"{fmt_de(kenn['e_jahr_kwh'])} kWh/a")
        e4.metric("Spez. Energie je m³ Zulauf", f"{kenn['spez_energie_kwh_m3']:.3f} kWh/m³",
                  delta=f"{kenn['spez_energie_kwh_m3'] / ref:.1f}× Ozonung", delta_color="off",
                  help=f"Vergleichswert Ozonung (KA Aachen-Sörs 2018): {ref:.3f} kWh/m³. "
                       "Die NF-Werte enthalten nur Hochdruckpumpe und GAK-Druckverlust (Untergrenze).")
        if prod_gak or konz_gak:
            k1, k2, k3, k4 = st.columns(4)
            k1.metric("Kohleverbrauch Produktwasser", f"{kenn['kohle_prod_g_m3']:.1f} g/m³")
            k2.metric("Kohleverbrauch Konzentrat", f"{kenn['kohle_konz_g_m3']:.1f} g/m³")
            k3.metric("Kohlebedarf pro Jahr", f"{kenn['kohle_gesamt_t_a']:.1f} t/a")
            k4.metric("Leistung GAK-Stufen", f"{kenn['p_gak_prod_kw'] + kenn['p_gak_konz_kw']:.1f} kW")
    st.divider()

    # Tabelle
    st.header("📊 Detaillierte Stoffbilanz")
    spalten = ["Substanzname", "Gruppe", "Einheit", "C_Rohwasser", "Bio_Elimination_%"]
    if ist_ozon:
        spalten.append("Ozon_Elimination_%")
    else:
        spalten.append("R_NF_%")
        if prod_gak:
            spalten.append("GAK_Elimination_%")
        if konz_gak:
            spalten.append("Konz_GAK_Elimination_%")
        spalten.append("C_Produktwasser")
    spalten += ["C_Endwasser", "Konzentrationsreduktion_%"]
    anzeige = df_res[spalten].copy()
    for c in anzeige.columns:
        if c not in ("Substanzname", "Gruppe", "Einheit"):
            anzeige[c] = anzeige[c].round(2)
    st.dataframe(anzeige, use_container_width=True, hide_index=True)

    # Diagramme
    st.subheader("📈 Konzentrationsverlauf und Frachtverteilung")
    stoffe_liste = df_res.loc[df_res["Gruppe"] != "Summenparameter", "Substanzname"].tolist()
    gewaehlt = st.selectbox("Substanz für das Verlaufsdiagramm", stoffe_liste,
                            index=stoffe_liste.index("PFOA") if "PFOA" in stoffe_liste else 0,
                            key="sel_verlauf")
    st.plotly_chart(plot_konzentrationsverlauf(df_res, gewaehlt, ist_ozon), use_container_width=True)
    if not ist_ozon:
        st.plotly_chart(plot_frachtverteilung(df_res), use_container_width=True)

    csv = df_res.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig")
    st.download_button("Ergebnisse als CSV herunterladen", data=csv,
                       file_name="abwasser_konzept_simulation.csv", mime="text/csv")

# ---------------------------------------------------------
# 6.2 VERGLEICH
# ---------------------------------------------------------
with tab_vgl:
    st.header("⚖️ Vergleich Ozonung ↔ Nanofiltration")
    st.caption("Beide Konzepte werden mit denselben Eingaben (Zulauf, Biologie) und den jeweiligen "
               "Verfahrensparametern aus der Seitenleiste gerechnet.")

    p_o = dict(p_basis)
    p_o["konzept"] = "ozon"
    p_n = dict(p_basis)
    p_n["konzept"] = "nf"
    r_o, r_n = sim.simuliere(p_o), sim.simuliere(p_n)
    ko, kn = r_o["kenn"], r_n["kenn"]

    zeilen = [
        ("Ø Arzneimittel-Reduktion (Endwasser)", "%", "avg_phac_red", 1),
        ("Ø PFAS-Reduktion (Endwasser)", "%", "avg_pfas_red", 1),
        ("PFAS-Summe (Endwasser)", "ng/L", "pfas_sum_end", 0),
        ("PFOA + PFOS (Endwasser)", "ng/L", "pfoa_pfos_end", 1),
        ("PFAS-Summe (Produktwasser)", "ng/L", "pfas_sum_prod", 0),
        ("CSB (Endwasser)", "mg/L", "csb_end", 1),
        ("BSB5 (Endwasser)", "mg/L", "bsb_end", 1),
        ("Spez. Energiebedarf je m³ Zulauf", "kWh/m³", "spez_energie_kwh_m3", 3),
        ("Jahresenergie", "kWh/a", "e_jahr_kwh", 0),
        ("Konzentrat (Abzug bzw. zur GAK)", "m³/h", "q_konz_aus", 1),
        ("Kohlebedarf", "t/a", "kohle_gesamt_t_a", 1),
        ("Betriebskosten Ozonung (pauschal)", "€/a", "opex_ozon_eur_a", 0),
    ]
    tabelle = pd.DataFrame({
        "Kennzahl": [z[0] for z in zeilen],
        "Einheit": [z[1] for z in zeilen],
        "Ozonung": [fmt_de(ko[z[2]], z[3]) for z in zeilen],
        "Nanofiltration": [fmt_de(kn[z[2]], z[3]) for z in zeilen],
    })
    st.dataframe(tabelle, use_container_width=True, hide_index=True)
    st.caption("Der Energiewert der Ozonung ist der Referenzwert der KA Aachen-Sörs (2018). Die NF-Werte sind "
               "eine Untergrenze (nur Hochdruckpumpe und GAK-Druckverlust).")

    d_o = r_o["df"].set_index("Substanzname")
    d_n = r_n["df"].set_index("Substanzname")
    stoffe = [s for s in d_o.index if d_o.loc[s, "Gruppe"] != "Summenparameter"]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=stoffe, y=[d_o.loc[s, "Konzentrationsreduktion_%"] for s in stoffe],
                         name="Ozonung", marker=dict(color="white", line=dict(color="black", width=1.5))))
    fig.add_trace(go.Bar(x=stoffe, y=[d_n.loc[s, "Konzentrationsreduktion_%"] for s in stoffe],
                         name="Nanofiltration",
                         marker=dict(color="#888888", line=dict(color="black", width=1.5))))
    fig.update_layout(
        barmode="group", template="plotly_white", height=420,
        title="Reduktion gegenüber dem Rohwasser (Endwasser)",
        yaxis=dict(title="Reduktion [%]", range=[0, 105], showline=True, linecolor="black"),
        xaxis=dict(showline=True, linecolor="black"), legend=dict(orientation="h", y=-0.2),
        margin=dict(l=50, r=25, t=60, b=60), font=dict(family="Arial", color="black"),
    )
    st.plotly_chart(fig, use_container_width=True)


    # =====================================================
    # KOSTENVERGLEICH NF UND OZONUNG
    # =====================================================

    st.divider()
    st.subheader("💶 Spezifische Behandlungskosten")

    # Einzelkosten NF aus dem Rechenkern
    kosten_nf = [
        ("Strom", kn["kosten_strom"]),
        ("Membran (annualisiert)", kn["kosten_membran"]),
        ("Antiscalant", kn["kosten_antiscalant"]),
        ("Reinigungschemikalien", kn["kosten_cip"]),
        ("Wartung NF", kn["kosten_wartung"]),
        ("Personal NF", kn["kosten_personal_nf"]),
        ("Aktivkohle", kn["kosten_aktivkohle"]),
        ("GAK-Betrieb", kn["kosten_gak_betrieb"]),
        ("Personal GAK", kn["kosten_personal_gak"]),
        ("Konzentratentsorgung", kn["kosten_konzentrat"]),
    ]

    kosten_nf_gesamt = kn["kosten_nf_gesamt_m3"]
    kosten_ozon_gesamt = ko["kosten_ozon_gesamt_m3"]

    # Kennzahlen
    c1, c2, c3 = st.columns(3)

    c1.metric(
        "Gesamtkosten Ozonung",
        f"{fmt_de(kosten_ozon_gesamt, 4)} €/m³"
    )

    c2.metric(
        "Gesamtkosten NF",
        f"{fmt_de(kosten_nf_gesamt, 4)} €/m³"
    )

    c3.metric(
        "Differenz NF - Ozonung",
        f"{fmt_de(kosten_nf_gesamt - kosten_ozon_gesamt, 4)} €/m³",
        delta=(
            f"{(kosten_nf_gesamt / kosten_ozon_gesamt - 1) * 100:.1f} %"
            if kosten_ozon_gesamt > 0 else "Nicht berechenbar"
        ),
        delta_color="inverse"
    )

    # Daten für gestapeltes Balkendiagramm
    fig_kosten = go.Figure()

    farben = [
        "#4C78A8", "#72A0C1", "#F2CF5B", "#ECA82C",
        "#B279A2", "#79706E", "#59A14F", "#8CD17D",
        "#BAB0AC", "#E45756"
    ]

    for i, (name, wert) in enumerate(kosten_nf):
        if wert > 0:
            fig_kosten.add_trace(go.Bar(
                name=name,
                x=["Nanofiltration"],
                y=[wert],
                marker=dict(
                    color=farben[i],
                    line=dict(color="black", width=0.6)
                ),
                hovertemplate=(
                    f"<b>{name}</b><br>"
                    "%{y:.5f} €/m³<extra></extra>"
                )
            ))

    # Ozonung als pauschaler Gesamt-OPEX
    fig_kosten.add_trace(go.Bar(
        name="Gesamt-OPEX Ozonung",
        x=["Ozonung"],
        y=[kosten_ozon_gesamt],
        marker=dict(
            color="white",
            line=dict(color="black", width=1.5)
        ),
        hovertemplate=(
            "<b>Gesamt-OPEX Ozonung</b><br>"
            "%{y:.5f} €/m³<extra></extra>"
        )
    ))

    fig_kosten.update_layout(
        barmode="stack",
        template="plotly_white",
        height=500,
        title="Spezifische Kosten von Nanofiltration und Ozonung",
        yaxis=dict(
            title="Behandlungskosten [€/m³ Zulauf]",
            rangemode="tozero",
            showline=True,
            linecolor="black",
            gridcolor="#eeeeee"
        ),
        xaxis=dict(
            showline=True,
            linecolor="black"
        ),
        legend=dict(
            orientation="h",
            y=-0.25,
            x=0.5,
            xanchor="center"
        ),
        margin=dict(l=70, r=30, t=70, b=110),
        font=dict(family="Arial", color="black")
    )

    # Gesamtkosten über den Balken anzeigen
    fig_kosten.add_annotation(
        x="Nanofiltration",
        y=kosten_nf_gesamt,
        text=f"{kosten_nf_gesamt:.4f} €/m³",
        showarrow=False,
        yshift=12,
        font=dict(size=13, color="black")
    )

    fig_kosten.add_annotation(
        x="Ozonung",
        y=kosten_ozon_gesamt,
        text=f"{kosten_ozon_gesamt:.4f} €/m³",
        showarrow=False,
        yshift=12,
        font=dict(size=13, color="black")
    )

    st.plotly_chart(fig_kosten, use_container_width=True)

    st.caption(
        "Bezugsgröße: m³ behandelter Zulauf. "
        "Die Ozonung wird mit einem pauschalen Gesamt-OPEX von 1,80 ct/m³ "
        "angesetzt. Die NF-Kosten enthalten annualisierte Membrankosten "
        "und die berechneten Betriebskosten. "
        "Alle Kostenbestandteile der NF werden einzeln dargestellt."
    )


