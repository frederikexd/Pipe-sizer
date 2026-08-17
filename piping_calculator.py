"""
piping_calculator.py - Pipe Sizing & Pressure Drop Calculator

Darcy-Weisbach pressure drop with iterative Colebrook-White friction factor.
Pipe dimensions per ASME B36.10M / AS 1074. DN and NPS sizing.
Standard fitting K-factors (Crane TP-410). Elevation head.

Applicable to AU / US / UK / EU jurisdictions.

References:
  - AS 4041: Pressure piping (Australia)
  - AS 1074: Steel tubes and tubulars (Australia)
  - AS/NZS 1477: PVC pipes (Australia/NZ)
  - AS 1432: Copper tubes (Australia)
  - AS 3500: Plumbing and drainage (velocity guidelines)
  - AS 2419: Fire hydrant installations (velocity guidelines)
  - ASME B36.10M: Welded and seamless wrought steel pipe (US)
  - Crane TP-410: Flow of fluids through valves, fittings, and pipe

Run standalone:
    pip install streamlit pandas
    streamlit run piping_calculator.py

--------------------------------------------------------------------
Pipe Sizing & Pressure Drop Calculator - Darcy-Weisbach flow analysis.
Copyright (C) 2026  Frederik Thio

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import streamlit as st
import math
import pandas as pd
from typing import Tuple

PIPE_DATA = {
    "DN15 (1/2\")":    {"dn": 15,  "nps": "1/2\"",   "od": 21.3,  "40": 15.8,  "80": 13.9},
    "DN20 (3/4\")":    {"dn": 20,  "nps": "3/4\"",   "od": 26.7,  "40": 20.9,  "80": 18.9},
    "DN25 (1\")":      {"dn": 25,  "nps": "1\"",     "od": 33.4,  "40": 26.6,  "80": 24.3},
    "DN32 (1-1/4\")":  {"dn": 32,  "nps": "1-1/4\"", "od": 42.2,  "40": 35.1,  "80": 32.5},
    "DN40 (1-1/2\")":  {"dn": 40,  "nps": "1-1/2\"", "od": 48.3,  "40": 40.9,  "80": 38.1},
    "DN50 (2\")":      {"dn": 50,  "nps": "2\"",     "od": 60.3,  "40": 52.5,  "80": 49.3},
    "DN65 (2-1/2\")":  {"dn": 65,  "nps": "2-1/2\"", "od": 73.0,  "40": 62.7,  "80": 59.0},
    "DN80 (3\")":      {"dn": 80,  "nps": "3\"",     "od": 88.9,  "40": 77.9,  "80": 73.7},
    "DN100 (4\")":     {"dn": 100, "nps": "4\"",     "od": 114.3, "40": 102.3, "80": 97.2},
    "DN150 (6\")":     {"dn": 150, "nps": "6\"",     "od": 168.3, "40": 154.1, "80": 146.3},
    "DN200 (8\")":     {"dn": 200, "nps": "8\"",     "od": 219.1, "40": 202.7, "80": 193.7},
    "DN250 (10\")":    {"dn": 250, "nps": "10\"",    "od": 273.1, "40": 254.5, "80": 242.9},
    "DN300 (12\")":    {"dn": 300, "nps": "12\"",    "od": 323.8, "40": 303.2, "80": 288.8},
}

SCHEDULES = ["40", "80"]

MATERIALS = {
    "Carbon Steel (AS 1074 / ASME B36.10M)":  0.046,
    "Stainless Steel (AS 1528)":               0.015,
    "Copper (AS 1432)":                        0.0015,
    "PVC (AS/NZS 1477)":                       0.0015,
    "PE/HDPE (AS/NZS 4130)":                   0.007,
    "Galvanized Steel (AS 1074)":              0.15,
    "Cast Iron":                                0.26,
    "Ductile Iron (AS/NZS 2280)":              0.03,
    "Concrete Lined":                           0.3,
    "GRP/FRP":                                  0.01,
}

FLUIDS = {
    "Water (20C)":               {"density": 998,   "viscosity": 0.001002},
    "Water (40C)":               {"density": 992,   "viscosity": 0.000653},
    "Water (60C)":               {"density": 983,   "viscosity": 0.000467},
    "Water (80C)":               {"density": 972,   "viscosity": 0.000355},
    "Air (20C, 1 atm)":         {"density": 1.204, "viscosity": 0.0000181},
    "Air (40C, 1 atm)":         {"density": 1.127, "viscosity": 0.0000190},
    "Hydraulic Oil (ISO 32)":    {"density": 860,   "viscosity": 0.032},
    "Hydraulic Oil (ISO 46)":    {"density": 870,   "viscosity": 0.046},
    "Hydraulic Oil (ISO 68)":    {"density": 880,   "viscosity": 0.068},
    "Glycol 30% (20C)":         {"density": 1040,  "viscosity": 0.0024},
    "Glycol 50% (20C)":         {"density": 1070,  "viscosity": 0.0038},
    "Diesel Fuel":               {"density": 850,   "viscosity": 0.0035},
    "Mine Process Water (~25C)": {"density": 1020,  "viscosity": 0.001},
    "Seawater (20C)":            {"density": 1025,  "viscosity": 0.00108},
}

FLOW_TO_M3H = {
    "L/s":      3.6,
    "m3/h":     1.0,
    "L/min":    0.06,
    "GPM (US)": 0.2271,
}

FITTINGS = {
    "90 Elbow (standard)":    0.9,
    "90 Elbow (long radius)": 0.6,
    "45 Elbow":               0.4,
    "Tee (thru-run)":         0.4,
    "Tee (branch)":           1.8,
    "Gate Valve (open)":      0.2,
    "Globe Valve (open)":     10.0,
    "Ball Valve (open)":      0.05,
    "Check Valve (swing)":    2.5,
    "Butterfly Valve (open)": 0.35,
    "Reducer (sudden)":       0.5,
    "Expansion (sudden)":     1.0,
    "Entry (sharp-edged)":    0.5,
    "Entry (bell-mouth)":     0.05,
    "Exit (submerged)":       1.0,
    "Strainer (Y-type)":      2.0,
    "Flexible coupling":      0.6,
}

VELOCITY_LIMITS = {
    "Water - suction (AS 3500)":       (0.5, 2.0),
    "Water - discharge":               (1.0, 3.5),
    "Water - gravity/drain":           (0.3, 1.5),
    "Water - fire services (AS 2419)": (1.0, 4.0),
    "Hydraulic Oil - pressure":        (2.0, 6.0),
    "Hydraulic Oil - suction":         (0.5, 1.5),
    "Hydraulic Oil - return":          (1.0, 4.0),
    "Air - low pressure":             (5.0, 25.0),
    "Air - compressed":               (6.0, 15.0),
    "Steam - low pressure":           (15.0, 25.0),
    "Steam - high pressure":          (25.0, 50.0),
    "Mining slurry":                   (1.5, 5.0),
}


def colebrook_white(Re, roughness_m, diameter_m):
    if Re <= 0:
        return 0.0
    if Re < 2300:
        return 64.0 / Re
    eD = roughness_m / diameter_m
    f = 0.02
    for _ in range(100):
        rhs = -2.0 * math.log10(eD / 3.7 + 2.51 / (Re * math.sqrt(f)))
        f_new = 1.0 / (rhs * rhs)
        if abs(f_new - f) < 1e-10:
            break
        f = f_new
    return f


def calc_segment(length, elevation, flow_m3h, fittings_k, pipe_id_m, roughness_m, density, viscosity):
    area = math.pi * (pipe_id_m / 2) ** 2
    Q = flow_m3h / 3600
    V = Q / area if area > 0 else 0
    Re = (density * V * pipe_id_m) / viscosity if viscosity > 0 else 0
    f = colebrook_white(Re, roughness_m, pipe_id_m)
    dp_friction = f * (length / pipe_id_m) * (density * V * V / 2) if pipe_id_m > 0 else 0
    dp_fittings = fittings_k * (density * V * V / 2)
    dp_elevation = density * 9.81 * elevation
    dp_total = dp_friction + dp_fittings + dp_elevation
    if Re < 2300:
        regime = "Laminar"
    elif Re < 4000:
        regime = "Transitional"
    else:
        regime = "Turbulent"
    return {
        "velocity": V, "reynolds": Re, "friction_factor": f, "regime": regime,
        "dp_friction_kpa": dp_friction / 1000, "dp_fittings_kpa": dp_fittings / 1000,
        "dp_elevation_kpa": dp_elevation / 1000, "dp_total_kpa": dp_total / 1000,
        "dp_total_pa": dp_total, "total_k": fittings_k,
    }


def velocity_status(v, limits):
    v_min, v_max = limits
    if v < v_min:
        return "LOW - risk of sedimentation", "warning"
    elif v > v_max:
        return "HIGH - risk of erosion/noise/water hammer", "error"
    return "OK", "ok"


def render_piping_tab():
    st.markdown("### Pipe Sizing & Pressure Drop Calculator")
    st.markdown(
        "Darcy-Weisbach with Colebrook-White friction factor. "
        "Pipe dimensions per ASME B36.10M / AS 1074. "
        "DN and NPS sizing. Applicable to AU, US, UK, and EU standards."
    )

    if "pipe_segments" not in st.session_state:
        st.session_state.pipe_segments = [
            {"name": "Supply line", "length": 15.0, "elevation": 0.0, "flow": 2.0,
             "fittings": {"Entry (sharp-edged)": 1, "90 Elbow (standard)": 3, "Gate Valve (open)": 1}},
            {"name": "Riser", "length": 8.0, "elevation": 8.0, "flow": 2.0,
             "fittings": {"90 Elbow (long radius)": 2}},
            {"name": "Branch to equipment", "length": 6.0, "elevation": 0.0, "flow": 2.0,
             "fittings": {"Tee (branch)": 1, "Ball Valve (open)": 1, "Exit (submerged)": 1}},
        ]

    st.markdown("#### System Configuration")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        pipe_size = st.selectbox("Pipe size", list(PIPE_DATA.keys()), index=7)
    with c2:
        schedule = st.selectbox("Schedule", SCHEDULES, index=0)
    with c3:
        material = st.selectbox("Material", list(MATERIALS.keys()), index=0)
    with c4:
        fluid = st.selectbox("Fluid", list(FLUIDS.keys()), index=0)

    pipe_info = PIPE_DATA[pipe_size]
    pipe_id_mm = pipe_info.get(schedule, 52.5)
    pipe_id_m = pipe_id_mm / 1000
    roughness_mm = MATERIALS[material]
    roughness_m = roughness_mm / 1000
    fluid_props = FLUIDS[fluid]

    fc1, _ = st.columns([1, 3])
    with fc1:
        flow_unit = st.selectbox("Flow unit", list(FLOW_TO_M3H.keys()), index=0)

    ic1, ic2, ic3, ic4 = st.columns(4)
    ic1.metric("Internal Diameter", f"{pipe_id_mm:.1f} mm")
    ic2.metric("Roughness", f"{roughness_mm} mm")
    ic3.metric("Density", f"{fluid_props['density']} kg/m3")
    ic4.metric("Viscosity", f"{fluid_props['viscosity']} Pa.s")

    st.markdown("#### Pipe Segments")
    segments = st.session_state.pipe_segments

    for i, seg in enumerate(segments):
        with st.expander(f"Segment: {seg['name']}", expanded=True):
            sc1, sc2, sc3, sc4 = st.columns([3, 1.5, 1.5, 1.5])
            with sc1:
                seg["name"] = st.text_input("Name", value=seg["name"], key=f"sn_{i}")
            with sc2:
                seg["length"] = st.number_input("Length (m)", value=seg["length"], step=0.5, min_value=0.0, key=f"sl_{i}")
            with sc3:
                seg["elevation"] = st.number_input("Elev. (m)", value=seg["elevation"], step=0.5, key=f"se_{i}")
            with sc4:
                seg["flow"] = st.number_input(f"Flow ({flow_unit})", value=seg["flow"], step=0.1, min_value=0.001, key=f"sf_{i}")

            st.markdown("**Fittings**")
            for ft_name in list(seg["fittings"].keys()):
                fc1, fc2, fc3 = st.columns([4, 1, 0.5])
                with fc1:
                    st.text(f"{ft_name}  (K={FITTINGS[ft_name]})")
                with fc2:
                    nq = st.number_input("Qty", value=seg["fittings"][ft_name], min_value=0, step=1, key=f"fq_{i}_{ft_name}", label_visibility="collapsed")
                    if nq == 0:
                        del seg["fittings"][ft_name]
                    else:
                        seg["fittings"][ft_name] = nq
                with fc3:
                    if st.button("x", key=f"fd_{i}_{ft_name}"):
                        del seg["fittings"][ft_name]
                        st.rerun()

            avail = [f for f in FITTINGS.keys() if f not in seg["fittings"]]
            if avail:
                af1, af2 = st.columns([3, 1])
                with af1:
                    nf = st.selectbox("Add", avail, key=f"fn_{i}", label_visibility="collapsed")
                with af2:
                    if st.button("+ Add", key=f"fa_{i}"):
                        seg["fittings"][nf] = 1
                        st.rerun()

            if len(segments) > 1 and st.button("Remove segment", key=f"sr_{i}"):
                segments.pop(i)
                st.rerun()

    if st.button("+ Add segment"):
        segments.append({"name": f"Segment {len(segments)+1}", "length": 5.0, "elevation": 0.0, "flow": segments[0]["flow"] if segments else 2.0, "fittings": {}})
        st.rerun()

    st.markdown("#### Velocity Reference")
    vc1, vc2 = st.columns(2)
    with vc1:
        vel_ref = st.selectbox("Application", list(VELOCITY_LIMITS.keys()), index=1)
    vel_limits = VELOCITY_LIMITS[vel_ref]
    with vc2:
        st.markdown(f"**Recommended:** {vel_limits[0]} - {vel_limits[1]} m/s")

    if st.button("Calculate Pressure Drop", type="primary"):
        conv = FLOW_TO_M3H[flow_unit]
        results = []
        total_dp = 0
        total_length = 0
        total_elev = 0

        for seg in segments:
            total_k = sum(FITTINGS.get(ft, 0) * qty for ft, qty in seg["fittings"].items())
            flow_m3h = seg["flow"] * conv
            res = calc_segment(seg["length"], seg["elevation"], flow_m3h, total_k, pipe_id_m, roughness_m, fluid_props["density"], fluid_props["viscosity"])
            res["name"] = seg["name"]
            res["flow_m3h"] = flow_m3h
            results.append(res)
            total_dp += res["dp_total_pa"]
            total_length += seg["length"]
            total_elev += seg["elevation"]

        st.markdown("---")
        st.markdown("### Results")

        tc1, tc2, tc3 = st.columns(3)
        tc1.metric("Total Pressure Drop", f"{total_dp/1000:.2f} kPa", delta=f"{total_dp/100000:.3f} bar | {total_dp*0.000145038:.2f} psi", delta_color="off")
        head = total_dp / (fluid_props["density"] * 9.81) if fluid_props["density"] > 0 else 0
        tc2.metric("Total Head Loss", f"{head:.2f} m")
        tc3.metric("System", f"{total_length:.1f} m | {total_elev:.1f} m elev.")

        st.markdown("#### Segment Breakdown")
        for res in results:
            vs, vt = velocity_status(res["velocity"], vel_limits)
            icon = {"ok": "✅", "warning": "⚠️", "error": "🔴"}[vt]
            st.markdown(f"**{res['name']}** - {icon} {res['velocity']:.2f} m/s ({vs})")
            r1, r2, r3, r4 = st.columns(4)
            r1.metric("Reynolds", f"{res['reynolds']:.0f}")
            r2.metric("Regime", res["regime"])
            r3.metric("Friction f", f"{res['friction_factor']:.5f}")
            r4.metric("Sum K", f"{res['total_k']:.2f}")
            d1, d2, d3, d4 = st.columns(4)
            d1.metric("dP friction", f"{res['dp_friction_kpa']:.2f} kPa")
            d2.metric("dP fittings", f"{res['dp_fittings_kpa']:.2f} kPa")
            d3.metric("dP elevation", f"{res['dp_elevation_kpa']:.2f} kPa")
            d4.metric("dP total", f"{res['dp_total_kpa']:.2f} kPa")
            st.markdown("---")

        warns = []
        for res in results:
            vs, vt = velocity_status(res["velocity"], vel_limits)
            if vt != "ok":
                warns.append(f"**{res['name']}**: {res['velocity']:.2f} m/s - {vs}")
        if warns:
            st.warning("**Velocity warnings:**\n\n" + "\n\n".join(warns))

        rows = []
        for res in results:
            vs, _ = velocity_status(res["velocity"], vel_limits)
            rows.append({"Segment": res["name"], "Flow (m3/h)": round(res["flow_m3h"], 3), "V (m/s)": round(res["velocity"], 3), "Status": vs.split("-")[0].strip(), "Re": int(round(res["reynolds"])), "Regime": res["regime"], "f": round(res["friction_factor"], 6), "Sum K": round(res["total_k"], 2), "dP Fric (kPa)": round(res["dp_friction_kpa"], 3), "dP Fit (kPa)": round(res["dp_fittings_kpa"], 3), "dP Elev (kPa)": round(res["dp_elevation_kpa"], 3), "dP Tot (kPa)": round(res["dp_total_kpa"], 3)})
        rows.append({"Segment": "TOTAL", "dP Tot (kPa)": round(total_dp / 1000, 3)})
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button("Download CSV", df.to_csv(index=False), "piping_results.csv", "text/csv")

        st.info(
            "**Method & Standards**\n\n"
            "Darcy-Weisbach (dP = f.L/D.rhoV2/2). Friction factor: Colebrook-White (turbulent) or 64/Re (laminar). "
            "Fitting losses: K-factor method per Crane TP-410. Elevation head: dP = rho.g.dh.\n\n"
            "Pipe dimensions per ASME B36.10M / AS 1074. Velocity guidelines per AS 3500, AS 2419, and industry practice. "
            "Material roughness from published references.\n\n"
            "**Applicable standards:** AS 4041 (AU), ASME B31 (US), BS EN 13480 (UK/EU)."
        )


if __name__ == "__main__":
    st.set_page_config(page_title="Pipe Sizing & Pressure Drop", page_icon="wrench", layout="wide")
    st.title("Pipe Sizing & Pressure Drop")
    st.caption("DN/NPS sizing | Darcy-Weisbach | AS 4041 / ASME B31 applicable")
    render_piping_tab()
