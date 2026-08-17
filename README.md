# Pipe Sizing & Pressure Drop Calculator

Darcy-Weisbach pressure drop with an iterative Colebrook-White friction
factor. Real pipe schedules, standard fitting K-factors, elevation head, and
velocity checking against application limits.

Frederik Thio.

## Run

```bash
pip install -r requirements.txt
streamlit run piping_calculator.py
```

## What it does

Build a system segment by segment — length, elevation change, flow rate, and
the fittings on that run. The tool reports, per segment and for the system:
velocity, Reynolds number, flow regime, friction factor, the pressure drop
split into friction / fittings / elevation, total head loss, and a velocity
warning where the design falls outside the recommended band.

Results export to CSV so the calculation is a documented record rather than a
number someone remembers.

## Method

**Pressure drop** — Darcy-Weisbach:

    ΔP = f · (L/D) · ρV²/2

**Friction factor** — Colebrook-White, solved iteratively for turbulent flow;
64/Re in the laminar regime:

    1/√f = −2·log₁₀( ε/3.7D + 2.51/(Re·√f) )

**Fitting losses** — K-factor method, coefficients per Crane TP-410:

    ΔP = ΣK · ρV²/2

**Elevation** — ΔP = ρ·g·Δh

Colebrook-White rather than a Moody chart read or a Swamee-Jain approximation:
it is the reference correlation, and iterating it costs nothing in software.
K-factors rather than equivalent-length, because equivalent-length depends on
the friction factor you are trying to find.

## Data

- **13 pipe sizes**, DN15–DN300 with NPS equivalents, Schedule 40 and 80,
  dimensions per ASME B36.10M / AS 1074
- **10 materials** with published absolute roughness — carbon steel, stainless,
  copper, PVC, PE/HDPE, galvanised, cast iron, ductile iron, concrete-lined, GRP
- **14 fluids** — water at four temperatures, air, hydraulic oil ISO 32/46/68,
  glycol mixes, diesel, mine process water, seawater
- **17 fittings** — elbows, tees, gate/globe/ball/butterfly/check valves,
  reducers, entries, exits, strainers, flexible couplings
- **12 velocity references** — water suction/discharge/gravity, fire services,
  hydraulic pressure/suction/return, compressed air, steam, mining slurry

Flow accepted in L/s, m³/h, L/min, or US GPM.

## Standards

Pipe dimensions per **ASME B36.10M** and **AS 1074**. Fitting coefficients from
**Crane TP-410**. Velocity guidance per **AS 3500** (plumbing and drainage) and
**AS 2419** (fire hydrant installations).

Applicable under **AS 4041** (Australia), **ASME B31** (US), and
**BS EN 13480** (UK/EU).

## Scope

This sizes and checks steady-state single-phase flow. It does not cover
transient analysis (water hammer, surge), two-phase or slurry rheology beyond
a velocity guideline, compressible flow at high pressure ratios, or heat
transfer along the run. Sizing output is an engineering input, not a
substitute for review by a qualified engineer against the governing code.

## License

Copyright (C) 2026 Frederik Thio.

This program is free software: you can redistribute it and/or modify it under
the terms of the **GNU Affero General Public License, version 3** or (at your
option) any later version. See [LICENSE](LICENSE) for the full text.

AGPL section 13 matters here: if you run a modified version of this program on
a server and let users interact with it over a network, you must offer those
users the source of your modified version.

Copyright is retained by the author, so the project can also be offered under
separate commercial terms in future.
