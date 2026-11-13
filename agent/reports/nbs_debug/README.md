# NBS debugging and check scripts (2026-11)

Scripts used while fixing and restructuring the NBS module. See `agent/NBS_INVARIANTS.md` ("Open follow-ups").

| Script | What it does | Run with |
|---|---|---|
| `run_tests.jl` | Runs every NBS/dynamics test file plus the `Upstream area`, `Sequencing` and `Sequence runoff matches watercourses` testsets, each isolated, so one failure cannot hide the others. (`Pkg.test()` also runs the whole suite now: the known-broken checks are wrapped or `@test_broken`.) | `julia --project=. agent/reports/nbs_debug/run_tests.jl` (from the repo root) |
| `setup.jl` | Helpers to rebuild the SUrbArea webapp's synthetic NBS case outside the webapp (DTM, green-roof placements, gauges, mm/h → m/min rain). Included by the two scripts below. | — |
| `run1.jl` | Runs both synthetic layouts (one roof, three roofs) through `fill_sequence`, and saves `run1.jld2`. | `julia --project=<SUrbArea>/WebInterface run1.jl` (from this folder) |
| `run2.jl` | Samples both runs every 10 min, prints the global water balance (rain = storage change + outflow), and compares gauge flows with and without NBS. Needs `run1.jld2`. | as above |
| `zero_vol.jl` | Counts zero-volume traps on the test terrains and the synthetic DTM (the basis for invariant S6). | as above |

**Inputs.** These come from the SUrbArea checkout, which is not in git: `<SUrbArea>/debug/synthetic_dtm_flattened.jld2` and `<SUrbArea>/debug/synthetic_userdefined*.geojson`. The SUrbArea checkout is assumed to sit next to this repository; otherwise set `SURBAREA_DIR`.

**Rain.** The runs use 120 mm/h for 1000 min, which matches the levels in the July gauge CSVs in `<SUrbArea>/debug/`.
