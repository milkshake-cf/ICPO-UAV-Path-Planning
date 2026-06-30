# SPSO-CPO-AGWO: UAV 3D Path Planning with Swarm Intelligence

MATLAB project extending the SPSO framework (Phung & Ha, 2021) with CPO (Abdel-Basset, 2024) and a custom AGWO variant. Research code for a first-year MS paper.

## Project

- **Stack**: MATLAB R2024b, requires Curve Fitting Toolbox
- **Root**: `d:\workspace\idea\SPSO\SPSO-main\`
- **Terrain**: `ChrismasTerrain.tif` (800×800 digital elevation model, 6 cylindrical threats)
- **Upstream**: https://github.com/duongpm/SPSO (Applied Soft Computing 2021)

## Commands

```bash
# Run a single algorithm (from project root)
/d/MATLAB/R2024b/bin/matlab -batch "SPSO_MAIN" -nodisplay -nosplash
/d/MATLAB/R2024b/bin/matlab -batch "CPO_MAIN_v2" -nodisplay -nosplash
/d/MATLAB/R2024b/bin/matlab -batch "AGWO_MAIN" -nodisplay -nosplash

# Run batch comparison (N_RUNS × 3 algorithms, ~3 min for 5 runs)
/d/MATLAB/R2024b/bin/matlab -batch "batch_compare" -nodisplay -nosplash
```

Adjust `N_RUNS` in `batch_compare.m` (line 10) and `MaxIt` per algorithm before running.

## Architecture

```
SPSO_MAIN.m          — SPSO entry: PSO with spherical vector encoding (nPop=500, MaxIt=200)
CPO_MAIN_v2.m        — Corrected CPO port (nPop=150): 4 defense strategies, clip bounds
AGWO_MAIN.m          — Adaptive GWO (nPop=150): adds pBest memory + adaptive exploration ratio
batch_compare.m      — Runs all 3 algorithms N_RUNS times, outputs stats + plots

CreateModel.m        — Loads terrain .tif, defines 6 cylindrical threats, start/end points
CreateRandomSolution.m — Random spherical vector (r, ψ, φ) for n=10 path nodes
SphericalToCart.m    — Converts spherical vectors → Cartesian (x, y, z) waypoints
MyCost.m             — Cost = b1·length + b2·threat + b3·altitude + b4·smooth (b=[5,1,10,1])
DistP2S.m            — Point-to-line-segment distance (collision check)
PlotModel.m          — 3D terrain + threat cylinder visualization
PlotSolution.m       — Plots smoothed path over terrain (uses Curve Fitting csaps)

paper_template.tex   — LaTeX paper draft with actual experiment data filled in
```

## Conventions

- **Path encoding**: spherical vectors (r, ψ, φ) — never operate directly on Cartesian coordinates during optimization
- **Cost = Inf**: collision or ground crash; initialization loops MUST retry until finding a valid solution (see `isInit` pattern in SPSO_MAIN.m)
- **Bounds**: clip to [VarMin, VarMax] for spherical components (not random reinit — that destroys search progress on this constrained problem)
- **Modifying SPSO framework files** (CreateModel, MyCost, SphericalToCart, etc.): don't — they're shared by all entry points
- **Adding a new algorithm**: copy the SPSO problem-definition block, swap the optimizer loop, keep the same CostFunction + SphericalToCart pipeline
- **Batch results** saved to `batch_comparison_results.mat`; plots to `convergence_comparison.png` and `boxplot_comparison.png`

## Notes

- CPO original MATLAB source is at `CPO_original/CPO/CPO.m` (from Abdel-Basset's MATLAB Drive). The Python-corrected `U2` fix is applied in CPO_MAIN_v2.m and AGWO_MAIN.m.
- `CPO_MAIN.m` is the user's earlier buggy attempt; use `CPO_MAIN_v2.m` instead.
- AGWO key innovations: `pBest` memory field on each agent, `explRatio = 0.7*(1-t/MaxIt)^0.5 + 0.3`, pBest-referenced position updates in all 4 defense strategies.
- Target journal: Drones (MDPI), IEEE Access, or Applied Soft Computing.
