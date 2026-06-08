# Robustness and Sensitivity Analysis



## What it computes

| Experiment | Question | Main output |
|---|---|---|
| **1. Parameter sweeps** | When does each population survive vs wash out as `D1`, `D2`, `D3` vary? | `outputs/sensitivity/sweep_*.csv`, `figures/sweep_*.png` |
| **2. Elasticity coefficients** | How sensitive is the coexistence equilibrium to each parameter (1 % perturbation)? | `outputs/sensitivity/elasticity_coefficients.csv`, `figures/elasticity_heatmap.png` |
| **3. Step‑size robustness** | Does the answer depend on `dt` or on Euler vs RK4? | `outputs/robustness/step_size.csv`, `figures/step_size_convergence.png` |
| **4. Initial‑condition robustness** | Does the long‑run outcome depend on the starting state? | `outputs/robustness/initial_conditions.csv`, `figures/initial_conditions_timeseries.png` |


## Quick start — run everything

Run everything from main root:

```bash
python3 robustness_sensitivity/run_all.py
```

This runs Experiments 1–4 and writes all CSVs and figures under
`robustness_sensitivity/outputs/`. To skip the figures:

```bash
python3 robustness_sensitivity/run_all.py --no-plots
```

To analyse a different baseline scenario from `data/baseline_scenarios.csv`:

```bash
python3 robustness_sensitivity/run_all.py --baseline prey_first_predator_P2
```

## Running each experiment on its own

```bash
# Experiments 1 + 2 (sensitivity)
python3 robustness_sensitivity/sensitivity.py

# Experiments 3 + 4 (robustness)
python3 robustness_sensitivity/robustness.py

# Figures (reads the CSVs produced above)
python3 robustness_sensitivity/plots.py
```

### Useful options

`sensitivity.py`
```
--baseline NAME     baseline scenario (default: coexistence_P3)
--params D1,D2,D3   which removal rates to sweep (Experiment A)
--rel-delta 0.01    relative perturbation for elasticity (default 1%)
--forward           use forward differences (course style) instead of central
--method rk4|euler  integrator used for the simulation cross-check
--dt 0.05 --t-end 400   integration step / horizon for the sweep cross-check
```

`robustness.py`
```
--baseline NAME     baseline scenario (default: coexistence_P3)
--dt0 0.4 --levels 6   coarsest step and number of halvings for Experiment C
--t-acc 2.0         short horizon for the convergence-order estimate
--t-end 400         long horizon for the regime / initial-condition study
--method rk4|euler  integrator for the initial-condition study
```

## Input

`data/baseline_scenarios.csv` — fixed, non‑random parameter sets:

```
scenario,a1,k1,a2,k2,a3,k3,D1,D2,D3,S0,x0,y0,z0
```

`a1,k1` are the Monod parameters of `f1(S)` (nutrient→prey), `a2,k2` of `f2(x)`
(prey→predator 1), `a3,k3` of `f3(y)` (predator 1→predator 2); `D1,D2,D3` are the
removal rates; `S0,x0,y0,z0` is the default initial state.

## Output files

```
outputs/
├── sensitivity/
│   ├── sweep_D1.csv  sweep_D2.csv  sweep_D3.csv   # regime + coords vs each rate
│   ├── sweep_all.csv                              # all sweeps stacked
│   └── elasticity_coefficients.csv                # 9×4 sensitivity matrix
├── robustness/
│   ├── step_size.csv                              # self-convergence + stability
│   └── initial_conditions.csv                     # outcome from each start
└── figures/
    ├── sweep_D1.png  sweep_D2.png  sweep_D3.png
    ├── elasticity_heatmap.png
    ├── step_size_convergence.png
    └── initial_conditions_timeseries.png
```



