# Mathematical-Modelling

## Real-time demo: Interactive Food Chain

An interactive, browser-based real-time demo lives in [`demo/`](demo/). 

Run:

```bash
python3 -m demo.app
# then open http://127.0.0.1:8000
```



## Equilibria and Stability Threshold Analysis

Equilibrium-analysis code is in [`equilibria_stability/`](equilibria_stability/).

Run:

```bash
python3 equilibria_stability/equilibria.py
```

This computes the equilibria `P0`, `P1`, `P2`, and `P3` for the project scenarios and writes:

```text
equilibria_stability/outputs/equilibrium_examples.csv
```


## Run experiment for robustness and sensitivity

Equilibrium-analysis code is in [`robustness_sensitivity/`](robustness_sensitivity/).

Run everything from main root:

```bash
python3 robustness_sensitivity/run_all.py
```

Run separate Experiment:

```bash
# Experiments 1 + 2 (sensitivity)
python3 robustness_sensitivity/sensitivity.py

# Experiments 3 + 4 (robustness)
python3 robustness_sensitivity/robustness.py

# Figures (reads the CSVs produced above)
python3 robustness_sensitivity/plots.py
```



