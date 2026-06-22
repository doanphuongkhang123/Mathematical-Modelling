# Main Simulations

This package runs the five required chemostat scenarios and produces the
simulation evidence used in the report and presentation.

Run from the repository root:

```bash
python -m main_simulations.run_all
```

Default configuration:

```text
method: RK4
dt: 0.01
initial state: [0.5, 0.2, 0.1, 0.05]
P0-P3 final time: 500
oscillatory final time: 2000
CSV sampling interval: every 10 solver steps
```

Outputs are written to `results/main_scenarios/`. Use `--help` to see
configuration overrides.

Run focused tests:

```bash
python -m unittest main_simulations.test_main_simulations
```
