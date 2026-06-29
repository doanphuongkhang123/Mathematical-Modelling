# Simulator Core

Numerical ODE solvers for the chemostat food-chain model.

## Model

```
dS/dt = 1 - S - f1(S)*x
dx/dt = x*(f1(S) - D1) - f2(x)*y
dy/dt = y*(f2(x) - D2) - f3(y)*z
dz/dt = z*(f3(y) - D3)
```

## Quick Start

```bash
cd Mathematical-Modelling

# Run a scenario (RK4 by default)
python3 -m simulator.runner --scenario coexistence_P3

# Run with both RK4 and Euler
python3 -m simulator.runner --scenario prey_only_P1 --euler

# Custom time settings
python3 -m simulator.runner --scenario washout_P0 --t-end 200 --dt 0.005

# Custom initial conditions
python3 -m simulator.runner --scenario coexistence_P3 --init 0.5,0.2,0.1,0.05
```

## Output

Results go to `results/main_scenarios/`:

- `{scenario}_rk4.csv` — trajectory (`t,S,x,y,z`)
- `{scenario}_euler.csv` — (only with `--euler`)
- `{scenario}_diagnostics.json` — final state, equilibrium distances, convergence check

## CLI Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--scenario` | required | Name from `scenario_parameters.csv` |
| `--dt` | 0.01 | Time step |
| `--t-end` | 100 | Final time |
| `--euler` | off | Also run Euler solver |
| `--init` | 0.5,0.2,0.1,0.05 | Initial conditions `S0,x0,y0,z0` |
| `--input` | `../equilibria_stability/data/scenario_parameters.csv` | Parameter CSV |
| `--output-dir` | `results/main_scenarios` | Output directory |

## Available Scenarios

From `scenario_parameters.csv`:

| Scenario | Expected behavior |
|----------|-------------------|
| `washout_P0` | S→1, all species wash out |
| `prey_only_P1` | S,x positive; y,z→0 |
| `prey_first_predator_P2` | S,x,y positive; z→0 |
| `coexistence_P3` | All four positive |
| `oscillatory_interior_response` | Oscillatory approach to P3 |
| `washout_high_prey_removal` | S→1, all wash out (high D1) |
| `washout_weak_prey_growth` | S→1, all wash out (weak growth) |
| `prey_only_predator_limited` | Prey only (weak predator) |
| `prey_first_predator_strong_y` | Prey + first predator |
| `coexistence_low_D3` | Coexistence (low D3) |
| `coexistence_high_D3` | Coexistence (high D3) |
| `top_predator_washout_high_D3` | Top predator washes out |
| `coexistence_alt_parameters` | Coexistence (alt params) |

## Running Tests

```bash
cd Mathematical-Modelling

# Run the complete simulator test suite
python3 -m unittest discover -s simulator/tests -t .

# Or run each module separately
python3 -m simulator.tests.test_rhs
python3 -m simulator.tests.test_solvers
python3 -m simulator.tests.test_diagnostics
```

## Using as a Library

```python
import sys
sys.path.insert(0, "Mathematical-Modelling")

from equilibria_stability.model import params_from_dict
from simulator.rhs import rhs
from simulator.solvers import rk4
from simulator.diagnostics import diagnose

params = params_from_dict({
    "a1": 4.8, "k1": 0.17, "a2": 4.3, "k2": 1.85,
    "a3": 4.6, "k3": 1.66, "D1": 0.72, "D2": 1.11, "D3": 0.665,
})

rhs_fn = lambda s: rhs(s, params)
times, states, info = rk4(rhs_fn, [0.5, 0.2, 0.1, 0.05], t_end=100, dt=0.01)
diag = diagnose(times, states, info, params)

print(diag["final_state"])
print(diag["equilibrium_proximity"])
```
