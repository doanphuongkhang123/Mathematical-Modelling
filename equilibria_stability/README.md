# Equilibria and Stability Threshold Analysis

This folder contains the coding and testing work for the equilibrium-analysis part of the chemostat project.

## What This Code Does

- Computes equilibria `P0`, `P1`, `P2`, and `P3`.
- Checks whether each equilibrium has biologically meaningful coordinates.
- Checks the washout stability threshold `f1(1)-D1 < 0`.
- Writes an equilibrium table to CSV.

## Model

The project uses the model that is internally consistent with the paper's Jacobian and equilibrium formulas:

```text
dS/dt = 1 - S - f1(S)x
dx/dt = x(f1(S) - D1) - f2(x)y
dy/dt = y(f2(x) - D2) - f3(y)z
dz/dt = z(f3(y) - D3)
```

The response functions are Monod functions:

```text
fi(u) = ai*u/(ki + u)
```

The displayed equation `(1.2)` in the paper appears inconsistent with its later analysis, so the implementation follows the model above.

## Input

The code reads:

```text
equilibria_stability/data/scenario_parameters.csv
```

Required columns:

```text
scenario,a1,k1,a2,k2,a3,k3,D1,D2,D3
```

Meanings:

- `a1,k1`: Monod parameters for nutrient-to-prey response `f1(S)`;
- `a2,k2`: Monod parameters for prey-to-first-predator response `f2(x)`;
- `a3,k3`: Monod parameters for first-predator-to-second-predator response `f3(y)`;
- `D1,D2,D3`: removal rates.

The parameter values are not experimental data and are not generated randomly by the script. They are fixed representative values chosen to reproduce the theoretical regimes discussed in the paper.

## Output

Run:

```bash
python3 equilibria_stability/equilibria.py
```

The script writes:

```text
equilibria_stability/outputs/equilibrium_examples.csv
```

Output columns:

```text
scenario,equilibrium,exists,S,x,y,z,P0_stable,f1(1)-D1
```

## Testing

Run the main code:

```bash
python3 equilibria_stability/equilibria.py
```

Run Python syntax checks:

```bash
python3 -m py_compile equilibria_stability/*.py
```

## Code Layout

```text
equilibria_stability/equilibria.py   command-line entry point
equilibria_stability/model.py        Monod response function and parameter objects
equilibria_stability/numerics.py     bisection root solver
equilibria_stability/analysis.py     P0-P3 formulas and P0 stability threshold
equilibria_stability/io_utils.py     CSV input/output and printed summaries
```

## Equilibrium Formulas

### P0: washout equilibrium

```text
P0 = (1, 0, 0, 0)
```

Basic stability threshold:

```text
P0 is locally stable if f1(1)-D1 < 0
```

### P1: prey-only equilibrium

```text
P1 = (S1, x1, 0, 0)
f1(S1) = D1
x1 = (1-S1)/D1
```

For Monod response:

```text
S1 = D1*k1/(a1-D1)
```

### P2: prey and first predator equilibrium

```text
P2 = (S2, x2, y2, 0)
f2(x2) = D2
S2 + f1(S2)*x2 = 1
y2 = x2*(f1(S2)-D1)/D2
```

For Monod response:

```text
x2 = D2*k2/(a2-D2)
```

### P3: coexistence equilibrium

```text
P3 = (S3, x3, y3, z3)
f3(y3) = D3
S3 + f1(S3)*x3 = 1
x3*(f1(S3)-D1) = f2(x3)*y3
z3 = y3*(f2(x3)-D2)/D3
```

For Monod response:

```text
y3 = D3*k3/(a3-D3)
```

The values `x3` and `S3` are solved numerically because the equations are coupled.
