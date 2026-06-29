# Main Simulation Summary

## Configuration

All scenarios use the common Monod model, RK4 with `dt=0.01`, and initial
state `(S,x,y,z)=(0.5,0.2,0.1,0.05)`. The four equilibrium scenarios run to
`t=500`. The oscillatory scenario runs to `t=2000` because its damping is very
slow. No run exploded or required a negative-state clamp.

## Numerical Results

| Scenario | Predicted state | Final state `(S,x,y,z)` | Distance to target | Interpretation |
|---|---|---|---:|---|
| `washout_P0` | `P0` | `(1.000000, ~0, ~0, ~0)` | `5.55e-15` | Prey cannot offset removal; the entire chain washes out. |
| `prey_only_P1` | `P1` | `(0.166667, 1.666667, ~0, ~0)` | `2.57e-14` | Prey survives, but neither predator can invade. |
| `prey_first_predator_P2` | `P2` | `(0.628666, 0.333334, 0.409333, ~0)` | `1.29e-6` | The first predator persists; the top predator washes out. |
| `coexistence_P3` | `P3` | `(0.055067, 0.804597, 0.280534, 0.081550)` | `1.94e-11` | All four state variables approach a positive equilibrium. |
| `oscillatory_interior_response` | `P3` | `(0.377832, 0.539609, 0.008618, 0.000146)` | `2.49e-3` | The instantaneous final state lies on a slowly damped cycle around `P3`. |

The oscillatory final point should not be interpreted alone. During
`t=1500-2000`, the mean state is only `4.66e-4` from `P3`. After `t=500`, the
trajectory contains 23 detected peaks in `y` and 20 in `z`. The largest
late-window amplitude is 42.5% of the corresponding amplitude in
`t=500-1000`, confirming that the cycle is damping.

## Biological Interpretation

The first four scenarios show a trophic survival hierarchy. When prey growth
cannot exceed removal, nutrient accumulates and all organisms disappear. Once
prey can survive, higher trophic levels persist only if their food supply
supports growth faster than their own removal rates. Thus each boundary
equilibrium describes the loss of one or more upper trophic levels.

At `P3`, positive nutrient and population values coexist. This is the
simulation counterpart of the paper's persistence discussion, although one
trajectory cannot by itself prove uniform persistence for every admissible
initial condition.

The oscillatory scenario demonstrates why a sufficiently long horizon is
essential. Early in the run, `y` and especially `z` become extremely small,
which could be mistaken for permanent washout. They later recover and
continue cycling around a positive interior equilibrium. The decreasing
amplitude is consistent with a stable `P3` close to a Hopf stability boundary;
the experiment does not locate or prove a Hopf bifurcation.

## Figure Guide

- `main_scenarios_overview.png`: comparison of all five regimes.
- `*_timeseries.png`: detailed four-variable plot for each scenario.
- `oscillatory_interior_response_zoom.png`: predator oscillations after
  `t=500`, with theoretical equilibrium values shown as dashed lines.

