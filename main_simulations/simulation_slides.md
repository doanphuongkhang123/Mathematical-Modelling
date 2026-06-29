# Main Simulation Slides

## Slide 1: Simulation Setup and Required Scenarios

**On slide**

- Monod chemostat model with state `(S,x,y,z)`
- RK4, `dt=0.01`, initial state `(0.5,0.2,0.1,0.05)`
- `t=500` for `P0-P3`; `t=2000` for slow oscillations
- Five regimes: washout, prey-only, two-level survival, coexistence,
  oscillatory coexistence
- Figure: `main_scenarios_overview.png`

**Speaker notes**

We used the same positive initial state for every scenario, so the differences
come from the parameters rather than the starting populations. The numerical
states were compared directly with equilibria computed by Person 2. No run
exploded or required negative-value correction.

## Slide 2: Washout Hierarchy from P0 to P2

**On slide**

| Regime | Long-run survivors | Distance to equilibrium |
|---|---|---:|
| `P0` | none | `5.55e-15` |
| `P1` | prey `x` | `2.57e-14` |
| `P2` | prey `x`, predator `y` | `1.29e-6` |

- Each upper trophic level survives only if growth exceeds removal.
- Figure: use the `P0-P2` panels from `main_scenarios_overview.png`.

**Speaker notes**

At `P0`, nutrient returns to one because no organism consumes it. At `P1`,
prey survives but cannot support the first predator. At `P2`, the first
predator survives, but the top predator washes out. This is the biological
meaning of the boundary equilibria and their invasion conditions.

## Slide 3: Full Coexistence at P3

**On slide**

```text
Computed: (0.055067, 0.804597, 0.280534, 0.081550)
Theory:   (0.055067, 0.804597, 0.280534, 0.081550)
Distance: 1.94e-11
```

- All trophic levels remain positive.
- Numerical evidence agrees with the interior equilibrium.
- Figure: `coexistence_P3_timeseries.png`

**Speaker notes**

This scenario illustrates coexistence and supports the persistence
interpretation. We should distinguish simulation evidence from proof:
uniform persistence is an analytical statement about all relevant interior
trajectories, not just this one run.

## Slide 4: Slowly Damped Oscillation Around P3

**On slide**

- 23 peaks in `y`, 20 peaks in `z` after `t=500`
- Late-window mean distance to `P3`: `4.66e-4`
- Late/middle amplitude ratio: at most `0.425`
- Conclusion: slowly damped oscillation, not predator washout
- Figure: `oscillatory_interior_response_zoom.png`

**Speaker notes**

Stopping at `t=100` would be misleading because the top predator becomes
extremely small before recovering. The long run reveals repeated cycles whose
amplitude decreases. This behavior is consistent with a stable `P3` near a
Hopf boundary, but we did not calculate the bifurcation point itself.
