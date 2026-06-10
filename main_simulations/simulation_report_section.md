# Simulation Setup, Results, and Biological Interpretation

## Simulation Setup

The numerical experiments were designed to reproduce the principal long-run
regimes implied by the equilibrium and stability analysis. We used the
internally consistent four-dimensional chemostat system

```text
dS/dt = 1 - S - f1(S)x,
dx/dt = x(f1(S)-D1) - f2(x)y,
dy/dt = y(f2(x)-D2) - f3(y)z,
dz/dt = z(f3(y)-D3),
```

with Monod response functions `fi(u)=ai*u/(ki+u)`. This version agrees with
the Jacobian and equilibrium formulas in the reference. The displayed
equation (1.2) in the paper contains inconsistent extra population factors,
so it was not used.

Five parameter sets were selected from
`equilibria_stability/data/scenario_parameters.csv`: washout `P0`, prey-only
`P1`, prey plus first predator `P2`, full coexistence `P3`, and an interior
response with slowly damped oscillations. Every experiment started from the
same positive state,

```text
(S0,x0,y0,z0) = (0.5,0.2,0.1,0.05).
```

The system was integrated with the classical fourth-order Runge-Kutta method
using time step `dt=0.01`. The `P0-P3` experiments were run to `t=500`. The
oscillatory experiment was extended to `t=2000` because a shorter horizon can
give a misleading impression of predator washout. Full-resolution solutions
were used for diagnostics, while every tenth point was written to CSV. The
implementation recorded minimum and maximum values, negative-state clamps,
explosion status, final states, and Euclidean distances from all existing
analytical equilibria.

The simulations produced no numerical explosion and required no negative
clamping. This agrees with the paper's theoretical statements that the
non-negative state space is invariant and that solutions are eventually
bounded. It also indicates that the selected time step was appropriate for
these experiments. Step-size robustness is addressed separately in the
group's robustness analysis.

Each equilibrium scenario was classified by computing the Euclidean distance
between the numerical final state and every analytical equilibrium that
exists for the same parameter set. The equilibrium with the smallest distance
was recorded as the nearest equilibrium and compared with the regime expected
from the parameter design. A tolerance of `10^-4` was used as the acceptance
criterion for the four equilibrium cases. This procedure is more reliable
than deciding from a graph alone, especially when one population is small or
when oscillatory transients remain visible.

The oscillatory case required a different criterion because its final state
depends on the phase at which integration stops. We compared the mean state
with `P3`, counted local maxima after `t=500`, and measured peak-to-trough
amplitudes in the windows `t=500-1000` and `t=1500-2000`. The trajectory was
classified as a damped oscillation only when the late mean was close to `P3`,
multiple cycles were present, and the late amplitude was below 75% of the
middle-window amplitude. This distinction prevents a temporary low predator
concentration from being incorrectly labeled as washout.

## Equilibrium Scenarios

In the washout scenario, the final state was approximately
`(1,0,0,0)`, with Euclidean distance `5.55e-15` from `P0`. The nutrient
returned to its normalized input concentration while every organism
disappeared. This occurs when prey growth at the available nutrient level
cannot compensate for prey removal. Since the predators depend indirectly or
directly on the prey, prey washout propagates through the whole chain. The
correct local stability condition is `f1(1)-D1<0`, as follows from Theorem 2.2
and the eigenvalues at `P0`. The opposite inequality printed in Section 5 of
the paper is a typographical error.

For the prey-only scenario, the computed state was
`(0.166667,1.666667,0,0)`. Its distance from `P1` was `2.57e-14`. Here the prey
can exploit the incoming nutrient, but the prey concentration is insufficient
for the first predator to overcome its removal rate. Consequently, both
predators wash out. The reduction in nutrient from one to approximately
`0.1667` also shows that prey persistence is associated with continued
nutrient consumption.

The `P2` scenario converged to
`(0.628666,0.333334,0.409333,0)`, only `1.29e-6` from the analytical `P2`
equilibrium. The first predator survives by consuming the prey, whereas the
second predator cannot invade and its concentration approaches zero. The
trajectory approaches `P2` through visible damped oscillations, consistent
with the complex stable eigenvalue pair associated with this parameter set.
This case demonstrates that transient oscillation does not necessarily imply
a sustained periodic orbit.

For the coexistence scenario, the solution converged rapidly to
`(0.055067,0.804597,0.280534,0.081550)`. The final distance from `P3` was
`1.94e-11`. All trophic levels remain positive: nutrient enters continuously,
the prey consumes nutrient, the first predator consumes prey, and the second
predator consumes the first predator. The numerical behavior supports the
paper's persistence interpretation for this selected parameter set and
initial state. However, a single trajectory is not a proof of uniform
persistence, which is a statement about all relevant interior solutions and
requires the analytical conditions developed in Section 3.

## Slowly Damped Oscillatory Scenario

The fifth scenario required special treatment. At intermediate times, the
upper predators can become extremely small, particularly `z`. If the
simulation were stopped near `t=100`, the trajectory could be misclassified
as prey-only. A longer integration shows repeated recovery and decline of the
predator populations.

The theoretical interior equilibrium is
`P3=(0.377311,0.540994,0.006908,0.001193)`. The instantaneous state at
`t=2000` was `(0.377832,0.539609,0.008618,0.000146)`, with distance
`2.49e-3` from `P3`. Since an oscillating trajectory should not be assessed
from one phase point, we also examined window averages and amplitudes. The
mean over `t=1500-2000` was only `4.66e-4` from `P3`. After `t=500`, 23 local
maxima were detected in `y` and 20 in `z`. For every state variable, the
amplitude in `t=1500-2000` was less than half the corresponding amplitude in
`t=500-1000`; the maximum ratio was `0.425`.

These observations identify a slowly damped oscillatory approach to `P3`.
They are qualitatively consistent with an interior equilibrium located close
to a Hopf stability boundary. They do not establish a Hopf bifurcation,
because locating such a bifurcation would require varying a designated
parameter, tracking eigenvalues, and verifying the transversality condition.

## Overall Interpretation

Together, the simulations provide a numerical picture of the equilibrium
hierarchy described by the analysis. Survival progresses from no organisms
at `P0`, through one and two surviving biological populations at `P1` and
`P2`, to complete coexistence at `P3`. The controlling biological mechanism
is whether each trophic level can grow rapidly enough on its food source to
offset its distinct removal rate.

The experiments also show why analytical equilibria and sufficiently long
simulations must be used together. Final-state comparison accurately
classifies stable equilibrium regimes, while time-window statistics are
needed for weakly damped oscillations. The agreement between the computed
states and analytical equilibria supports the correctness of the simulator
and gives concrete biological meaning to the local stability, washout,
persistence, and bifurcation ideas developed in the reference.
