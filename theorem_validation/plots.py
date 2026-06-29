"""Figures for the theorem-validation experiments (TV1--TV4).

Reads the CSVs written by the validate_*.py scripts and, where a phase portrait
or time series is needed, recomputes trajectories directly through ``_common``
(the same model/solver as the rest of the project).  Pure matplotlib.
"""

from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _common as C


def _read(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _f(v: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


# --- TV1: total mass attracted into the dissipativity band --------------------

def plot_tv1(scenario: str = "coexistence_P3") -> None:
    path = os.path.join(C.DATA_DIR, "tv1_total_mass_series.csv")
    if not os.path.exists(path):
        return
    rows = _read(path)
    names = [c for c in rows[0].keys() if c != "t"]
    t = [_f(r["t"]) for r in rows]

    params = C.params_from_dict(C.read_scenarios()[scenario])
    d_max = max(1.0, params.D1, params.D2, params.D3)
    d_min = min(1.0, params.D1, params.D2, params.D3)
    lower, upper = 1.0 / d_max, 1.0 / d_min

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axhspan(lower, upper, color="green", alpha=0.12,
               label=r"absorbing band $[1/D_{\max},\,1/D_{\min}]$")
    ax.axhline(lower, color="green", lw=1, ls="--")
    ax.axhline(upper, color="green", lw=1, ls="--")
    for name in names:
        ax.plot(t, [_f(r[name]) for r in rows], lw=1.2, label=f"start {name}")
    ax.set_yscale("log")
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"total mass $T=S+x+y+z$")
    ax.set_title("Experiment 1 -- every start is attracted into the bounded band "
                 r"$[1/D_{\max},\,1/D_{\min}]$")
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out = os.path.join(C.FIG_DIR, "experiment1_boundedness.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out}")


# --- TV2: eigenvalue spectra of P0..P3 for one scenario -----------------------

def plot_tv2(scenario: str = "coexistence_P3") -> None:
    params = C.params_from_dict(C.read_scenarios()[scenario])
    eqs = C.compute_all_equilibria(params)
    names = [n for n in ("P0", "P1", "P2", "P3") if eqs[n] is not None]

    fig, axes = plt.subplots(1, len(names), figsize=(3.0 * len(names), 3.2))
    if len(names) == 1:
        axes = [axes]
    for ax, name in zip(axes, names):
        eigs = C.eigenvalues(eqs[name], params)
        kind, max_re = C.classify_eigs(eigs)
        re = [e.real for e in eigs]
        im = [e.imag for e in eigs]
        stable = max_re < 0
        ax.axvline(0, color="grey", lw=1)
        ax.scatter(re, im, c=("tab:green" if stable else "tab:red"), s=40, zorder=3)
        ax.set_title(f"{name}: {kind}", fontsize=9)
        ax.set_xlabel("Re")
        ax.set_ylabel("Im")
        ax.grid(True, alpha=0.3)
    fig.suptitle(f"Experiment 2 -- local-stability spectra ({scenario}): "
                 r"left half-plane $\Rightarrow$ stable", y=1.02)
    fig.tight_layout()
    out = os.path.join(C.FIG_DIR, "experiment2_stability.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# --- Experiment 4: Hopf bifurcation (eigenvalues, Phi, amplitude) --------------

def plot_tv4_hopf(scenario: str = "oscillatory_interior_response", param: str = "a1") -> None:
    es = os.path.join(C.DATA_DIR, "tv4_hopf_eigen_sweep.csv")
    amp = os.path.join(C.DATA_DIR, "tv4_hopf_amplitude.csv")
    if not (os.path.exists(es) and os.path.exists(amp)):
        return
    erows = _read(es)
    arows = _read(amp)
    mu = [_f(r[param]) for r in erows]
    max_re = [_f(r["max_Re"]) for r in erows]
    phi = [_f(r["Phi_RH"]) for r in erows]
    amu = [_f(r[param]) for r in arows]
    amp_z = [_f(r["amp_z_late"]) for r in arows]

    # critical mu* where max Re crosses 0
    mu_star = None
    for (m0, r0), (m1, r1) in zip(zip(mu, max_re), zip(mu[1:], max_re[1:])):
        if r0 < 0 <= r1:
            mu_star = m0 + (0 - r0) * (m1 - m0) / (r1 - r0)
            break

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

    ax = axes[0]
    ax.plot(mu, max_re, "b.-", lw=1, label=r"$\max\,\mathrm{Re}\,\lambda$")
    ax.axhline(0, color="grey", lw=1)
    if mu_star is not None:
        ax.axvline(mu_star, color="red", ls="--", lw=1, label=fr"$\mu^*\approx{mu_star:.2f}$")
    ax.set_xlabel(fr"bifurcation parameter $\mu={param}$")
    ax.set_ylabel(r"$\max\,\mathrm{Re}\,\lambda$ at $P_3$")
    ax.set_title("(a) eigenvalue crossing")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(mu, phi, "m.-", lw=1)
    ax.axhline(0, color="grey", lw=1)
    if mu_star is not None:
        ax.axvline(mu_star, color="red", ls="--", lw=1)
    ax.set_xlabel(fr"$\mu={param}$")
    ax.set_ylabel(r"$\Phi=d_1d_2d_3-d_3^2-d_1^2d_4$")
    ax.set_title("(b) Routh--Hurwitz Hopf quantity")
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(amu, amp_z, "ko-", lw=1.2)
    if mu_star is not None:
        ax.axvline(mu_star, color="red", ls="--", lw=1, label=fr"$\mu^*\approx{mu_star:.2f}$")
        ax.legend(fontsize=8)
    ax.set_xlabel(fr"$\mu={param}$")
    ax.set_ylabel(r"late-window amplitude of $z$")
    ax.set_title("(c) limit-cycle bifurcation diagram")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Experiment 4 -- Hopf bifurcation of the coexistence equilibrium $P_3$", y=1.02)
    fig.tight_layout()
    out = os.path.join(C.FIG_DIR, "experiment4_hopf.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_tv4_persistence(persist: str = "coexistence_P3",
                         contrast: str = "top_predator_washout_high_D3") -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    labels = ["S", "x", "y", "z"]
    for ax, scenario, title in [(axes[0], persist, "persistent (P3 exists)"),
                                (axes[1], contrast, "non-persistent (z washes out)")]:
        params = C.params_from_dict(C.read_scenarios()[scenario])
        sim = C.integrate((1.0, 0.5, 0.3, 0.1), params, dt=0.02, t_end=400.0,
                          method="rk4", sample_every=10)
        for i, lab in enumerate(labels):
            ax.plot(sim.times, [max(s[i], 1e-12) for s in sim.trajectory], lw=1.2, label=lab)
        ax.set_yscale("log")
        ax.set_xlabel("time $t$")
        ax.set_title(f"{scenario}\n{title}", fontsize=9)
        ax.legend(fontsize=8, ncol=4)
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("component (log scale)")
    fig.suptitle("Experiment 3 -- uniform persistence vs. washout of a trophic level", y=1.02)
    fig.tight_layout()
    out = os.path.join(C.FIG_DIR, "experiment3_persistence.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    C.ensure_dirs()
    print(f"Writing report figures to {C.FIG_DIR}")
    # The report presents four experiments (boundedness, stability, persistence,
    # Hopf); only their figures are generated here.
    plot_tv1()              # Experiment 1 -- boundedness / dissipativity
    plot_tv2()              # Experiment 2 -- local stability
    plot_tv4_persistence()  # Experiment 3 -- uniform persistence
    plot_tv4_hopf()         # Experiment 4 -- Hopf bifurcation


if __name__ == "__main__":
    main()
