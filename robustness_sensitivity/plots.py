from __future__ import annotations
import argparse
import csv
import math
import matplotlib.pyplot as plt
import os
from typing import Optional
import chemostat as cx
from sensitivity import DEFAULT_DATA, DEFAULT_STATE0, OUTPUT_LABELS, read_baselines
from robustness import DEFAULT_INITIAL_CONDITIONS

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SENS_DIR = os.path.join(THIS_DIR, "outputs", "sensitivity")
ROB_DIR = os.path.join(THIS_DIR, "outputs", "robustness")
FIG_DIR = os.path.join(THIS_DIR, "outputs", "figures")


def _read_csv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _f(value: str) -> float:
    if value is None or value == "" or value.lower() == "nan":
        return float("nan")
    return float(value)



def plot_sweeps(plt):
    for param in ("D1", "D2", "D3"):
        path = os.path.join(SENS_DIR, f"sweep_{param}.csv")
        if not os.path.exists(path):
            continue
        rows = _read_csv(path)
        values = [_f(r["value"]) for r in rows]
        theory = {c: [_f(r[f"theory_{c}"]) for r in rows] for c in ("S", "x", "y", "z")}
        sim_z = [_f(r["sim_z"]) for r in rows]
        p3 = [r["P3_exists"] == "True" for r in rows]

        fig, ax = plt.subplots(figsize=(7.0, 4.3))
        ax.plot(values, theory["S"], label="S (nutrient)", color="blue")
        ax.plot(values, theory["x"], label="x (prey)", color="green")
        ax.plot(values, theory["y"], label="y (predator 1)", color="orange")
        ax.plot(values, theory["z"], label="z (predator 2)", color="red", linewidth=2)
        ax.plot(values, sim_z, "k.", markersize=4, label="z (simulated)")


        threshold = next((values[i] for i in range(len(p3))
                          if p3[i] and i + 1 < len(p3) and not p3[i + 1]), None)
        if threshold is not None:
            ax.axvline(threshold, color="grey", linestyle="--", linewidth=1)
            ax.text(threshold, ax.get_ylim()[1] * 0.92, f"  z washes out\n  {param}≈{threshold:.3g}",
                    fontsize=8, color="grey", va="top")

        ax.set_xlabel(param)
        ax.set_ylabel("equilibrium coordinate")
        ax.set_title(f"Sensitivity sweep: equilibrium vs {param}")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        out = os.path.join(FIG_DIR, f"sweep_{param}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"  wrote {out}")


def plot_elasticity(plt):
    """Heatmap of the elasticity matrix S(output, param)."""
    path = os.path.join(SENS_DIR, "elasticity_coefficients.csv")
    if not os.path.exists(path):
        return
    rows = _read_csv(path)
    params = [r["param"] for r in rows]
    matrix = [[_f(r[label]) for label in OUTPUT_LABELS] for r in rows]

    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    vmax = max((abs(v) for row in matrix for v in row if v == v), default=1.0)
    im = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(OUTPUT_LABELS)))
    ax.set_xticklabels(OUTPUT_LABELS)
    ax.set_yticks(range(len(params)))
    ax.set_yticklabels(params)
    ax.set_xlabel("equilibrium output")
    ax.set_ylabel("parameter")
    ax.set_title("Elasticity  S(output, param)\nat the coexistence equilibrium P3")
    for i in range(len(params)):
        for j in range(len(OUTPUT_LABELS)):
            v = matrix[i][j]
            if v == v:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                        color="black" if abs(v) < 0.6 * vmax else "white")
    fig.colorbar(im, ax=ax, label="relative sensitivity")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "elasticity_heatmap.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out}")


def plot_step_size(plt):
    """Log-log self-convergence plot with reference slopes for Euler/RK4."""
    path = os.path.join(ROB_DIR, "step_size.csv")
    if not os.path.exists(path):
        return
    rows = _read_csv(path)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    colors = {"rk4": "blue", "euler": "red"}
    for method in ("rk4", "euler"):
        pts = [(_f(r["dt"]), _f(r["diff_to_half_step"])) for r in rows
               if r["method"] == method and _f(r["diff_to_half_step"]) == _f(r["diff_to_half_step"])
               and _f(r["diff_to_half_step"]) > 0]
        if not pts:
            continue
        dts, diffs = zip(*pts)
        ax.loglog(dts, diffs, "o-", color=colors[method], label=f"{method.upper()} self-convergence")

    # Reference slopes (order 1 and order 4) anchored to the plot.
    dt_ref = [r for r in rows if _f(r["diff_to_half_step"]) == _f(r["diff_to_half_step"])]
    if dt_ref:
        anchors = sorted({_f(r["dt"]) for r in dt_ref})
        d0, d1 = anchors[0], anchors[-1]
        base = max(_f(r["diff_to_half_step"]) for r in dt_ref if _f(r["dt"]) == d0)
        ax.loglog([d0, d1], [base, base * (d1 / d0) ** 1], "k--", linewidth=1, label="slope 1 (Euler)")
        ax.loglog([d0, d1], [base, base * (d1 / d0) ** 4], "k:", linewidth=1, label="slope 4 (RK4)")

    ax.set_xlabel("step size dt")
    ax.set_ylabel("|| u(dt) − u(dt/2) ||  at t_acc")
    ax.set_title("Step-size robustness / order of convergence")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "step_size_convergence.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out}")


def plot_initial_conditions(plt, baseline: str):
    """Time series from several initial conditions, showing convergence."""
    baselines = read_baselines(DEFAULT_DATA)
    if baseline not in baselines:
        return
    params = cx.params_from_dict(baselines[baseline])
    equilibrium = cx.compute_all_equilibria(params)
    target = next((equilibrium[k] for k in ("P3", "P2", "P1", "P0")
                   if equilibrium[k] is not None), None)

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0), sharex=True)
    names = ["S (nutrient)", "x (prey)", "y (predator 1)", "z (predator 2)"]
    for state0 in DEFAULT_INITIAL_CONDITIONS:
        sim = cx.integrate(state0, params, dt=0.02, t_end=80.0, method="rk4", sample_every=10)
        for k in range(4):
            ax = axes[k // 2][k % 2]
            ax.plot(sim.times, [s[k] for s in sim.trajectory], linewidth=1,
                    label=f"start {tuple(round(v, 2) for v in state0)}")
    for k in range(4):
        ax = axes[k // 2][k % 2]
        if target is not None:
            ax.axhline(target[k], color="grey", linestyle="--", linewidth=1)
        ax.set_ylabel(names[k])
        ax.grid(True, alpha=0.3)
    axes[1][0].set_xlabel("time t")
    axes[1][1].set_xlabel("time t")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.suptitle(f"Initial-condition robustness ({baseline}):\nall starts reach the same equilibrium (dashed line)",
                 y=0.99)
    fig.legend(handles, labels, fontsize=7, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    out = os.path.join(FIG_DIR, "initial_conditions_timeseries.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    parser = argparse.ArgumentParser(description="Plot the robustness/sensitivity results.")
    parser.add_argument("--baseline", default="coexistence_P3")
    args = parser.parse_args()

    

    os.makedirs(FIG_DIR, exist_ok=True)
    print(f"Writing figures to {FIG_DIR}")
    plot_sweeps(plt)
    plot_elasticity(plt)
    plot_step_size(plt)
    plot_initial_conditions(plt, args.baseline)


if __name__ == "__main__":
    main()
