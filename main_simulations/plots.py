"""Plotting helpers for the main chemostat simulations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt


COMPONENTS = ("S", "x", "y", "z")
LABELS = ("S (nutrient)", "x (prey)", "y (predator 1)", "z (predator 2)")
COLORS = ("#1f77b4", "#2ca02c", "#f28e2b", "#d62728")


def _plot_panel(ax, times, states, component, target, title=None):
    ax.plot(
        times,
        [state[component] for state in states],
        color=COLORS[component],
        linewidth=1.1,
    )
    ax.axhline(target[component], color="#666666", linestyle="--", linewidth=1.0)
    ax.set_ylabel(LABELS[component])
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.25)


def create_scenario_plot(name: str, data: dict, output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), sharex=True)
    for component, ax in enumerate(axes.flat):
        _plot_panel(
            ax,
            data["times"],
            data["states"],
            component,
            data["target"],
        )
    axes[1, 0].set_xlabel("time t")
    axes[1, 1].set_xlabel("time t")
    fig.suptitle(f"{name}: {data['spec'].description}\n(dashed line: theoretical equilibrium)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(output_dir / f"{name}_timeseries.png", dpi=160)
    plt.close(fig)


def create_oscillatory_zoom(data: dict, output_dir: Path) -> None:
    selected = [
        (time, state)
        for time, state in zip(data["times"], data["states"])
        if time >= 500.0
    ]
    times = [item[0] for item in selected]
    states = [item[1] for item in selected]
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 6.2), sharex=True)
    for component in (2, 3):
        axes[component - 2].plot(
            times,
            [state[component] for state in states],
            color=COLORS[component],
            linewidth=1.0,
        )
        axes[component - 2].axhline(
            data["target"][component],
            color="#666666",
            linestyle="--",
            linewidth=1.0,
        )
        axes[component - 2].set_ylabel(LABELS[component])
        axes[component - 2].grid(True, alpha=0.25)
    axes[1].set_xlabel("time t")
    fig.suptitle("Oscillatory scenario: slowly decreasing predator-cycle amplitude")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_dir / "oscillatory_interior_response_zoom.png", dpi=160)
    plt.close(fig)


def create_overview(plot_data: dict, output_dir: Path) -> None:
    fig, axes = plt.subplots(5, 1, figsize=(10.0, 13.0))
    for ax, (name, data) in zip(axes, plot_data.items()):
        every = max(1, len(data["times"]) // 5000)
        times = data["times"][::every]
        states = data["states"][::every]
        for component in range(4):
            ax.plot(
                times,
                [state[component] for state in states],
                color=COLORS[component],
                linewidth=0.9,
                label=COMPONENTS[component],
            )
        ax.set_title(f"{name}: expected {data['spec'].expected_equilibrium}", fontsize=10)
        ax.set_ylabel("state")
        ax.grid(True, alpha=0.2)
    axes[-1].set_xlabel("time t")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle("Main chemostat simulation scenarios", fontsize=14)
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    fig.savefig(output_dir / "main_scenarios_overview.png", dpi=160)
    plt.close(fig)


def create_all_plots(plot_data: dict, output_dir: Path) -> None:
    for name, data in plot_data.items():
        create_scenario_plot(name, data, output_dir)
    create_oscillatory_zoom(plot_data["oscillatory_interior_response"], output_dir)
    create_overview(plot_data, output_dir)
