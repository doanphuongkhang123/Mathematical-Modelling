#!/usr/bin/env python3
"""Run and summarize the five main chemostat scenarios."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EQUILIBRIA_DIR = PROJECT_ROOT / "equilibria_stability"
for import_path in (PROJECT_ROOT, EQUILIBRIA_DIR):
    if str(import_path) not in sys.path:
        sys.path.insert(0, str(import_path))

from analysis import compute_all_equilibria  # noqa: E402
from model import PARAMETER_NAMES, Params, params_from_dict  # noqa: E402
from simulator.rhs import rhs  # noqa: E402
from simulator.solvers import rk4  # noqa: E402


DEFAULT_INPUT = PROJECT_ROOT / "equilibria_stability" / "data" / "scenario_parameters.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "main_scenarios"
DEFAULT_STATE = (0.5, 0.2, 0.1, 0.05)
COMPONENTS = ("S", "x", "y", "z")


@dataclass(frozen=True)
class ScenarioSpec:
    name: str
    expected_equilibrium: str
    description: str
    oscillatory: bool = False


SCENARIOS = (
    ScenarioSpec("washout_P0", "P0", "Complete organism washout"),
    ScenarioSpec("prey_only_P1", "P1", "Prey survives; both predators wash out"),
    ScenarioSpec("prey_first_predator_P2", "P2", "Prey and first predator survive"),
    ScenarioSpec("coexistence_P3", "P3", "All trophic levels coexist"),
    ScenarioSpec(
        "oscillatory_interior_response",
        "P3",
        "Slowly damped oscillation around coexistence",
        oscillatory=True,
    ),
)


def parse_state(raw: str) -> tuple[float, float, float, float]:
    values = tuple(float(value.strip()) for value in raw.split(","))
    if len(values) != 4:
        raise argparse.ArgumentTypeError("--init requires exactly four values: S0,x0,y0,z0")
    if any(not math.isfinite(value) or value < 0.0 for value in values):
        raise argparse.ArgumentTypeError("initial values must be finite and non-negative")
    return values  # type: ignore[return-value]


def validate_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive, got {value}")


def read_scenarios(path: Path) -> dict[str, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {}
        for row in csv.DictReader(handle):
            rows[row["scenario"]] = {key: float(row[key]) for key in PARAMETER_NAMES}
    missing = [spec.name for spec in SCENARIOS if spec.name not in rows]
    if missing:
        raise ValueError(f"Missing required scenarios in {path}: {', '.join(missing)}")
    return rows


def sample_indices(length: int, every: int) -> list[int]:
    if length <= 0:
        return []
    if every <= 0:
        raise ValueError("sample_every must be positive")
    indices = list(range(0, length, every))
    if indices[-1] != length - 1:
        indices.append(length - 1)
    return indices


def euclidean(left: Iterable[float], right: Iterable[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def window_stats(
    times: list[float],
    states: list[list[float]],
    start: float,
    end: float,
) -> tuple[list[float], list[float]]:
    window = [state for time, state in zip(times, states) if start <= time <= end]
    if not window:
        raise ValueError(f"No trajectory samples in window [{start}, {end}]")
    means = [sum(state[i] for state in window) / len(window) for i in range(4)]
    amplitudes = [
        max(state[i] for state in window) - min(state[i] for state in window)
        for i in range(4)
    ]
    return means, amplitudes


def count_peaks(values: list[float], prominence: float = 1e-8) -> int:
    peaks = 0
    for index in range(1, len(values) - 1):
        rise = values[index] - values[index - 1]
        fall = values[index] - values[index + 1]
        if rise > prominence and fall > prominence:
            peaks += 1
    return peaks


def oscillation_criteria(
    late_mean_distance: float,
    peak_count_y: int,
    peak_count_z: int,
    amplitude_ratio: float,
) -> dict[str, bool]:
    return {
        "late_mean_near_P3": late_mean_distance < 1e-3,
        "enough_peaks": max(peak_count_y, peak_count_z) >= 5,
        "amplitude_decreased": amplitude_ratio < 0.75,
    }


def write_trajectory(
    path: Path,
    times: list[float],
    states: list[list[float]],
    every: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    indices = sample_indices(len(times), every)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("t", "S", "x", "y", "z"))
        for index in indices:
            writer.writerow(
                (format(times[index], ".12g"),)
                + tuple(format(value, ".12g") for value in states[index])
            )


def analyze_scenario(
    spec: ScenarioSpec,
    params: Params,
    times: list[float],
    states: list[list[float]],
    solver_info: dict,
) -> tuple[dict, dict]:
    equilibria = compute_all_equilibria(params)
    final = states[-1]
    distances = {
        name: (euclidean(final, equilibrium) if equilibrium is not None else None)
        for name, equilibrium in equilibria.items()
    }
    existing_distances = {name: value for name, value in distances.items() if value is not None}
    nearest = min(existing_distances, key=existing_distances.get)  # type: ignore[arg-type]
    target = equilibria[spec.expected_equilibrium]
    if target is None:
        raise ValueError(f"{spec.name}: expected equilibrium {spec.expected_equilibrium} does not exist")

    tail_start = 0.9 * times[-1]
    tail_mean, tail_amplitude = window_stats(times, states, tail_start, times[-1])
    peak_y = 0
    peak_z = 0
    damping_ratio = None
    oscillation_checks = None

    if spec.oscillatory:
        middle_mean, middle_amplitude = window_stats(times, states, 500.0, 1000.0)
        late_mean, late_amplitude = window_stats(times, states, 1500.0, 2000.0)
        after_500 = [state for time, state in zip(times, states) if time >= 500.0]
        peak_y = count_peaks([state[2] for state in after_500])
        peak_z = count_peaks([state[3] for state in after_500])
        active = [i for i in range(4) if middle_amplitude[i] > 1e-12]
        ratios = [late_amplitude[i] / middle_amplitude[i] for i in active]
        damping_ratio = max(ratios)
        late_mean_distance = euclidean(late_mean, target)
        criteria = oscillation_criteria(
            late_mean_distance,
            peak_y,
            peak_z,
            damping_ratio,
        )
        oscillation_checks = {
            "late_mean_distance_to_P3": late_mean_distance,
            "peak_count_y_after_t500": peak_y,
            "peak_count_z_after_t500": peak_z,
            "maximum_late_to_middle_amplitude_ratio": damping_ratio,
            **criteria,
            "middle_mean": dict(zip(COMPONENTS, middle_mean)),
            "middle_amplitude": dict(zip(COMPONENTS, middle_amplitude)),
            "late_mean": dict(zip(COMPONENTS, late_mean)),
            "late_amplitude": dict(zip(COMPONENTS, late_amplitude)),
        }
        outcome = (
            "damped oscillation around P3"
            if all(criteria.values())
            else "oscillation criteria not met"
        )
    else:
        outcome = spec.description

    diagnostics = {
        "scenario": spec.name,
        "method": "rk4",
        "initial_state": dict(zip(COMPONENTS, states[0])),
        "final_state": dict(zip(COMPONENTS, final)),
        "min_state": {
            component: min(state[i] for state in states)
            for i, component in enumerate(COMPONENTS)
        },
        "max_state": {
            component: max(state[i] for state in states)
            for i, component in enumerate(COMPONENTS)
        },
        "tail_mean": dict(zip(COMPONENTS, tail_mean)),
        "tail_amplitude": dict(zip(COMPONENTS, tail_amplitude)),
        "equilibrium_distances": distances,
        "expected_equilibrium": spec.expected_equilibrium,
        "nearest_equilibrium": nearest,
        "target_equilibrium": dict(zip(COMPONENTS, target)),
        "target_distance": euclidean(final, target),
        "negative_clamps": solver_info["total_clamps"],
        "exploded": solver_info["exploded"],
        "steps": solver_info["steps"],
        "outcome": outcome,
        "oscillation_checks": oscillation_checks,
    }
    summary = {
        "scenario": spec.name,
        "expected_equilibrium": spec.expected_equilibrium,
        "nearest_equilibrium": nearest,
        "outcome": outcome,
        **{f"final_{name}": final[i] for i, name in enumerate(COMPONENTS)},
        **{f"target_{name}": target[i] for i, name in enumerate(COMPONENTS)},
        "target_distance": euclidean(final, target),
        **{f"tail_mean_{name}": tail_mean[i] for i, name in enumerate(COMPONENTS)},
        **{f"tail_amplitude_{name}": tail_amplitude[i] for i, name in enumerate(COMPONENTS)},
        "peak_count_y": peak_y,
        "peak_count_z": peak_z,
        "damping_ratio": damping_ratio,
        "negative_clamps": solver_info["total_clamps"],
        "exploded": solver_info["exploded"],
    }
    return diagnostics, summary


def write_summary(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: format(value, ".12g") if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--init", type=parse_state, default=DEFAULT_STATE)
    parser.add_argument("--equilibrium-t-end", type=float, default=500.0)
    parser.add_argument("--oscillatory-t-end", type=float, default=2000.0)
    parser.add_argument("--sample-every", type=int, default=10)
    parser.add_argument("--no-plots", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    validate_positive("dt", args.dt)
    validate_positive("equilibrium_t_end", args.equilibrium_t_end)
    validate_positive("oscillatory_t_end", args.oscillatory_t_end)
    if args.oscillatory_t_end < 2000.0:
        raise ValueError("oscillatory_t_end must be at least 2000 for the fixed damping windows")
    if args.sample_every <= 0:
        raise ValueError("sample_every must be positive")

    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    parameter_rows = read_scenarios(input_path)
    summaries = []
    plot_data = {}

    print(f"Writing main simulation outputs to {output_dir}")
    for spec in SCENARIOS:
        raw_params = parameter_rows[spec.name]
        params = params_from_dict(raw_params)
        t_end = args.oscillatory_t_end if spec.oscillatory else args.equilibrium_t_end
        print(f"  {spec.name}: RK4 dt={args.dt:g}, t_end={t_end:g}")
        times, states, solver_info = rk4(
            lambda state, current=params: rhs(state, current),
            list(args.init),
            t_end=t_end,
            dt=args.dt,
        )
        trajectory_path = output_dir / f"{spec.name}_rk4.csv"
        write_trajectory(trajectory_path, times, states, args.sample_every)
        diagnostics, summary = analyze_scenario(spec, params, times, states, solver_info)
        with (output_dir / f"{spec.name}_diagnostics.json").open("w", encoding="utf-8") as handle:
            json.dump(diagnostics, handle, indent=2)
        summaries.append(summary)
        plot_data[spec.name] = {
            "spec": spec,
            "times": times,
            "states": states,
            "target": list(diagnostics["target_equilibrium"].values()),
        }

    write_summary(output_dir / "scenario_summary.csv", summaries)
    manifest = {
        "method": "rk4",
        "dt": args.dt,
        "initial_state": dict(zip(COMPONENTS, args.init)),
        "equilibrium_t_end": args.equilibrium_t_end,
        "oscillatory_t_end": args.oscillatory_t_end,
        "sample_every": args.sample_every,
        "input_csv": os.path.relpath(input_path, PROJECT_ROOT),
        "scenarios": [
            {**asdict(spec), "parameters": parameter_rows[spec.name]}
            for spec in SCENARIOS
        ],
    }
    with (output_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    if not args.no_plots:
        from .plots import create_all_plots

        create_all_plots(plot_data, output_dir)

    print("Completed all five main scenarios.")


if __name__ == "__main__":
    main()
