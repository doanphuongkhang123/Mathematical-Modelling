"""CSV input/output and console formatting for equilibrium-analysis outputs."""

from __future__ import annotations

import csv
import os

from analysis import compute_all_equilibria, p0_stability
from model import PARAMETER_NAMES, State, params_from_dict


ScenarioRows = list[tuple[str, dict[str, float]]]


def read_scenario_csv(path: str) -> ScenarioRows:
    examples: ScenarioRows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            examples.append(
                (
                    row["scenario"],
                    {name: float(row[name]) for name in PARAMETER_NAMES},
                )
            )
    return examples


def write_equilibrium_csv(examples: ScenarioRows, output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("scenario", "equilibrium", "exists", "S", "x", "y", "z", "P0_stable", "f1(1)-D1"))
        for scenario_name, raw_params in examples:
            params = params_from_dict(raw_params)
            p0_is_stable, threshold = p0_stability(params)
            for equilibrium_name, state in compute_all_equilibria(params).items():
                if state is None:
                    writer.writerow((scenario_name, equilibrium_name, "no", "", "", "", "", p0_is_stable, threshold))
                else:
                    writer.writerow((scenario_name, equilibrium_name, "yes", *state, p0_is_stable, threshold))


def format_state(state: State | None) -> str:
    if state is None:
        return "does not exist"
    return "(" + ", ".join(f"{value:.8g}" for value in state) + ")"


def print_summary(examples: ScenarioRows) -> None:
    for scenario_name, raw_params in examples:
        params = params_from_dict(raw_params)
        p0_is_stable, threshold = p0_stability(params)
        print(f"\n{scenario_name}")
        print(f"  P0 threshold f1(1)-D1 = {threshold:.8g}; P0 stable? {p0_is_stable}")
        for equilibrium_name, state in compute_all_equilibria(params).items():
            print(f"  {equilibrium_name}: {format_state(state)}")
