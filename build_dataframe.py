#!/usr/bin/env python3
"""Build one Birch–Murnaghan summary row per LAMMPS experiment.

Examples (run from repository root):
    python build_dataframe.py --experiment experiment000011
    python build_dataframe.py --experiment experiment000011 experiment000012
    python build_dataframe.py --experiment experiment000011 --output_filename bm_results.pkl

Uses models/bm_eos.py without modifying it. The current EOS class supports
pristine B2 only; failed/defective fits are recorded with fit_status/error.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from models.bm_eos import LammpsExperiment


def read_experiment_file(experiment_path):
    lines = Path(experiment_path).read_text().splitlines()
    names = [line.strip() for line in lines if line.strip().startswith("simulation")]
    if not names:
        raise ValueError(f"No simulations in {experiment_path}")
    return names


def read_simulation_metadata(simulation_dir):
    """Read measurable metadata; don't guess nominal defect settings."""
    input_path = simulation_dir / "inputc"
    data_path = simulation_dir / "atoms_positions.data"
    result = {
        "temperature_K": np.nan,
        "n_atoms": np.nan,
        "n_Ni": np.nan,
        "n_Ti": np.nan,
        "x_Ni": np.nan,
        "defect_type": "unknown",
        "defect_concentration_percent": np.nan,
    }

    if input_path.is_file():
        text = input_path.read_text()
        match = re.search(r"^\s*fix\s+\S+\s+all\s+nvt\s+temp\s+([\d.eE+-]+)", text, re.M)
        if match:
            result["temperature_K"] = float(match.group(1))

    if not data_path.is_file():
        return result

    lines = data_path.read_text().splitlines()
    for line in lines:
        match = re.match(r"^\s*(\d+)\s+atoms\s*$", line)
        if match:
            result["n_atoms"] = int(match.group(1))
            break

    start = next((i for i, line in enumerate(lines) if line.strip().startswith("Atoms")), None)
    if start is None:
        return result

    counts = {1: 0, 2: 0}
    for line in lines[start + 1:]:
        parts = line.split()
        if len(parts) >= 5:
            try:
                int(parts[0]); species = int(parts[1]); float(parts[2]); float(parts[3]); float(parts[4])
            except ValueError:
                continue
            if species in counts:
                counts[species] += 1
    total = counts[1] + counts[2]
    if total:
        result.update(n_atoms=total, n_Ni=counts[1], n_Ti=counts[2], x_Ni=counts[1] / total)
    return result


def summarize_metadata(repo_root, simulation_names):
    metadata = [read_simulation_metadata(repo_root / "simulations" / name)
                for name in simulation_names]
    row = {"n_simulations": len(simulation_names),
           "simulation_ids": simulation_names,
           "temperature_K": np.nan, "x_Ni": np.nan,
           "defect_type": "unknown", "defect_concentration_percent": np.nan}
    for key in ("temperature_K", "x_Ni"):
        values = [m[key] for m in metadata if pd.notna(m[key])]
        if len(values) == len(metadata) and np.allclose(values, values[0], atol=1e-8):
            row[key] = float(values[0])
    # An equiatomic composition alone does not establish absence of antisites.
    # Nominal defect type/concentration are not saved in inputc or thermo.npz.
    return row


def build_dataframe(experiments, repo_root=None, discard_fraction=0.5):
    repo_root = Path(repo_root).resolve() if repo_root else Path(__file__).resolve().parent
    rows = []
    for experiment_name in experiments:
        name = Path(experiment_name).name
        row = {"experiment": name, "fit_status": "failed", "fit_error": None}
        try:
            simulation_names = read_experiment_file(repo_root / "experiments" / name)
            row.update(summarize_metadata(repo_root, simulation_names))
            experiment = LammpsExperiment(name, repo_root=repo_root)
            fit = experiment.fit_birch_murnaghan(discard_fraction=discard_fraction)
            row.update({key: float(fit[key]) for key in (
                "E0_eV_per_atom", "V0_A3_per_atom", "a0_A", "B0_eV_per_A3",
                "B0_GPa", "B0_prime")})
            row["discard_fraction"] = discard_fraction
            row["fit_status"] = "ok"
        except (FileNotFoundError, ValueError, RuntimeError, TypeError) as exc:
            row["fit_error"] = str(exc)
            print(f"WARNING: {name}: {exc}")
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", nargs="+", required=True,
                        help="One or more experiment manifest names")
    parser.add_argument("--repo_root", default=None)
    parser.add_argument("--output_filename", default=None)
    parser.add_argument("--discard_fraction", type=float, default=0.5)
    args = parser.parse_args()
    root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parent
    df = build_dataframe(args.experiment, root, args.discard_fraction)
    if len(args.experiment) == 1:
        match = re.search(r"(\d+)$", args.experiment[0])
        default_output_name = (
            f"dataframe{match.group(1)}.pkl"
            if match
            else "bm_experiments.pkl"
    )
    else:
        default_output_name = "bm_experiments.pkl"

    output_name = args.output_filename or default_output_name
    output_path = root / "dataframes" / output_name
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(output_path)
    print(f"Saved {len(df)} experiment rows to {output_path}")
    print(df[["experiment", "temperature_K", "x_Ni", "fit_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
