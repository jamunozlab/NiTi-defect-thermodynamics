#!/usr/bin/env python3

import argparse
from pathlib import Path

from utils.lammps_etl_utils import etl


def read_experiment_file(experiment_path):
    simulation_names = []

    with open(experiment_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("simulation"):
                simulation_names.append(line)

    return simulation_names


def etl_experiment(
    experiment_name,
    repo_root=None,
    input_filename="inputc",
    output_filename="log.lammps",
    verbose=False,
):
    # Automatically detect repo root from this file location
    if repo_root is None:
        repo_root = Path(__file__).resolve().parent
    else:
        repo_root = Path(repo_root)

    experiment_path = repo_root / "experiments" / experiment_name
    simulations_root = repo_root / "simulations"

    if not experiment_path.is_file():
        raise FileNotFoundError(f"Experiment file not found: {experiment_path}")

    simulation_names = read_experiment_file(experiment_path)

    if verbose:
        print(f"Experiment: {experiment_name}")
        print(f"Repo root: {repo_root}")
        print(f"Found {len(simulation_names)} simulations.")

    output_paths = []

    for sim_name in simulation_names:
        simulation_path = simulations_root / sim_name

        try:
            output_path = etl(
                simulation_path=str(simulation_path),
                input_filename=input_filename,
                output_filename=output_filename,
                verbose=verbose,
            )
            output_paths.append(output_path)

        except (FileNotFoundError, ValueError) as err:
            print(f"WARNING: Skipping {sim_name}")
            print(f"  {err}")

    print(f"\nDone. Wrote {len(output_paths)} npz files.\n")

    return output_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create compressed thermo.npz summaries for all simulations in a LAMMPS experiment."
    )

    parser.add_argument(
        "experiment_name",
        type=str,
        help="Experiment filename, e.g. experiment000000",
    )

    parser.add_argument(
        "--repo_root",
        type=str,
        default=None,
        help="Root path of the NiTi repository. Default: directory containing this script.",
    )

    parser.add_argument(
        "--input_filename",
        type=str,
        default="inputc",
        help="Input filename inside each simulation directory (retained for ETL compatibility).",
    )

    parser.add_argument(
        "--output_filename",
        type=str,
        default="log.lammps",
        help="LAMMPS log filename inside each simulation directory.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print details while processing.",
    )

    args = parser.parse_args()

    etl_experiment(
        experiment_name=args.experiment_name,
        repo_root=args.repo_root,
        input_filename=args.input_filename,
        output_filename=args.output_filename,
        verbose=args.verbose,
    )
