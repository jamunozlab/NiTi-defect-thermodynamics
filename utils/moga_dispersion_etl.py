#!/usr/bin/env python3

import os
import re
import argparse
import subprocess
import numpy as np
import yaml


def read_inputc(input_path):
    params = {}

    with open(input_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if "=" in line:
                key, value = line.split("=", 1)
            else:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                key, value = parts

            params[key.strip()] = value.strip()

    a_val = float(params["a_val"])

    atomic_masses = np.array([
        float(x)
        for x in re.split(r"[,\s]+", params["atomic_masses"])
        if x
    ])

    return a_val, atomic_masses


def read_best_solutions(generation_output_path):
    best_solutions = []

    with open(generation_output_path, "r") as f:
        for line in f:
            if "Best solution" not in line:
                continue

            match = re.search(r"\[([^\]]+)\]", line)
            if match is None:
                continue

            solution = [
                float(x)
                for x in re.split(r"[,\s]+", match.group(1).strip())
                if x
            ]

            best_solutions.append(solution)

    return np.array(best_solutions)


def read_band_yaml(band_yaml_path):
    with open(band_yaml_path, "r") as f:
        data = yaml.safe_load(f)

    q_distance = []
    frequencies = []

    for qpoint in data["phonon"]:
        q_distance.append(qpoint["distance"])
        frequencies.append([
            band["frequency"]
            for band in qpoint["band"]
        ])

    return np.array(q_distance), np.array(frequencies)


def write_solution_file(solution, output_path):
    np.savetxt(output_path, np.asarray(solution)[None, :])


def generate_band_yaml_for_solution(
    simulation_path,
    solution,
    generation_index,
    dispersion_script=None,
    solution_filename="current_solution.txt",
    verbose=False,
):
    simulation_path = os.path.abspath(simulation_path)

    if dispersion_script is None:
        repo_root = os.path.abspath(os.path.join(simulation_path, "..", ".."))
        dispersion_script = os.path.join(
            repo_root,
            "executables",
            "dispersion_generator_from_solution.py",
        )

    if not os.path.isfile(dispersion_script):
        raise FileNotFoundError(f"Missing dispersion script: {dispersion_script}")

    solution_path = os.path.join(simulation_path, solution_filename)
    write_solution_file(solution, solution_path)

    command = [
        "python",
        dispersion_script,
        "--simulation_path", simulation_path,
        "--solution_file", solution_path,
        "--generation", str(generation_index),
    ]

    if verbose:
        print("Running:", " ".join(command))

    subprocess.run(command, cwd=simulation_path, check=True)

    generation_tag = f"{generation_index:06d}"
    band_yaml_path = os.path.join(
        simulation_path,
        f"band_generation_{generation_tag}.yaml",
    )

    if not os.path.isfile(band_yaml_path):
        raise FileNotFoundError(f"Expected band YAML not found: {band_yaml_path}")

    return band_yaml_path


def etl_simulation_dispersion_data(
    simulation_path,
    input_filename="inputc",
    generation_output_filename="generation_output_mod.txt",
    output_filename="moga_dispersion_data.npz",
    dispersion_script=None,
    verbose=False,
):
    simulation_path = os.path.abspath(simulation_path)

    input_path = os.path.join(simulation_path, input_filename)
    generation_output_path = os.path.join(
        simulation_path,
        generation_output_filename,
    )
    output_path = os.path.join(simulation_path, output_filename)

    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Missing input file: {input_path}")

    if not os.path.isfile(generation_output_path):
        raise FileNotFoundError(
            f"Missing generation output file: {generation_output_path}"
        )

    a_val, atomic_masses = read_inputc(input_path)
    best_solutions = read_best_solutions(generation_output_path)

    if best_solutions.size == 0:
        raise ValueError(f"No best solutions found in {generation_output_path}")

    all_q_distances = []
    all_frequencies = []

    for gen_idx, solution in enumerate(best_solutions):
        if verbose:
            print(f"Generation {gen_idx}: solution = {solution}")

        band_yaml_path = generate_band_yaml_for_solution(
            simulation_path=simulation_path,
            solution=solution,
            generation_index=gen_idx,
            dispersion_script=dispersion_script,
            verbose=verbose,
        )

        q_distance, frequencies = read_band_yaml(band_yaml_path)

        all_q_distances.append(q_distance)
        all_frequencies.append(frequencies)

    all_q_distances = np.array(all_q_distances)
    all_frequencies = np.array(all_frequencies)

    np.savez_compressed(
        output_path,
        a_val=a_val,
        atomic_masses=atomic_masses,
        generations=np.arange(len(best_solutions)),
        best_solutions=best_solutions,
        q_distances=all_q_distances,
        frequencies=all_frequencies,
    )

    if verbose:
        print(f"Wrote {output_path}")
        print(f"  a_val = {a_val}")
        print(f"  atomic_masses = {atomic_masses}")
        print(f"  best_solutions shape = {best_solutions.shape}")
        print(f"  q_distances shape = {all_q_distances.shape}")
        print(f"  frequencies shape = {all_frequencies.shape}")

    return output_path


def read_experiment_file(experiment_path):
    simulation_names = []

    with open(experiment_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("simulation"):
                simulation_names.append(line)

    return simulation_names


def etl_experiment_dispersion_data(
    experiment_name,
    repo_root=None,
    input_filename="inputc",
    generation_output_filename="generation_output_mod.txt",
    output_filename="moga_dispersion_data.npz",
    dispersion_script=None,
    verbose=False,
):
    user = os.environ["USER"]

    if repo_root is None:
        repo_root = f"/home/{user}/MOGA-Phonons"

    repo_root = os.path.abspath(repo_root)

    experiment_path = os.path.join(repo_root, "experiments", experiment_name)
    simulations_root = os.path.join(repo_root, "simulations")

    if not os.path.isfile(experiment_path):
        raise FileNotFoundError(f"Experiment file not found: {experiment_path}")

    simulation_names = read_experiment_file(experiment_path)

    if verbose:
        print(f"Experiment: {experiment_name}")
        print(f"Repository root: {repo_root}")
        print(f"Found {len(simulation_names)} simulations.")

    output_paths = []

    for sim_name in simulation_names:
        simulation_path = os.path.join(simulations_root, sim_name)

        try:
            output_path = etl_simulation_dispersion_data(
                simulation_path=simulation_path,
                input_filename=input_filename,
                generation_output_filename=generation_output_filename,
                output_filename=output_filename,
                dispersion_script=dispersion_script,
                verbose=verbose,
            )
            output_paths.append(output_path)

        except Exception as err:
            print(f"WARNING: Skipping {sim_name}")
            print(f"  {err}")

    print(f"\nDone. Wrote {len(output_paths)} dispersion npz files.\n")

    return output_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Generate q_distance and frequency arrays for all MOGA best "
            "solutions in an experiment."
        )
    )

    parser.add_argument(
        "experiment_name",
        help="Experiment filename, e.g. experiment000002",
    )

    parser.add_argument(
        "--repo_root",
        default=None,
        help="Root path of MOGA-Phonons. Default: /home/$USER/MOGA-Phonons",
    )

    parser.add_argument(
        "--input_filename",
        default="inputc",
        help="Input filename inside each simulation directory.",
    )

    parser.add_argument(
        "--generation_output_filename",
        default="generation_output_mod.txt",
        help="Generation output filename inside each simulation directory.",
    )

    parser.add_argument(
        "--output_filename",
        default="moga_dispersion_data.npz",
        help="Output npz filename written inside each simulation directory.",
    )

    parser.add_argument(
        "--dispersion_script",
        default=None,
        help=(
            "Path to dispersion_generator_from_solution.py. "
            "Default: <repo_root>/executables/dispersion_generator_from_solution.py"
        ),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print detailed progress.",
    )

    args = parser.parse_args()

    etl_experiment_dispersion_data(
        experiment_name=args.experiment_name,
        repo_root=args.repo_root,
        input_filename=args.input_filename,
        generation_output_filename=args.generation_output_filename,
        output_filename=args.output_filename,
        dispersion_script=args.dispersion_script,
        verbose=args.verbose,
    )
