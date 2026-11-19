#!/usr/bin/env python3
"""Build one Birch-Murnaghan summary row per LAMMPS experiment.

Examples:
    python build_dataframe.py --experiment experiment000012
    python build_dataframe.py --experiment experiment000011 experiment000012

Read nominal settings from simulation_parameters.json where available. For older
simulations, fall back to inputc for temperature and atoms_positions.data for
composition. Never infer antisite concentrations from composition alone.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from models.bm_eos import LammpsExperiment


def read_experiment_file(path):
    names = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip().startswith('simulation')]
    if not names:
        raise ValueError(f'No simulations in {path}')
    return names


def read_simulation_metadata(simulation_dir):
    result = {'temperature_K': np.nan, 'n_atoms': np.nan, 'n_Ni': np.nan,
              'n_Ti': np.nan, 'x_Ni': np.nan, 'lattice_parameter_A': np.nan,
              'system_size': np.nan, 'ntimesteps': np.nan, 'random_seed': np.nan,
              'reduction_percentage': np.nan, 'atom_type_to_reduce': np.nan,
              'defect_percentage': np.nan, 'atom_type_to_remove': np.nan,
              'vacancy_percentage': np.nan, 'metadata_source': 'legacy'}
    json_path = simulation_dir / 'simulation_parameters.json'
    if json_path.is_file():
        with json_path.open() as f:
            params = json.load(f)
        if not isinstance(params, dict):
            raise ValueError(f'Metadata is not an object: {json_path}')
        for field in ('system_size', 'ntimesteps', 'random_seed',
                      'reduction_percentage', 'atom_type_to_reduce',
                      'defect_percentage', 'atom_type_to_remove', 'vacancy_percentage'):
            if field in params and params[field] is not None:
                result[field] = params[field]
        if 'temperature' in params:
            result['temperature_K'] = float(params['temperature'])
        if 'lattice_parameter' in params:
            result['lattice_parameter_A'] = float(params['lattice_parameter'])
        result['metadata_source'] = 'json'
    else:
        input_path = simulation_dir / 'inputc'
        if input_path.is_file():
            match = re.search(r'^\s*fix\s+\S+\s+all\s+nvt\s+temp\s+([\d.eE+-]+)',
                              input_path.read_text(), re.M)
            if match:
                result['temperature_K'] = float(match.group(1))

    data_path = simulation_dir / 'atoms_positions.data'
    if not data_path.is_file():
        return result
    lines = data_path.read_text().splitlines()
    for line in lines:
        match = re.match(r'^\s*(\d+)\s+atoms\s*$', line)
        if match:
            result['n_atoms'] = int(match.group(1))
            break
    start = next((i for i, line in enumerate(lines)
                  if line.strip().startswith('Atoms')), None)
    if start is None:
        return result
    counts = {1: 0, 2: 0}
    for line in lines[start + 1:]:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            int(parts[0])
            species = int(parts[1])
            for value in parts[2:5]:
                float(value)
        except ValueError:
            continue
        if species in counts:
            counts[species] += 1
    total = sum(counts.values())
    if total:
        result.update(n_atoms=total, n_Ni=counts[1], n_Ti=counts[2],
                      x_Ni=counts[1] / total)
    return result


def constant_or_nan(metadata, key, *, atol=1e-8):
    """Only report a scalar if all simulations supply a consistent value."""
    values = [m[key] for m in metadata]
    if any(pd.isna(v) for v in values):
        return np.nan
    try:
        numeric = np.asarray(values, dtype=float)
        if np.allclose(numeric, numeric[0], rtol=0, atol=atol):
            return float(numeric[0])
    except (ValueError, TypeError):
        if all(v == values[0] for v in values):
            return values[0]
    return np.nan


def summarize_metadata(repo_root, simulation_names):
    metadata = [read_simulation_metadata(repo_root / 'simulations' / name)
                for name in simulation_names]
    row = {'n_simulations': len(simulation_names),
           'simulation_ids': simulation_names,
           'metadata_complete': all(m['metadata_source'] == 'json' for m in metadata),
           'n_metadata_json': sum(m['metadata_source'] == 'json' for m in metadata)}
    for field in ('temperature_K', 'x_Ni', 'n_atoms', 'n_Ni', 'n_Ti',
                  'system_size', 'ntimesteps', 'random_seed',
                  'reduction_percentage', 'atom_type_to_reduce',
                  'defect_percentage', 'atom_type_to_remove', 'vacancy_percentage'):
        row[field] = constant_or_nan(metadata, field)
    lattice = [m['lattice_parameter_A'] for m in metadata]
    known = [float(v) for v in lattice if pd.notna(v)]
    row['lattice_parameter_min_A'] = min(known) if len(known) == len(metadata) else np.nan
    row['lattice_parameter_max_A'] = max(known) if len(known) == len(metadata) else np.nan
    row['n_lattice_parameters'] = len(set(known)) if len(known) == len(metadata) else np.nan

    # These are *nominal* input settings; an integer number of sites is changed.
    antisite = row['defect_percentage']
    vacancy = row['vacancy_percentage']
    if pd.isna(antisite) or pd.isna(vacancy):
        row['defect_type'] = 'unknown'
        row['defect_concentration_percent'] = np.nan
    elif antisite > 0 and vacancy > 0:
        row['defect_type'] = 'antisite+vacancy'
        row['defect_concentration_percent'] = np.nan
    elif antisite > 0:
        row['defect_type'] = 'antisite'
        row['defect_concentration_percent'] = antisite
    elif vacancy > 0:
        row['defect_type'] = 'vacancy'
        row['defect_concentration_percent'] = vacancy
    else:
        row['defect_type'] = 'none'
        row['defect_concentration_percent'] = 0.0
    return row


def build_dataframe(experiments, repo_root=None, discard_fraction=0.5):
    root = Path(repo_root).resolve() if repo_root else Path(__file__).resolve().parent
    rows = []
    for experiment_name in experiments:
        name = Path(experiment_name).name
        row = {'experiment': name, 'fit_status': 'failed', 'fit_error': None,
               'discard_fraction': discard_fraction}
        try:
            names = read_experiment_file(root / 'experiments' / name)
            row.update(summarize_metadata(root, names))
            experiment = LammpsExperiment(name, repo_root=root)
            fit = experiment.fit_birch_murnaghan(discard_fraction=discard_fraction)
            row.update({key: float(fit[key]) for key in (
                'E0_eV_per_atom', 'V0_A3_per_atom', 'a0_A', 'B0_eV_per_A3',
                'B0_GPa', 'B0_prime')})
            row['fit_status'] = 'ok'
        except (FileNotFoundError, ValueError, RuntimeError, TypeError,
                KeyError, json.JSONDecodeError) as exc:
            row['fit_error'] = str(exc)
            print(f'WARNING: {name}: {exc}')
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', nargs='+', required=True)
    parser.add_argument('--repo_root', default=None)
    parser.add_argument('--output_filename', default=None)
    parser.add_argument('--discard_fraction', type=float, default=0.5)
    args = parser.parse_args()
    root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parent
    df = build_dataframe(args.experiment, root, args.discard_fraction)
    if len(args.experiment) == 1:
        match = re.search(r'(\d+)$', Path(args.experiment[0]).name)
        default_name = f'dataframe{match.group(1)}.pkl' if match else 'bm_experiments.pkl'
    else:
        default_name = 'bm_experiments.pkl'
    output_path = root / 'dataframes' / (args.output_filename or default_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(output_path)
    print(f'Saved {len(df)} experiment rows to {output_path}')
    print(df[['experiment', 'temperature_K', 'x_Ni', 'defect_type',
              'defect_concentration_percent', 'fit_status']].to_string(index=False))


if __name__ == '__main__':
    main()
