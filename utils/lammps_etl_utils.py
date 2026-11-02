#!/usr/bin/env python
"""Extract LAMMPS thermodynamic output into a compressed NumPy archive."""

from pathlib import Path
import numpy as np


def get_thermo_scalars_from_lammps_log_file(scalar, simulation_path, output_filename='log.lammps'):
    """Return a NumPy array for a requested LAMMPS thermo column.

    Reads the thermo table associated with the last ``run`` command, excluding
    preceding minimization tables. Supported scalars: Step, Temp, TotEng, Press.
    """
    columns = ('Step', 'Temp', 'TotEng', 'Press')
    if scalar not in columns:
        raise ValueError(f"scalar must be one of {columns}; got {scalar!r}")

    log_path = Path(simulation_path) / output_filename
    tables = []
    current_rows = None
    indices = None
    in_run = False

    with log_path.open('r') as log_file:
        for line in log_file:
            stripped = line.strip()
            parts = stripped.split()

            if parts and parts[0] == 'run' and len(parts) >= 2:
                in_run = True
                current_rows = None
                indices = None
                continue

            if in_run and all(name in parts for name in columns):
                indices = [parts.index(name) for name in columns]
                current_rows = []
                tables.append(current_rows)
                continue

            if current_rows is None or indices is None:
                continue

            if stripped.startswith('Loop time of'):
                current_rows = None
                indices = None
                in_run = False
                continue

            try:
                values = [float(parts[i]) for i in indices]
            except (ValueError, IndexError):
                # Non-data lines terminate the current thermo table.
                current_rows = None
                indices = None
                continue

            current_rows.append(values)

    if not tables or not tables[-1]:
        raise ValueError(f"No LAMMPS run thermo table with {columns} found in {log_path}")

    data = np.asarray(tables[-1], dtype=float)
    result = data[:, columns.index(scalar)]
    return result.astype(np.int64) if scalar == 'Step' else result


def etl(simulation_path, input_filename='inputc', output_filename='log.lammps', verbose=False):
    """Save the last LAMMPS run's Step, Temp, TotEng, Press as thermo.npz.

    input_filename is retained for compatibility with the AMReX ETL signature.
    """
    scalars = ('Step', 'Temp', 'TotEng', 'Press')
    data = {
        name: get_thermo_scalars_from_lammps_log_file(
            name, simulation_path, output_filename=output_filename
        )
        for name in scalars
    }
    output_path = Path(simulation_path) / 'thermo.npz'
    np.savez_compressed(output_path, **data)
    if verbose:
        print(f"Saved {len(data['Step'])} thermo records to {output_path}")
    return output_path
