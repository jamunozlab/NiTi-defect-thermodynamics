# NiTi Defect Thermodynamics

This repository contains a lightweight workflow for generating, running, extracting, and analyzing molecular dynamics (MD) simulations of B2 NiTi using LAMMPS and a modified embedded-atom method (MEAM) potential. It uses **UTEP (Utilities To Execute Pipelines)** to organize parameter sweeps into reproducible simulation directories and experiment manifests.

The workflow has four layers:

1. **Experiment generation:** create B2 NiTi structures and LAMMPS input files for a parameter sweep.
2. **Execution:** run the simulations using LAMMPS and MPI.
3. **ETL (extract, transform, load):** convert LAMMPS thermodynamic logs to compressed NumPy archives.
4. **Analysis:** load all simulation archives for an experiment into a pandas DataFrame and visualize thermodynamic time series.

## Requirements

- Python 3 with NumPy
- pandas and Matplotlib for the analysis notebook
- Jupyter Notebook or JupyterLab to run the analysis notebook
- LAMMPS compiled with MEAM support
- An MPI implementation (such as Open MPI) for the current local launcher
- A Unix-like shell (the current workflow has been tested on macOS)

## Repository layout

```text
NiTi-defect-thermodynamics/
├── new_experiment.py                    # Define and launch parameter sweeps
├── etl_experiment_lammps.py             # Extract data for a complete experiment
├── UTEP_LAMMPS_experiment_analysis.ipynb # DataFrame analysis and plotting
├── utils/
│   ├── config.py                        # Local paths, executable, linked files
│   ├── default_dictionaries.py          # Default simulation parameters
│   ├── new_simulation_from_default.py
│   ├── write_input_from_default.py      # LAMMPS input and B2 structure generation
│   ├── write_run_sbatch_from_default.py
│   └── lammps_etl_utils.py              # Extract LAMMPS thermodynamic time series
├── executables/
│   ├── NiTi.library.meam                # MEAM library parameters
│   ├── NiTi.meam                        # MEAM alloy parameters
│   └── lmp_mpi                          # Local symbolic link, NOT a bundled binary
├── simulations/
│   └── .gitkeep
└── experiments/
    └── .gitkeep
```

## LAMMPS executable: important symlink setup

The `executables/lmp_mpi` entry is a **symbolic link**, not a copy of LAMMPS. Git records the link target; it **does not upload the LAMMPS binary**. In the original macOS setup, the link points to `/opt/homebrew/bin/lmp_mpi`. That path may not exist on another machine, so **recreate the link locally after cloning**.

First, install LAMMPS with MEAM support and an MPI implementation. Check where the executable is installed:

```bash
command -v lmp_mpi
```

From the repository root, replace the tracked link with one pointing to your installed executable:

```bash
rm executables/lmp_mpi
ln -s "$(command -v lmp_mpi)" executables/lmp_mpi
```

Verify the link and executable:

```bash
ls -l executables/lmp_mpi
./executables/lmp_mpi -h
```

The LAMMPS help output should list MEAM support. For an MPI smoke test (adjust the process count to your machine):

```bash
mpirun -np 4 ./executables/lmp_mpi -h > mpi_test.log 2>&1
```

**Note:** Replacing the tracked symlink changes your local Git working tree. This is expected for a machine-specific path; do not commit that local symlink change unless you intentionally want to update the repository's default link. Alternatively, configure `executable` and `executables_path` in `utils/config.py` to reference a suitable local executable directory.

## Configure the workflow

Edit `utils/config.py` to set the local paths for `simulations/`, `experiments/`, and `executables/`. The current configuration uses absolute paths, so update them after cloning. Ensure that the executable name matches the file or symlink in `executables/` and that `other_files` includes the two MEAM parameter files:

```python
other_files = ["NiTi.library.meam", "NiTi.meam"]
```

Do **not** include `atoms_positions.data` in `other_files`: it is generated separately for each simulation.

Simulation defaults are defined in the `lammps` dictionary in `utils/default_dictionaries.py`. The current options include:

| Parameter | Meaning |
| --- | --- |
| `temperature` | MD target temperature (K) |
| `timestep` | MD timestep (ps with LAMMPS `metal` units) |
| `ntimesteps` | Number of MD integration steps |
| `logstep` | Interval between thermodynamic log records |
| `tdamp` | Thermostat damping time (ps) |
| `system_size` | Number of B2 conventional cells along each box direction |
| `lattice_parameter` | Cubic B2 lattice parameter (Å) |
| `atom_type_to_reduce` | Species converted to the other species (1 = Ni, 2 = Ti) |
| `reduction_percentage` | Fraction of all sites converted (e.g., `0.01` = 1% of sites) |
| `defect_percentage` | Nominal antisite percentage (e.g., `1.0` = 1%) |
| `atom_type_to_remove` | Species selected for vacancies (1 = Ni, 2 = Ti) |
| `vacancy_percentage` | Percentage of the selected species removed (e.g., `1.0` = 1%) |
| `random_seed` | Seed for reproducible random selections |

The B2 structure generator creates two sites per cubic unit cell, with ideal sublattice positions `(0, 0, 0)` and `(1/2, 1/2, 1/2)`. It applies composition changes, antisite assignments, and vacancy removal in that order, then writes `atoms_positions.data` in each simulation directory.

## 1. Define and run an experiment

Edit the section between `#### SCRIPT SHOULD BE BELOW ####` and `#### SCRIPT SHOULD BE ABOVE ####` in `new_experiment.py`. For example, to sweep temperature and lattice parameter, use nested loops and pass the overrides to `new_simulation_from_default()`:

```python
experiment_description = (
    "NiTi temperature and lattice-parameter sweep\n"
    "Temperatures: 200, 300, 400 K\n"
    "Lattice parameters: 2.95, 3.00, 3.05 Angstrom\n"
)

for lattice_parameter in [2.95, 3.00, 3.05]:
    for temperature in [200, 300, 400]:
        sim_params_dict = {
            "temperature": temperature,
            "lattice_parameter": lattice_parameter,
        }
        sbatch_options_dict = {"nproc": 4}

        simulation_id = new_simulation_from_default(
            default_input_dictionary,
            simulations_path,
            executables_path,
            executable,
            sim_params_dict,
            sbatch_options_dict,
            other_files=other_files,
            submit_job=submit_job,
        )
        simulation_ids.append(simulation_id)

experiment_id = find_next_available_file("experiment", experiments_path)
```

From the repository root, create simulation directories **without running LAMMPS**:

```bash
python new_experiment.py lammps --test
```

To create and execute a fresh set of simulations:

```bash
python new_experiment.py lammps
```

The second command creates **new** simulation directories; it does not resume the test directories. The current macOS runner launches simulations sequentially, with each simulation using the configured number of MPI processes.

Each simulation directory contains its generated LAMMPS input (`inputc`), structure (`atoms_positions.data`), run script (`runc.sh`), linked MEAM files, and output files after execution. The `experiments/` directory contains a manifest listing the simulation IDs belonging to each sweep.

## 2. Extract LAMMPS thermodynamic data (ETL)

After running an experiment, use `etl_experiment_lammps.py` from the repository root. For example:

```bash
python etl_experiment_lammps.py experiment000009
```

For detailed processing messages:

```bash
python etl_experiment_lammps.py experiment000009 --verbose
```

The experiment-level script reads `experiments/experiment000009`, identifies the associated simulation directories, and calls `etl()` from `utils/lammps_etl_utils.py` for each simulation. It reads the thermodynamic table associated with the **last LAMMPS `run` command** in each `log.lammps` file, rather than the preceding minimization table.

Each successfully processed simulation receives **one** compressed archive:

```text
simulations/simulationXXXXXX/thermo.npz
```

The archive contains four NumPy arrays, each with one value per recorded timestep:

| Array key | Quantity | Units (`metal`) |
| --- | --- | --- |
| `Step` | MD timestep index | dimensionless |
| `Temp` | Temperature | K |
| `TotEng` | Total energy | eV |
| `Press` | Pressure | bar |

The ETL script identifies columns from the log header, so their positions in the LAMMPS thermo table need not be fixed. If a simulation's log is missing or has no suitable thermodynamic table, the experiment-level driver prints a warning and continues with the remaining simulations. Re-running ETL overwrites an existing `thermo.npz` for a successfully processed simulation.

**Command-line options:**

```bash
python etl_experiment_lammps.py --help
```

`--repo_root` specifies a different repository location; `--input_filename` retains the compatible ETL argument (currently unused by the log parser); and `--output_filename` specifies the **input LAMMPS log filename**, defaulting to `log.lammps`. The compressed output filename remains `thermo.npz`.

You can also inspect one archive directly:

```python
import numpy as np

with np.load("simulations/simulationXXXXXX/thermo.npz") as data:
    for key in data.files:
        print(key, data[key].shape)
```

## 3. Analyze an experiment in Jupyter

Open `UTEP_LAMMPS_experiment_analysis.ipynb` with Jupyter Notebook or JupyterLab. Run it with the working directory set to the repository root, or change `REPO_ROOT` to an absolute path.

The notebook defines a lightweight class, `LammpsExperiment`, that reads the experiment manifest, loads each simulation's `thermo.npz`, and combines the time series into a single pandas DataFrame. Each row represents one recorded timestep in one simulation, with columns:

```text
simulation_id    Step    Temp    TotEng    Press
```

To load an experiment:

```python
from pathlib import Path

experiment = LammpsExperiment(
    "experiment000009",
    repo_root=Path.cwd(),
)

print(experiment.df.head())
```

To plot temperature, total energy, or pressure against timestep for **all** simulations in the experiment:

```python
experiment.plot("Temp")
experiment.plot("TotEng")
experiment.plot("Press")
```

To plot only a selected simulation (replace the identifier with one listed in your experiment manifest):

```python
experiment.plot("Temp", simulation_ids="simulationXXXXXX")
```

The `plot()` method also accepts a list of simulation IDs or an existing Matplotlib `ax`. The notebook skips missing or inconsistent archives with a warning; it raises an error if no valid archives are available.

**Note:** `LammpsExperiment` is defined in the notebook, not in an independently importable Python module. Execute the class-definition cell before using the examples above.

## Git and generated data

The repository tracks `.gitkeep` placeholders so `simulations/` and `experiments/` exist after cloning. The `.gitignore` rules exclude generated contents of those directories while retaining the placeholders. The `executables/` directory and MEAM parameter files are tracked.

Simulation outputs, experiment manifests, and generated `thermo.npz` archives are local data and are not tracked under the current `.gitignore` rules. Local MPI test logs should generally not be committed; optionally add `mpi_test.log` to `.gitignore`.
