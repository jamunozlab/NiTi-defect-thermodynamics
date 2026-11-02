# NiTi Defect Thermodynamics

This repository implements a lightweight, reproducible workflow for generating, running, extracting, and analyzing molecular-dynamics (MD) simulations of B2 NiTi using LAMMPS and a modified embedded-atom method (MEAM) potential. **UTEP (Utilities To Execute Pipelines)** organizes parameter sweeps into simulation directories and experiment manifests.

The workflow consists of five stages:

1. **Generate:** define parameter sweeps and generate structures and LAMMPS input files.
2. **Execute:** run LAMMPS simulations locally with MPI.
3. **Extract (ETL):** convert thermodynamic output in `log.lammps` to `thermo.npz`.
4. **Fit:** use `models/bm_eos.py` to fit a third-order Birch–Murnaghan equation of state to energy–volume data.
5. **Build a DataFrame:** run `build_dataframe.py` on one or more experiments and save their EOS summaries as a pandas pickle.

## Requirements

- Python 3, NumPy, pandas, SciPy, and Matplotlib
- LAMMPS compiled with MEAM support
- MPI implementation (for example, Open MPI)
- A Unix-like shell (tested on macOS)
- Jupyter Notebook or JupyterLab, **optional**, for interactive analysis

A working Conda environment can be activated with:

```bash
conda env list
conda activate moga
python -c "import numpy, pandas, scipy, matplotlib; print('Python dependencies available')"
```

The `moga` environment is an example of an existing environment that has been used successfully; it is **not** a required environment name. To install the Python dependencies into an appropriate Conda environment:

```bash
conda install numpy pandas scipy matplotlib
```

## Repository layout

```text
NiTi-defect-thermodynamics/
├── new_experiment.py                     # Define and run parameter sweeps
├── etl_experiment_lammps.py              # Convert LAMMPS logs to thermo.npz
├── build_dataframe.py                    # Build experiment-level EOS DataFrame
├── UTEP_LAMMPS_experiment_analysis.ipynb # Interactive analysis notebook
├── models/
│   └── bm_eos.py                         # LammpsExperiment and BM EOS fit
├── dataframes/                           # Generated pandas pickle files
├── utils/
│   ├── config.py                         # Machine-specific paths and executable
│   ├── default_dictionaries.py           # Default simulation parameters
│   ├── new_simulation_from_default.py
│   ├── write_input_from_default.py       # LAMMPS input and B2 structure
│   ├── write_run_sbatch_from_default.py
│   └── lammps_etl_utils.py               # Thermodynamic data extraction
├── executables/
│   ├── NiTi.library.meam
│   ├── NiTi.meam
│   └── lmp_mpi                            # Local executable symlink
├── simulations/                          # Generated simulation directories
└── experiments/                          # Experiment manifests
```

## Configure LAMMPS and UTEP

The `executables/lmp_mpi` entry is a **symbolic link**, not a bundled executable. In the original macOS setup it points to `/opt/homebrew/bin/lmp_mpi`. After cloning on another machine, install LAMMPS with MEAM support and recreate the link as necessary:

```bash
command -v lmp_mpi
rm executables/lmp_mpi
ln -s "$(command -v lmp_mpi)" executables/lmp_mpi
./executables/lmp_mpi -h
```

Check that the help output includes MEAM. Replacing a tracked symlink may appear as a local Git modification; avoid committing a machine-specific path unintentionally. Alternatively, configure the executable location in `utils/config.py`.

Edit `utils/config.py` to set local paths to the `simulations/`, `experiments/`, and `executables/` directories. The MEAM files must be available to each simulation:

```python
other_files = ["NiTi.library.meam", "NiTi.meam"]
```

Do not add `atoms_positions.data` to `other_files`: it is generated for each simulation.

Simulation defaults are defined in `utils/default_dictionaries.py`:

| Parameter | Meaning |
| --- | --- |
| `temperature` | Target temperature (K) |
| `timestep` | Integration timestep (ps in LAMMPS `metal` units) |
| `ntimesteps` | Number of integration steps |
| `logstep` | Thermodynamic output interval |
| `tdamp` | Thermostat damping time (ps) |
| `system_size` | Number of B2 cubic cells per box edge |
| `lattice_parameter` | Cubic B2 lattice parameter (Å) |
| `atom_type_to_reduce` | Species converted to the other species (1 = Ni, 2 = Ti) |
| `reduction_percentage` | Fraction of all sites converted |
| `defect_percentage` | Nominal antisite percentage |
| `atom_type_to_remove` | Species selected for vacancies (1 = Ni, 2 = Ti) |
| `vacancy_percentage` | Percentage of the selected species removed |
| `random_seed` | Random seed for reproducibility |

The B2 structure has two sites per cubic cell at `(0, 0, 0)` and `(1/2, 1/2, 1/2)`. Structure generation applies composition changes, antisite assignments, and vacancy removal in that order.

## 1. Generate and run an experiment

Edit the section between `#### SCRIPT SHOULD BE BELOW ####` and `#### SCRIPT SHOULD BE ABOVE ####` in `new_experiment.py`. For example:

```python
experiment_description = (
    "Pristine B2 NiTi lattice-parameter sweep\n"
    "Temperature: 400 K\n"
    "Lattice parameters: 2.81 to 3.20 Angstrom\n"
)

lattice_parameters = [round(2.81 + 0.01 * i, 2) for i in range(40)]

for lattice_parameter in lattice_parameters:
    sim_params_dict = {
        "temperature": 400,
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

From the repository root, generate simulation directories without executing LAMMPS:

```bash
python new_experiment.py lammps --test
```

To generate **new** directories and execute the simulations:

```bash
python new_experiment.py lammps
```

The current local launcher runs simulations sequentially; each simulation can use multiple MPI processes. The `experiments/experimentXXXXXX` manifest lists the associated simulation IDs. Each simulation directory contains `inputc`, `atoms_positions.data`, `runc.sh`, linked potential files, and (after execution) LAMMPS outputs.

## 2. Extract thermodynamic data (ETL)

After running an experiment:

```bash
python etl_experiment_lammps.py experiment000011
```

For additional messages:

```bash
python etl_experiment_lammps.py experiment000011 --verbose
```

The ETL code reads the thermodynamic table associated with the **last LAMMPS `run`**, not the preceding minimization, and writes:

```text
simulations/simulationXXXXXX/thermo.npz
```

Each archive contains:

| Array | Quantity | Units |
| --- | --- | --- |
| `Step` | Timestep index | dimensionless |
| `Temp` | Temperature | K |
| `TotEng` | Total energy of the simulation cell | eV |
| `Press` | Pressure | bar |

The ETL driver skips missing or invalid logs with a warning. Re-running ETL overwrites existing `thermo.npz` files for successfully processed simulations. Use `python etl_experiment_lammps.py --help` for additional options.

## 3. Fit a Birch–Murnaghan EOS

`models/bm_eos.py` provides the importable `LammpsExperiment` class. It loads the `thermo.npz` archives referenced by an experiment manifest, and fits a third-order Birch–Murnaghan equation to the average total energy per atom versus volume per atom.

For example, in Python or a Jupyter notebook:

```python
from models.bm_eos import LammpsExperiment

experiment = LammpsExperiment("experiment000011", repo_root=".")
fit = experiment.fit_birch_murnaghan(discard_fraction=0.5)
print(fit["a0_A"], fit["B0_GPa"], fit["B0_prime"])
experiment.plot_birch_murnaghan()
```

The default `discard_fraction=0.5` excludes the first half of each simulation's thermodynamic samples before averaging. The fit returns `E0_eV_per_atom`, `V0_A3_per_atom`, `a0_A`, `B0_eV_per_A3`, `B0_GPa`, and `B0_prime`, as well as fit parameters and covariance.

**Current limitations:** The EOS fitting method assumes cubic, pristine B2 supercells with exactly `2 * system_size**3` atoms. It is not yet generalized to vacancy-containing structures. The fit is to the finite-temperature **mean internal energy**, not the Helmholtz free energy; its fitted minimum and curvature should be interpreted accordingly.

## 4. Build and save an experiment-level DataFrame

`build_dataframe.py` runs from the **Terminal**; VS Code and Jupyter are not required. It accepts one or more experiment manifest names, runs the BM fit for each experiment, and saves **one row per experiment**.

For one experiment:

```bash
python build_dataframe.py --experiment experiment000011
```

Default output:

```text
dataframes/dataframe000011.pkl
```

For several experiments:

```bash
python build_dataframe.py --experiment experiment000011 experiment000010
```

Default output:

```text
dataframes/bm_experiments.pkl
```

Specify a custom output filename or a different initial fraction to discard:

```bash
python build_dataframe.py \
    --experiment experiment000011 experiment000010 \
    --output_filename bm_selected.pkl \
    --discard_fraction 0.5
```

The DataFrame includes:

| Column | Description |
| --- | --- |
| `experiment` | Experiment manifest name |
| `n_simulations`, `simulation_ids` | Number and IDs of simulations in the manifest |
| `temperature_K` | Target temperature parsed from `inputc`, if consistent across simulations |
| `x_Ni` | Ni atomic fraction calculated from the structure files, if consistent |
| `defect_type` | Currently `unknown` |
| `defect_concentration_percent` | Currently missing (`NaN`) |
| `E0_eV_per_atom`, `V0_A3_per_atom`, `a0_A` | Fitted EOS parameters |
| `B0_eV_per_A3`, `B0_GPa`, `B0_prime` | Fitted EOS modulus and pressure derivative |
| `discard_fraction` | Fraction of initial samples discarded |
| `fit_status`, `fit_error` | Fit success/failure and diagnostic message |

If a fit fails, the builder records `fit_status="failed"` and an error message instead of silently dropping the experiment. The builder does **not** rerun LAMMPS or perform ETL: `thermo.npz` files must already exist.

**Defect metadata is a work in progress.** Composition can be calculated from atom types in `atoms_positions.data`, but composition alone does not identify antisites or the original nominal defect percentage. The current script therefore leaves `defect_type` and `defect_concentration_percent` unspecified rather than guessing. A future revision will preserve these parameters explicitly when simulations are generated.

### Read the DataFrame in a notebook

```python
import pandas as pd

df = pd.read_pickle("dataframes/bm_experiments.pkl")
print(df[["experiment", "temperature_K", "x_Ni", "a0_A", "B0_GPa", "fit_status"]])
```

The `.pkl` format preserves pandas data types, including lists of simulation IDs. Only load pickle files from trusted sources.

## 5. Inspect thermodynamic time series

The `LammpsExperiment` class also provides time-series plotting:

```python
from models.bm_eos import LammpsExperiment

experiment = LammpsExperiment("experiment000011", repo_root=".")
experiment.plot("Temp")
experiment.plot("TotEng")
experiment.plot("Press")
```

`plot()` accepts `simulation_ids` (a string or list) and an existing Matplotlib `ax`. The repository also contains `UTEP_LAMMPS_experiment_analysis.ipynb` for interactive exploration.

## Troubleshooting

If Python reports `ModuleNotFoundError: No module named 'pandas'`, activate a Conda environment containing the dependencies or install them there:

```bash
conda env list
conda activate moga
python -c "import pandas, scipy, numpy, matplotlib"
```

To check syntax **without running the script**:

```bash
python -m py_compile build_dataframe.py
```

`py_compile` catches syntax errors, but does not execute imports or detect runtime errors. For available command-line arguments:

```bash
python build_dataframe.py --help
```

## Git and generated data

Simulation directories, experiment manifests, and thermodynamic archives are generated data. The repository uses `.gitignore` rules to exclude generated simulation and experiment contents while retaining `.gitkeep` placeholders. The `dataframes/` directory holds generated pickle outputs; check `git status` and `.gitignore` before committing any generated data. MEAM potential files are tracked, but the LAMMPS executable itself is not bundled.
