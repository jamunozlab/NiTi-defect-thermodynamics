#!/usr/bin/env python

"""
Auxiliary function to write sbatch files
"""

def generate_run_sbatch_lines(simulation_path, executable, options_dict):
    lines = []
    lines.append("#!/bin/bash" + str('\n'))
    for key, value in options_dict.items():
        lines.append("#SBATCH " + str(key) + " " + value + str('\n'))
    lines.append("export SLURM_CPU_BIND=\"cores\"" + str('\n'))
    first  = "srun " + simulation_path + executable + " "
    lines.append(first + 'inputc' + " > " + 'output.txt')
    return lines

def generate_run_sbatch_from_default(**kwargs):
    slurm_options_dictionary = { '-A' : 'm3766_g',
                          '-C' : 'gpu',
                          '-q' : 'shared',
                          '-t' : '06:00:00',
                          '-c' : '32',
             '--gpus-per-task=1' : '',
  #         '--ntasks-per-node' : '4',
                          '-n' : '1'
                         }
    for key, value in kwargs.items():
        slurm_options_dictionary[key] = value 
    return slurm_options_dictionary
    
def write_run_sbatch_from_default(simulation_path, executable, **kwargs):
    slurm_options_dict = generate_run_sbatch_from_default(**kwargs)
    run_sbatch_lines = generate_run_sbatch_lines(simulation_path, executable, options_dict=slurm_options_dict)
    f = open(simulation_path+"runc.sbatch", "w")
    f.writelines(run_sbatch_lines)
    f.close()

def generate_run_sbatch_from_default_jakar(**kwargs):
    """
    Default Slurm options for Jakar.
    Keyword arguments can override any default.
    """

    slurm_options_dictionary = {
        "--job-name": "moga",
        "--account": "jakar_partner_physj",
        "--qos": "jakar_medium_physj",
        "--partition": "medium",
        "--nodes": "1",
        "--ntasks-per-node": "1",
        "--cpus-per-task": "2",
        "--time": "06:00:00",
        "--output": "%j.log",
        "--error": "%j.err",
    }

    for key, value in kwargs.items():
        slurm_options_dictionary[key] = value

    return slurm_options_dictionary


def generate_run_sbatch_lines_jakar(
    simulation_path="/Users/jamunoz/Documents/GitHub/MOGA-Phonons",
    executable="pygad_module_2_6_50_mod_parallel.py",
    conda_path="/scratch/$USER/miniconda3/etc/profile.d/conda.sh",
    conda_env="moga",
    modules=("gnu12", "openmpi4/4.1.5"),
    options_dict=None,
):
    """
    Generate a Jakar-compatible runc.sbatch script.
    """

    if options_dict is None:
        options_dict = generate_run_sbatch_from_default_jakar()

    lines = []

    lines.append("#!/bin/bash\n")

    for key, value in options_dict.items():
        if value == "":
            lines.append(f"#SBATCH {key}\n")
        else:
            lines.append(f"#SBATCH {key}={value}\n")

    lines.append("\n")
    lines.append('echo "========================================"\n')
    lines.append('echo "Job started on $(hostname)"\n')
    lines.append('echo "Date: $(date)"\n')
    lines.append('echo "Working dir: $(pwd)"\n')
    lines.append('echo "SLURM_NTASKS: $SLURM_NTASKS"\n')
    lines.append('echo "========================================"\n')
    lines.append("\n")

    lines.append("ulimit -s unlimited\n")
    lines.append("\n")

    lines.append("# Load MPI\n")
    for module in modules:
        lines.append(f"module load {module}\n")
    lines.append("\n")

    lines.append("# Activate conda\n")
    lines.append(f"source {conda_path}\n")
    lines.append(f"conda activate {conda_env}\n")
    lines.append("\n")

    lines.append("# Diagnostics\n")
    lines.append("which python\n")
    lines.append("python --version\n")
    lines.append("which mpirun\n")
    lines.append('python -c "import numpy, mpi4py; print(\\"Env OK\\")"\n')
    lines.append("\n")

    lines.append("# Go to your code\n")
    lines.append(f"cd {simulation_path}\n")
    lines.append("\n")

    lines.append('echo "Running MOGA..."\n')
    lines.append("\n")

    lines.append(
        f"python {executable} > output.txt\n"
    )
    lines.append("\n")

    lines.append('echo "========================================"\n')
    lines.append('echo "Job finished at $(date)"\n')
    lines.append('echo "========================================"\n')

    return lines


def write_run_sbatch_from_default_jakar(
    simulation_path="/Users/jamunoz/Documents/GitHub/MOGA-Phonons",
    executable="pygad_module_2_6_50_mod_parallel.py",
    output_filename="runc.sbatch",
    **kwargs,
):
    """
    Write a Jakar-compatible runc.sbatch file.
    """

    slurm_options_dict = generate_run_sbatch_from_default_jakar(**kwargs)

    run_sbatch_lines = generate_run_sbatch_lines_jakar(
        simulation_path=simulation_path,
        executable=executable,
        options_dict=slurm_options_dict,
    )

    with open(f"{simulation_path.rstrip('/')}/{output_filename}", "w") as f:
        f.writelines(run_sbatch_lines)

def write_run_sbatch_from_default_mac(
    simulation_path,
    executable,
    **kwargs
):
    """
    Write a local MPI execution script for LAMMPS on macOS.
    """

    from pathlib import Path

    # Default execution options
    options = {
        "nproc": 4
    }
    options.update(kwargs)

    nproc = int(options["nproc"])

    simulation_path = Path(simulation_path).resolve()

    lines = [
        "#!/bin/bash\n",
        "set -e\n",
        f'cd "{simulation_path}"\n',
        f"mpirun -np {nproc} ./{executable} -in inputc > output.txt 2>&1\n"
    ]

    output_file = simulation_path / "runc.sh"

    with open(output_file, "w") as f:
        f.writelines(lines)

    output_file.chmod(0o755)