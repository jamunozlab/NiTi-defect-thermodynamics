#!/usr/bin/env python

"""
This script is responsible for iterating over the variables that will be modified in an experiment,
creating a directory for each combination of interest with all the files necessary to run an individual
simulation in that directory, and submitting the job to the Slurm worload manager.

The user must provide code to modify the variables of interest in sim_params_dict and sbatch_options_dict,
as well as a description of the experiment. The function new_simulation_from_default should be called
once for each simulation. The script should be executed from the command line. 
"""

import argparse

from utils.new_simulation_from_default import new_simulation_from_default
from utils.dir_structure_utils import find_next_available_file
from utils.config import simulations_path, executables_path, executable, experiments_path, other_files
import utils.default_dictionaries as default_dictionaries

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='new_experiment',
                    description='Create new experiment. User must provide the name of a dictionary' +
                                    'in default_dictionaries.py with input parameters and default values.',
                    epilog='')
    
    parser.add_argument('default_input_dictionary_name') # positional argument
    parser.add_argument('-t', '--test', action='store_false')  # on/off flag
    
    submit_job=False

    args = parser.parse_args()
    submit_job = args.test
    default_input_dictionary = getattr(
        default_dictionaries,
        args.default_input_dictionary_name
    )
    
    simulation_ids = []


    #### SCRIPT SHOULD BE BELOW ####


    # Ni vacancy-count sweep in initially equiatomic B2 NiTi.
    # One simulation for each vacancy count: 0, 1, ..., 8.

    system_size = 20
    temperature = 300
    lattice_parameter = 3.00

    n_Ni_sites = system_size**3

    experiment_description = (
        "NiTi LAMMPS Ni-sublattice vacancy-count sweep\n"
        "Initial composition: equiatomic B2 NiTi (50% Ni, 50% Ti)\n"
        f"Temperature: {temperature} K\n"
        f"Lattice parameter: {lattice_parameter:.2f} Angstrom\n"
        f"System size: {system_size} x {system_size} x {system_size}\n"
        "Vacancy species: Ni (atom type 1)\n"
        "Vacancy counts: 0 through 8 inclusive\n"
        "No composition adjustment or antisite disorder\n"
        "MPI processes per simulation: 4\n"
    )

    for n_vacancies in range(9):

        if n_vacancies == 0:
            vacancy_percentage = 0.0
        else:
            vacancy_percentage = (
                100.0 * (n_vacancies + 0.5) / n_Ni_sites
            )

        sim_params_dict = {
            "temperature": temperature,
            "lattice_parameter": lattice_parameter,
            "system_size": system_size,
            "atom_type_to_reduce": 1,
            "reduction_percentage": 0.0,
            "defect_percentage": 0.0,
            "atom_type_to_remove": 1,
            "vacancy_percentage": vacancy_percentage,
            "random_seed": 12345,
        }

        sbatch_options_dict = {
            "nproc": 4,
        }

        print(
            f"Creating simulation with {n_vacancies} Ni vacancies "
            f"({vacancy_percentage:.6f}% of Ni sites)"
        )

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

    experiment_id = find_next_available_file(
        "experiment", experiments_path
    )


    #### SCRIPT SHOULD BE ABOVE ####

    lines = experiment_description + str('\n')
    for simulation_id in simulation_ids:
        lines = lines + simulation_id + str('\n')

    with open(experiments_path + experiment_id, "w") as f:
        f.writelines(lines)
    
    print(experiment_id)
    print(simulation_ids)
    
