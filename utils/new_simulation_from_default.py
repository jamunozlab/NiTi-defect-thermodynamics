#!/usr/bin/env python

"""
Sets up a new simulation given a simulations path, executable path and filename,
and input and sbatch parameters to change from a default.
"""

import os
from .dir_structure_utils import find_next_available_dir
from .write_input_from_default import write_input_from_default
from .write_run_sbatch_from_default import write_run_sbatch_from_default_mac

def new_simulation_from_default(default, simulations_path, executables_path, executable, sim_params_dict, 
                                sbatch_options_dict, other_files=[], submit_job=True):
    os.chdir(simulations_path)
    simulation_id = find_next_available_dir('simulation', simulations_path)
    simulation_path = simulations_path + simulation_id + '/'
    os.mkdir(simulation_path)
    executable_path = executables_path + executable
    new_executable_path = simulation_path + executable
    command = 'ln -s ' + executable_path + ' ' + new_executable_path
    os.system(command)
    for other_file in other_files:
        other_path = executables_path + other_file
        new_other_path = simulation_path + other_file
        command = 'ln -s ' + other_path + ' ' + new_other_path
        os.system(command)
    write_input_from_default(default, simulation_path, **sim_params_dict)
    write_run_sbatch_from_default_mac(
        simulation_path=simulation_path,
        executable=executable,
        **sbatch_options_dict
    )
    print('New simulation is in ' + simulation_path)
    
    #if submit_job:
    #    os.chdir(simulation_path)
    #    os.system('sbatch ./runc.sbatch')

    if submit_job:
        import subprocess
        subprocess.run(
        ["bash", "./runc.sh"],
        cwd=simulation_path,
        check=True
    )

    return simulation_id