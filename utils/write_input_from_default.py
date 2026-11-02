#!/usr/bin/env python

"""
Auxiliary functions to write input files
"""

import numpy as np
from pathlib import Path

def generate_input_from_default(default, **kwargs):
    simulation_parameters_dictionary = default.copy()
    for key, value in kwargs.items():
        simulation_parameters_dictionary[key] = value 
        
    return simulation_parameters_dictionary

def generate_input_lines(param_dict):

    temperature = param_dict["temperature"]
    ntimesteps = param_dict["ntimesteps"]
    tdamp = param_dict["tdamp"]
    logstep = param_dict["logstep"]

    input_text = f"""units         metal
dimension     3
boundary      p p p
atom_style    atomic
read_data     "atoms_positions.data"

# -------------------- INTERATOMIC POTENTIAL -------------------------

pair_style    meam
pair_coeff    * * NiTi.library.meam Ni Ti NiTi.meam Ni Ti

mass 1 58.690
mass 2 47.880

############################ EQUILIBRATION ############################

reset_timestep 0
timestep 0.005

minimize 1e-6 1e-8 10000 100000
write_data minimized_structure.data

dump 4a all custom 500 dump.NiTi.* id type x y z fx fy fz
dump_modify 4a pbc yes sort id

############################ MD ######################################

reset_timestep 0

# Set temperature
velocity all create {temperature} 12345 mom yes rot no

reset_timestep 0

fix 3 all nvt temp {temperature} {temperature} {tdamp}

# Setting thermo
thermo {logstep}

# Run
run {ntimesteps}

# SIMULATION DONE
print "All done"
"""

    return input_text.splitlines(keepends=True)

def write_bcc_poscar(simulation_path, a_val, system_name="Fe"):
    lines = [
        f"For VASP. System= {system_name}\n",
        f"{a_val:.16f}\n",
        "1.0000000000000000 0.0000000000000000 0.0000000000000000\n",
        "0.0000000000000000 1.0000000000000000 0.0000000000000000\n",
        "0.0000000000000000 0.0000000000000000 1.0000000000000000\n",
        "2\n",
        "Direct\n",
        "0.0000000000000000 0.0000000000000000 0.0000000000000000\n",
        "0.5000000000000000 0.5000000000000000 0.5000000000000000\n",
    ]

    with open(simulation_path + "POSCAR", "w") as f:
        f.writelines(lines)

# Defect functions adapted from UTEP-tedehamo/utils/gen_defect_sim.py.
# Their original signatures, return values, and selection rules are preserved.
def generate_atom_type_array(total_atoms):
    """Generates an atom type array for an initial B2 ordered structure."""
    atom_type_array = np.zeros(total_atoms, dtype=int)
    for i in range(total_atoms):
        atom_type_array[i] = 2 if (i + 1) % 2 == 0 else 1
    num_type_1 = np.sum(atom_type_array == 1)
    num_type_2 = np.sum(atom_type_array == 2)
    return atom_type_array, num_type_1, num_type_2


def reduce_concentration(atom_type_array, atom_type_to_reduce, reduction_percentage):
    """Change abs(reduction_percentage) of all sites from the specified type to the other.

    As in the original active implementation, the sign does not reverse the direction.
    """
    atom_type_to_reduce = int(atom_type_to_reduce)
    atom_type_array = np.array(atom_type_array)
    if atom_type_to_reduce not in [1, 2]:
        raise ValueError("Only supports atom types 1 and 2.")
    other_type = 2 if atom_type_to_reduce == 1 else 1
    from_type = atom_type_to_reduce
    to_type = other_type
    indices_to_modify = np.where(atom_type_array == from_type)[0]
    num_to_change = int(len(atom_type_array) * abs(reduction_percentage))
    if num_to_change > len(indices_to_modify):
        raise ValueError("Too few atoms to change.")
    if num_to_change == 0:
        print("No atoms changed due to small reduction_percentage.")
        return atom_type_array, np.sum(atom_type_array == 1), np.sum(atom_type_array == 2)
    selected_indices = np.random.choice(indices_to_modify, num_to_change, replace=False)
    atom_type_array[selected_indices] = to_type
    return atom_type_array, np.sum(atom_type_array == 1), np.sum(atom_type_array == 2)


def generate_defect_positions(total_atoms, atom_type_array, defect_percentage):
    """Flip a percentage of each currently present species (original algorithm)."""
    total_defects = int(total_atoms * defect_percentage / 100)  # original, unused
    indices_type_1 = np.where(atom_type_array == 1)[0]
    indices_type_2 = np.where(atom_type_array == 2)[0]
    num_defects_1 = int(len(indices_type_1) * defect_percentage / 100)
    num_defects_2 = int(len(indices_type_2) * defect_percentage / 100)
    defect_indices_1 = np.random.choice(indices_type_1, num_defects_1, replace=False)
    defect_indices_2 = np.random.choice(indices_type_2, num_defects_2, replace=False)
    modified_atom_type_array = atom_type_array.copy()
    modified_atom_type_array[defect_indices_1] = 2
    modified_atom_type_array[defect_indices_2] = 1
    num_type_1 = np.sum(modified_atom_type_array == 1)
    num_type_2 = np.sum(modified_atom_type_array == 2)
    return modified_atom_type_array, defect_indices_1, defect_indices_2, num_type_1, num_type_2


def introduce_vacancies(atom_type_array, atom_type_to_remove, vacancy_percentage):
    """Select vacancy sites from one species; percentage is on a 0-100 scale.

    Return the original-array indices to remove, so coordinates stay aligned.
    """
    atom_type_to_remove = int(atom_type_to_remove)
    if atom_type_to_remove not in (1, 2):
        raise ValueError("atom_type_to_remove must be 1 (Ni) or 2 (Ti)")
    if not 0 <= vacancy_percentage <= 100:
        raise ValueError("vacancy_percentage must be between 0 and 100")
    candidates = np.where(np.asarray(atom_type_array) == atom_type_to_remove)[0]
    count = int(len(candidates) * vacancy_percentage / 100)
    return np.random.choice(candidates, count, replace=False)


def write_b2_lammps_data(simulation_path, system_size, lattice_parameter, atom_types=None, vacancy_indices=None):
    """Write a B2 NiTi supercell: type 1 = Ni, type 2 = Ti."""
    system_size = int(system_size)
    lattice_parameter = float(lattice_parameter)
    if system_size < 1 or lattice_parameter <= 0:
        raise ValueError("system_size and lattice_parameter must be positive")

    total_atoms = 2 * system_size**3
    if atom_types is None:
        atom_types = np.tile([1, 2], system_size**3)
    atom_types = np.asarray(atom_types)
    if len(atom_types) != total_atoms or not np.all(np.isin(atom_types, [1, 2])):
        raise ValueError("atom_types must contain exactly 2*system_size**3 entries of type 1 or 2")

    basis = np.eye(3) * lattice_parameter
    base_atoms = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * lattice_parameter
    positions = []
    for i in range(system_size):
        for j in range(system_size):
            for k in range(system_size):
                base_position = np.array([i, j, k])
                cart_position = np.dot(basis.T, base_position)
                for atom in base_atoms:
                    positions.append(cart_position + atom)

    if vacancy_indices is not None:
        vacancy_indices = np.asarray(vacancy_indices, dtype=int)
        if np.any(vacancy_indices < 0) or np.any(vacancy_indices >= total_atoms):
            raise ValueError("vacancy_indices contain an out-of-range index")
        if len(np.unique(vacancy_indices)) != len(vacancy_indices):
            raise ValueError("vacancy_indices must be unique")
        keep = np.ones(total_atoms, dtype=bool)
        keep[vacancy_indices] = False
        atom_types = atom_types[keep]
        positions = np.asarray(positions)[keep]

    num_atoms = len(atom_types)
    box_length = system_size * lattice_parameter
    output_path = Path(simulation_path) / "atoms_positions.data"
    with open(output_path, "w") as f:
        f.write("B2 NiTi structure\n\n")
        f.write(f"{num_atoms} atoms\n2 atom types\n\n")
        for axis in "xyz":
            f.write(f"0.0 {box_length:.8f} {axis}lo {axis}hi\n")
        f.write("\nMasses\n\n1 58.6934\n2 47.867\n")
        f.write("\nAtoms # atomic\n\n")
        for atom_id, (atom_type, position) in enumerate(zip(atom_types, positions), start=1):
            x, y, z = position
            f.write(f"{atom_id} {int(atom_type)} {x:.8f} {y:.8f} {z:.8f}\n")
    return output_path

def write_input_from_default(default, simulation_path, **kwargs):
    simulation_parameters_dictionary = generate_input_from_default(default, **kwargs)
    input_lines = generate_input_lines(param_dict=simulation_parameters_dictionary)

    with open(simulation_path + "inputc", "w") as f:
        f.writelines(input_lines)

    system_size = int(simulation_parameters_dictionary["system_size"])
    total_atoms = 2 * system_size**3
    atom_types = simulation_parameters_dictionary.get("atom_types")

    vacancy_indices = None
    if atom_types is None:
        # Restore NumPy's random state afterwards to avoid affecting other code.
        previous_random_state = np.random.get_state()
        try:
            np.random.seed(int(simulation_parameters_dictionary.get("random_seed", 12345)))
            atom_types, _, _ = generate_atom_type_array(total_atoms)
            atom_types, _, _ = reduce_concentration(
                atom_types,
                simulation_parameters_dictionary.get("atom_type_to_reduce", 1),
                simulation_parameters_dictionary.get("reduction_percentage", 0.0),
            )
            atom_types, _, _, num_type_1, num_type_2 = generate_defect_positions(
                total_atoms,
                atom_types,
                simulation_parameters_dictionary.get("defect_percentage", 0.0),
            )
            vacancy_indices = introduce_vacancies(
                atom_types,
                simulation_parameters_dictionary.get("atom_type_to_remove", 1),
                simulation_parameters_dictionary.get("vacancy_percentage", 0.0),
            )
        finally:
            np.random.set_state(previous_random_state)
        print(f"Ni atoms: {num_type_1}; Ti atoms: {num_type_2}")

    write_b2_lammps_data(
        simulation_path,
        system_size,
        simulation_parameters_dictionary["lattice_parameter"],
        atom_types=atom_types,
        vacancy_indices=vacancy_indices,
    )