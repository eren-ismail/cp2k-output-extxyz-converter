import os
import re

def parse_cp2k_output_to_xyz(input_file, output_file):
    """
    Converts CP2K output file to ASE Extended XYZ format with lattice vectors and atomic forces.
    
    Args:
        input_file (str): Path to the CP2K output file.
        output_file (str): Path to the ASE Extended XYZ output file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Extract total energy
    energy = None
    for line in lines:
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = (float(line.split()[-1]))*27.211386
            break

    # Locate the lattice vectors
    lattice_vectors = []
    for line in lines:
        if "CELL| Vector" in line:
            vector = list(map(float, line.split()[4:7]))
            lattice_vectors.append(vector)

    if len(lattice_vectors) != 3:
        raise ValueError("Could not find all three lattice vectors.")

    # Locate the atomic coordinates section
    coordinates_start = None
    for idx, line in enumerate(lines):
        if "MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom" in line or "MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM" in line:
            coordinates_start = idx + 2  # Skip the header
            break

    if coordinates_start is None:
        raise ValueError("Could not find atomic coordinates section in the file.")

    # Parse atomic coordinates
    atoms = []
    for line in lines[coordinates_start:]:
        if not line.strip():
            break
        match = re.match(r"\s*\d+\s+\d+\s+(\w+)\s+\d+\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", line)
        if match:
            atom_type = match.group(1)
            x, y, z = map(float, match.groups()[1:])
            atoms.append([atom_type, x, y, z])

    # Locate the atomic forces section
    forces_start = None
    for idx, line in enumerate(lines):
        if "ATOMIC FORCES in [a.u.]" in line:
            forces_start = idx + 2  # Forces start two lines after the header
            break

    if forces_start is None:
        raise ValueError("Could not find atomic forces section in the file.")

    # Parse atomic forces
    forces = []
    for line in lines[forces_start:]:
        if not line.strip():
            break
        match = re.match(r"\s*\d+\s+\d+\s+\w+\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", line)
        if match:
            fx, fy, fz = map(float, match.groups())
            forces.append([51.422 * fx, 51.422 * fy, 51.422 * fz])

    if len(atoms) != len(forces):
        raise ValueError("Mismatch between number of atoms and forces.")

    # Combine atoms and forces
    for i in range(len(atoms)):
        atoms[i].extend(forces[i])

    # Write the Extended XYZ file
    with open(output_file, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(
            f"Properties=species:S:1:pos:R:3:force:R:3 energy={energy} "
            f"Lattice=\"{' '.join(map(str, lattice_vectors[0]))} {' '.join(map(str, lattice_vectors[1]))} {' '.join(map(str, lattice_vectors[2]))}\""
            f' pbc="T T T"\n'
        )
        for atom in atoms:
            f.write(
                f"{atom[0]} {atom[1]:.12f} {atom[2]:.12f} {atom[3]:.12f} {atom[4]:.12f} {atom[5]:.12f} {atom[6]:.12f}\n"
            )

def process_all_subdirectories(base_directory):
    """
    Processes all subdirectories in the given base directory to convert CP2K output files
    and then combines all output.xyz files into a single XYZ file.
    
    Args:
        base_directory (str): Path to the base directory containing subdirectories.
    """
    # First loop: Create output.xyz files in each subdirectory
    for root, dirs, files in os.walk(base_directory):
        for file in files:
            if file == "output.out":
                input_path = os.path.join(root, file)
                output_path = os.path.join(root, "output.xyz")
                try:
                    print(f"Processing: {input_path}")
                    parse_cp2k_output_to_xyz(input_path, output_path)
                    print(f"Converted: {output_path}")
                except Exception as e:
                    print(f"Error processing {input_path}: {e}")

    # Second loop: Combine all output.xyz files into a single file
    combined_output_path = os.path.join(base_directory, "total_output_unit.xyz")
    with open(combined_output_path, 'w') as combined_file:
        for root, dirs, files in os.walk(base_directory):
            for file in files:
                if file == "output.xyz":
                    output_path = os.path.join(root, file)
                    try:
                        print(f"Combining: {output_path}")
                        with open(output_path, 'r') as single_xyz_file:
                            combined_file.write(single_xyz_file.read())
                        print(f"Added: {output_path}")
                    except Exception as e:
                        print(f"Error combining {output_path}: {e}")

# Example usage
process_all_subdirectories('.')  # Uncomment this line to process the current working directory and combine files
