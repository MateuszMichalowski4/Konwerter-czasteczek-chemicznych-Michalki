import argparse
import glob
import sys
import time

def convert_charge(v3000_charge):
    charge_map = {0: 0, +3: 1, +2: 2, +1: 3, -1: 5, -2: 6, -3: 7}
    return charge_map.get(v3000_charge, 0)

def convert_cfg(v3000_cfg, bond_type):

    if bond_type == 1:
        cfg_map = {0: 0, 1: 1, 2: 4, 3: 6}
    elif bond_type == 2:
        cfg_map = {0: 0, 2: 3}
    else:
        cfg_map = {0: 0}
    return cfg_map.get(v3000_cfg, 0)

def parse_counts_line(line):
    parts = line.split()
    counts = {
        'num_atoms': int(parts[3]),
        'num_bonds': int(parts[4]),
        'num_sgroups': int(parts[5]),
        'num_3d_constraints': int(parts[6]),
        'chiral': int(parts[7]) == 1
    }
    return counts



def parse_v3000_atom(line):
    try:
        parts = line.split()
        atom = {
            'index': int(parts[2]),
            'element': parts[3],
            'x': float(parts[4]),
            'y': float(parts[5]),
            'z': float(parts[6]),
            'charge': 0,
            'stereo_parity': 0,
            'hydrogen_count': 0,
            'valence': 0,
            'isotopic_mass': 0,
            'radical': 0,
            'substitution_count': 0,
            'inversion_retention': 0,
            'exact_change_flag': 0,
            'atom_atom_mapping': 0,
            'stereo_box': 0,
            'unsaturation_flag': 0,
            'ring_bond_count': 0
        }

        for part in parts[7:]:
            if '=' in part:
                key, value = part.split('=')
                if key == 'CHG':
                    atom['charge'] = convert_charge(int(value))
                elif key == 'CFG':
                    atom['stereo_parity'] = int(value)
                elif key == 'HCOUNT':
                    atom['hydrogen_count'] = int(value)
                elif key == 'VAL':
                    atom['valence'] = int(value)
                elif key == 'MASS':
                    atom['isotopic_mass'] = int(value)
                elif key == 'RAD':
                    atom['radical'] = int(value)
                elif key == 'SUBST':
                    atom['substitution_count'] = int(value)
                elif key == 'INVRET':
                    atom['inversion_retention'] = int(value)
                elif key == 'EXACHG':
                    atom['exact_change_flag'] = int(value)
                elif key == 'AAMAP':
                    atom['atom_atom_mapping'] = int(value)
                elif key == 'STBOX':
                    atom['stereo_box'] = int(value)
                elif key == 'UNSAT':
                    atom['unsaturation_flag'] = int(value)
                elif key == 'RBCNT':
                    atom['ring_bond_count'] = int(value)

        return atom
    except IndexError as e:
        print(f"Error parsing atom line due to missing data: {line} - {str(e)}")
        return None
    except ValueError as e:
        print(f"Error parsing atom line due to incorrect data type: {line} - {str(e)}")
        return None





def parse_v3000_bond(line):
    try:
        parts = line.split()
        bond_type = int(parts[3])
        bond_cfg = 0
        bond = {
            'index': int(parts[2]),
            'type': bond_type,
            'start_atom': int(parts[4]),
            'end_atom': int(parts[5]),
            'cfg': bond_cfg,
            'topology': 0,
            'reacting_center': 0
        }
        for part in parts[6:]:
            if '=' in part:
                key, value = part.split('=')
                if key == 'CFG':
                    bond['cfg'] = convert_cfg(int(value), bond_type)
                elif key == 'TOPO':
                    bond['topology'] = int(value)
                elif key == 'RXCTR':
                    bond['reacting_center'] = int(value)
        return bond
    except IndexError as e:
        print(f"Error during bond line analysis due to missing data: {line} - {str(e)}")
        return None
    except ValueError as e:
        print(f"Error during bond line analysis due to invalid data type: {line} - {str(e)}")
        return None




def parse_properties(lines, start_index):
    properties = {}
    i = start_index
    while i < len(lines) and 'M  END' not in lines[i]:
        line = lines[i].strip()
        if line.startswith('M  V30'):
            parts = line.split()
            key = parts[1]
            if key not in properties:
                properties[key] = []
            if len(parts) > 3:
                count = int(parts[2])
                for j in range(count):
                    index = 3 + 2 * j
                    properties[key].append({
                        'atom_no': int(parts[index]),
                        'value': int(parts[index + 1])
                    })
        i += 1
    return properties


def parse_mol_file(file_path):
    molecules = []
    molecule = None
    in_atom_block = False
    in_bond_block = False
    in_property_block = False

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if 'M  V30 BEGIN CTAB' in line:
            molecule = {
                'atoms': [],
                'bonds': [],
                'properties': {},
                'counts': {}
            }
        elif 'M  V30 END CTAB' in line and molecule is not None:
            molecules.append(molecule)
            molecule = None
        elif molecule is not None:
            if 'M  V30 COUNTS' in line:
                molecule['counts'] = parse_counts_line(line)
            elif 'M  V30 BEGIN ATOM' in line:
                in_atom_block = True
            elif 'M  V30 END ATOM' in line:
                in_atom_block = False
            elif 'M  V30 BEGIN BOND' in line:
                in_bond_block = True
            elif 'M  V30 END BOND' in line:
                in_bond_block = False
            elif 'M  V30 BEGIN' in line:
                in_property_block = True
            elif 'M  V30 END' in line:
                in_property_block = False
            elif in_atom_block:
                molecule['atoms'].append(parse_v3000_atom(line))
            elif in_bond_block:
                molecule['bonds'].append(parse_v3000_bond(line))
            elif in_property_block:
                molecule['properties'].update(parse_properties(lines, lines.index(line)))

    return molecules


def mol_to_sdf(molecule):
    num_atoms = len(molecule['atoms'])
    num_bonds = len(molecule['bonds'])
    num_atom_lists = 0
    chiral_flag = 1 if molecule['counts']['chiral'] else 0
    num_stext_entries = 0

    obsolete = 0
    num_properties = 999

    header = f"{num_atoms:3d}{num_bonds:3d}{num_atom_lists:3d}{obsolete:3d}{chiral_flag:3d}{num_stext_entries:3d}"
    header += f"{obsolete:3d}{obsolete:3d}{obsolete:3d}{obsolete:3d}{num_properties:3d} V2000"

    lines = []
    lines.append(molecule.get('name', 'Untitled'))
    lines.append('KCCW V3000 to SDF conversion')
    lines.append('')
    lines.append('')
    lines.append(header)

    for atom in molecule['atoms']:
        atom_line = f"{atom['x']:10.4f}{atom['y']:10.4f}{atom['z']:10.4f} {atom['element']:>2} 0 {atom['charge']:2d}  {atom['stereo_parity']:2d}  {atom['hydrogen_count']:2d}  {atom['stereo_box']:2d} {atom['valence']:2d} 0  0  0  0  0 {atom['exact_change_flag']:2d}"
        lines.append(atom_line)

    for bond in molecule['bonds']:
        bond_line = f"{bond['start_atom']:3d}{bond['end_atom']:3d}{bond['type']:3d} {bond['cfg']:2d} 0 {bond['topology']:2d} {bond['reacting_center']:2d}"
        lines.append(bond_line)

    lines.append('M  END')
    return '\n'.join(lines)
def convert_charge1(sdf_charge):
    charge_map = {0: 0, 1: +3, 2: +2, 3: +1, 5: -1, 6: -2, 7: -3}
    return charge_map.get(sdf_charge, 0)

def convert_cfg1(sdf_cfg, bond_type):
    if bond_type == 1:
        cfg_map = {0: 0, 1: 1, 2: 4, 3: 6}
    elif bond_type == 2:
        cfg_map = {0: 0, 2: 3}
    else:
        cfg_map = {0: 0}
    return cfg_map.get(sdf_cfg, 0)

def parse_sdf_atom(line):
    try:
        parts = line.split()
        atom = {
            'x': float(parts[0]),
            'y': float(parts[1]),
            'z': float(parts[2]),
            'element': parts[3],
            'mass_difference': int(parts[4]) if len(parts) > 4 else 0,
            'charge': int(parts[5]) if len(parts) > 5 else 0,
            'stereo_parity': int(parts[6]) if len(parts) > 6 else 0,
            'hydrogen_count': int(parts[7]) if len(parts) > 7 else 0,
            'stereo_care_box': int(parts[8]) if len(parts) > 8 else 0,
            'valence': int(parts[9]) if len(parts) > 9 else 0,
            'h0_designator': int(parts[10]) if len(parts) > 10 else 0,
            'atom_atom_mapping_number': int(parts[11]) if len(parts) > 11 else 0,
            'inversion_retention_flag': int(parts[12]) if len(parts) > 12 else 0,
            'exact_change_flag': int(parts[13]) if len(parts) > 13 else 0
        }
        return atom
    except IndexError as e:
        print(f"Error parsing atom line due to missing data: {line} - {str(e)}")
        return None
    except ValueError as e:
        print(f"Error parsing atom line due to incorrect data type: {line} - {str(e)}")
        return None

def parse_sdf_bond(line):
    try:
        parts = line.split()
        bond = {
            'start_atom': int(parts[0]),
            'end_atom': int(parts[1]),
            'type': int(parts[2]),
            'stereo': int(parts[3]) if len(parts) > 3 else 0,
            'topology': int(parts[4]) if len(parts) > 4 else 0,
            'reacting_center': int(parts[5]) if len(parts) > 5 else 0
        }
        return bond
    except IndexError as e:
        print(f"Error parsing bond line due to missing data: {line} - {str(e)}")
        return None
    except ValueError as e:
        print(f"Error parsing bond line due to incorrect data type: {line} - {str(e)}")
        return None


def parse_sdf_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()

    molecules_data = content.split('$$$$\n')
    molecules = []

    for molecule_data in molecules_data:
        if molecule_data.strip() == "":
            continue

        lines = molecule_data.split('\n')
        molecule = {
            'atoms': [],
            'bonds': [],
            'name': '',
            'counts': {}
        }

        atom_block = False
        bond_block = False
        name_next_line = False
        atom_lines = []
        bond_lines = []

        for line in lines:
            if 'V2000' in line:
                molecule['counts'] = {
                    'num_atoms': int(line[0:3].strip()),
                    'num_bonds': int(line[3:6].strip()),
                    'num_atom_lists': int(line[6:9].strip()),
                    'chiral': int(line[12:15].strip()) == 1,
                    'num_sgroups': int(line[15:18].strip())
                }
                atom_block = True
                bond_block = False
                atom_count = molecule['counts']['num_atoms']
                bond_count = molecule['counts']['num_bonds']
            elif atom_block and atom_count > 0:
                atom_lines.append(line)
                atom_count -= 1
                if atom_count == 0:
                    atom_block = False
                    bond_block = True
            elif bond_block and bond_count > 0:
                bond_lines.append(line)
                bond_count -= 1
                if bond_count == 0:
                    bond_block = False
            elif line.strip().startswith('> <NAME>'):
                name_next_line = True
            elif name_next_line:
                molecule['name'] = line.strip()
                name_next_line = False

        for line in atom_lines:
            atom = parse_sdf_atom(line)
            if atom:
                molecule['atoms'].append(atom)

        for line in bond_lines:
            bond = parse_sdf_bond(line)
            if bond:
                molecule['bonds'].append(bond)

        molecules.append(molecule)

    return molecules


def sdf_to_v3000(molecule):
    if 'counts' not in molecule or 'num_atoms' not in molecule['counts'] or 'num_bonds' not in molecule['counts']:
        print(f"Missing counts in molecule data: {molecule}")
        return ""

    v3000_lines = []
    v3000_lines.append(molecule.get('name', 'Untitled'))
    v3000_lines.append('KCCW SDF to V3000 conversion')
    v3000_lines.append("\n 0 0 0        0  0       999 V3000")
    v3000_lines.append("M  V30 BEGIN CTAB")

    num_atoms = molecule['counts']['num_atoms']
    num_bonds = molecule['counts']['num_bonds']
    num_sgroups = molecule['counts'].get('num_sgroups', 0)
    num_3d_constraints = 0
    chiral = 1 if molecule['counts'].get('chiral', False) else 0

    v3000_lines.append(f"M  V30 COUNTS {num_atoms} {num_bonds} {num_sgroups} {num_3d_constraints} {chiral}")

    v3000_lines.append("M  V30 BEGIN ATOM")
    for i, atom in enumerate(molecule['atoms'], start=1):
        parts = [f"M  V30 {i} {atom['element']} {atom['x']} {atom['y']} {atom['z']}"]
        if atom.get('charge', 0) != 0:
            parts.append(f"CHG={convert_charge1(atom['charge'])}")
        if atom.get('stereo_parity', 0) != 0:
            parts.append(f"CFG={atom['stereo_parity']}")
        if atom.get('hydrogen_count', 0) != 0:
            parts.append(f"HCOUNT={atom['hydrogen_count']}")
        if atom.get('valence', 0) != 0:
            parts.append(f"VAL={atom['valence']}")
        if atom.get('exact_change_flag', 0) != 0:
            parts.append(f"EXACHG={atom['exact_change_flag']}")
        v3000_lines.append(' '.join(parts))
    v3000_lines.append("M  V30 END ATOM")


    v3000_lines.append("M  V30 BEGIN BOND")
    for i, bond in enumerate(molecule['bonds'], start=1):
        parts = [f"M  V30 {i} {bond['type']} {bond['start_atom']} {bond['end_atom']}"]
        if bond.get('stereo', 0) != 0:
            parts.append(f"CFG={convert_cfg1(bond['stereo'], bond['type'])}")
        if bond.get('topology', 0) != 0:
            parts.append(f"TOPO={bond['topology']}")
        if bond.get('reacting_center', 0) != 0:
            parts.append(f"RXCTR={bond['reacting_center']}")
        v3000_lines.append(' '.join(parts))
    v3000_lines.append("M  V30 END BOND")


    v3000_lines.append("M  V30 END CTAB")
    v3000_lines.append("M  END")
    return '\n'.join(v3000_lines)


def convert_files(input_files, output_path, conversion_type):
    start_time = time.time()
    with open(output_path, 'w') as output_file:
        for file_path in input_files:
            if conversion_type == 'sdf_to_v3000':
                molecules = parse_sdf_file(file_path)
                for molecule in molecules:
                    v3000_data = sdf_to_v3000(molecule)
                    output_file.write(v3000_data + '\n$$$$\n')
                    sys.stderr.write(f"Written MOLv3000 data for {file_path} to {output_path}\n")
            elif conversion_type == 'v3000_to_sdf':
                molecules = parse_mol_file(file_path)
                for molecule in molecules:
                    sdf_data = mol_to_sdf(molecule)
                    output_file.write(sdf_data + '\n$$$$\n')
                    sys.stderr.write(f"Written SDF data for {file_path} to {output_path}\n")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken to convert files: {elapsed_time:.4f} seconds")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert between V3000 and SDF formats.')
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # Parser for conversion
    convert_parser = subparsers.add_parser('convert', help='Convert files between formats')
    convert_parser.add_argument('-i', '--input', required=True, help='Input file pattern, e.g., *.mol or *.sdf')
    convert_parser.add_argument('-o', '--output', required=True, help='Output file path for the converted data')
    convert_parser.add_argument('-t', '--type', required=True, choices=['v3000_to_sdf', 'sdf_to_v3000'],
                                help='Conversion type: v3000_to_sdf or sdf_to_v3000')

    # Add help command
    help_parser = subparsers.add_parser('help', help='Show this help message and exit')

    args = parser.parse_args()

    if args.command == 'help':
        parser.print_help()
    elif args.command == 'convert':
        input_files = glob.glob(args.input)
        if not input_files:
            sys.stderr.write("No files matched the pattern.\n")
        else:
            convert_files(input_files, args.output, args.type)


