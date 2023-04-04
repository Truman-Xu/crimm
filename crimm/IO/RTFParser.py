import warnings

def skip_line(x):
    """Return if the line should be skipped (empty lines or comment only)"""
    return (x.startswith('!') or x.startswith('*') or x.strip() == '')

def comment_parser(line):
    """Parse and separate data fields and comments"""
    if line.find('!') == -1:
        # line has no comments
        return line, None
    # return fields string, comments
    line, desc = line.split('!', maxsplit = 1)[:2]
    return line, desc.strip()

def mass_parser(line):
    """Parse lines with keyword MASS"""
    field_str, desc = comment_parser(line)
    # fields: [key, deprecated_entry, atom_type, atomic_mass]
    key, _, atom_type, mass = field_str.split()
    return atom_type, mass, desc

def decl_parser(line):
    """Parse line with keyword DECLare"""
    field_str = comment_parser(line)[0]
    atom = field_str.split()[-1]
    return atom

def defa_parser(line):
    """Parse default terminal patches"""
    field_str = comment_parser(line)[0]
    fields = field_str.split()
    # find where these two keyword located in the list
    first_i = fields.index('FIRS')
    last_i = fields.index('LAST')
    first_patch = fields[first_i+1]
    last_patch = fields[last_i+1]
    return first_patch, last_patch

def auto_parser(line):
    """Parse keyword AUTOGEN"""
    field_str = comment_parser(line)[0]
    autogen_ops = field_str.split()[1:]
    return autogen_ops

def resi_parser(line):
    """Parse the keyword RESI for a new residue"""
    field_str, desc = comment_parser(line)
    # [key, residue_name, total_charge]
    key, res_name, total_charge = field_str.split()
    return res_name, total_charge, desc

def atom_parser(line):
    """Parse the keyword ATOM"""
    field_str = comment_parser(line)[0]
    key, atom_name, atom_type, charge = field_str.split()
    return atom_name, atom_type, charge

def pairwise_parser(line):
    """Parser for lines with pairwise atoms separated by space"""
    field_str = comment_parser(line)[0]
    # The first field is the keyword
    # followed by atom type in pairs
    fields = field_str.split()[1:]
    # combine each pair into list of tuples
    if len(fields) % 2 == 1:
        raise ValueError(
            f'Odd number of atoms! Cannot group into pairs: {line}'
        )
    atom_pairs = list(zip(fields[::2], fields[1::2]))
    return atom_pairs

def quad_parser(line):
    """Parse IMPR keyword for group of 4 atoms"""
    field_str = comment_parser(line)[0]
    fields = field_str.split()[1:]
    # combine each pair into list of tuples
    if len(fields) % 4 != 0:
        raise ValueError(
            f'Invalid length of topology specification: {line}\n'
            'Multiples  of 4 atoms required'
        )
    ##TODO: use for loop to combine into list of tuples
    quad = (tuple(fields[:4]), tuple(fields[4:]))
    return quad

def octa_parser(line):
    """Parse CMAP keyword for a pair of 4 atoms"""
    field_str = comment_parser(line)[0]
    fields = field_str.split()[1:]
    if len(fields) != 8:
        raise ValueError(
            f'Invalid length of topology specification: {line}\n'
            '8 atoms required'
        )
    # combine into tuple of tuples with shape (2,4)
    octa = (tuple(fields[:4]), tuple(fields[4:]))
    return octa

def ic_parser(line):
    """
    Parse internal coordinate entries. 9 fields total

    fields: [key, I,J,K,L, R(I(J/K)), T(I(JK/KJ)), PHI, T(JKL), R(KL)]
    IJKL are atoms, star (*) on atom K indicate that it is an improper structure
    Strucutres specified for both chain and improper (branch)

        R(I(J/K)): dist of I-J/I-K
        T(I(JK/KJ)): angle of I-J-K/I-K-J
        PHI: Dihedral
        T(JKL): angel of J-K-L
        R(KL): dist of K-L
    """
    field_str = comment_parser(line)[0]
    fields = field_str.split()[1:]
    if len(fields) != 9:
        raise ValueError(
            f'Invalid length of topology specification: {line}\n'
            '9 fields required'
        )
    i, j, k, l, r_ij, t_ijk, phi, t_jkl, r_kl = fields
    if k.startswith('*'):
        is_improper = True
        first_angle = 'R(I-K)'
        dihe = 'T(I-K-J)'
    else:
        is_improper = False
        first_angle = 'R(I-J)'
        dihe = 'T(I-J-K)'

    return (i, j, k.lstrip('*'), l), {
        'improper': is_improper,
        # 'I': i, 'J': j, 'K': k, 'L': l,
        first_angle: float(r_ij), 
        dihe: float(t_ijk),
        'Phi': float(phi),
        'T(J-K-L)': float(t_jkl),
        'R(K-L)': float(r_kl)
    }

def delete_parser(line):
    """
    Parse the keyword DELETE. Keyword has to be followed by the type of the 
    data to be deleted and the value
    """
    field_str = comment_parser(line)[0]
    fields = field_str.split()[1:]
    if len(fields) != 2:
        raise ValueError(
            f'Invalid length of topology specification: {line}\n'
            '2 fields required: key, data'
        )
    return (fields[0], fields[1])

class RTFParser:
    """A parser class to load rtf (residue topology files) into dictionary.
    Parser is initialized with RTF file from file path. If any lines from the 
    file are not parsed, they will be stored in the unparsed_lines"""
    def __init__(self, file_path):
        self.rtf_version = None
        self.topo_dict = None
        self.decl_peptide_atoms = []
        self.default_patchs = {'FIRST':None,'LAST':None}
        self.default_autogen = None
        self.unparsed_lines = []
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [l for l in f.readlines() if not skip_line(l.strip())]
        self._parse_lines(lines)
        n_unparsed = len(self.unparsed_lines)
        if n_unparsed > 0:
            warnings.warn(
                f"Failed to parse {n_unparsed} lines from {file_path}"
            )

    def _parse_lines(self, lines):
        ##TODO: Need to handle keywords DIHE, ANGLE, and PATCH
        self.rtf_version = '.'.join(lines[0].strip().split())
        self.topo_dict = {}
        mass_dict = {}
        for l in lines[1:]:
            if l.startswith('MASS'):
                symbol, mass, desc = mass_parser(l)
                mass_dict.update({symbol: (float(mass), desc)})
            elif l.startswith('DECL'):
                atom = decl_parser(l)
                self.decl_peptide_atoms.append(atom)
            elif l.startswith('DEFA'):
                first_patch, last_patch = defa_parser(l)
                self.default_patchs['FIRST'] = first_patch
                self.default_patchs['LAST'] = last_patch
            elif l.startswith('AUTO'):
                self.default_autogen = auto_parser(l)
            elif l.startswith('RESI') or l.startswith('PRES'):
                res_name, total_charge, desc = resi_parser(l)
                if res_name not in self.topo_dict:
                    # index first group
                    cur_group_i = -1
                    cur_res = {
                        'desc': desc,
                        'total_charge': float(total_charge),
                        'atoms':{
                        },
                        'bonds':{
                            'single':[],
                            'double':[],
                            'triple':[],
                            'aromatic':[]
                        },
                        'impropers':[],
                        'cmap':[],
                        'ic':{},
                        'is_patch': l.startswith('PRES')
                    }
                self.topo_dict.update({res_name: cur_res})
            elif l.startswith('GROUP'):
                # Update group number
                cur_group_i += 1
                cur_atom_group = {cur_group_i: {}}
                cur_res['atoms'].update(cur_atom_group)
            elif l.startswith('ATOM'):
                if cur_group_i == -1:
                    # if no GROUP keyword exist for patch, create a single group
                    cur_group_i = 0
                    cur_atom_group = {cur_group_i: {}}
                    cur_res['atoms'].update(cur_atom_group)
                atom_name, atom_type, atom_charge = atom_parser(l)
                cur_atom_dict = {
                    atom_name: 
                    {
                        'atom_type': atom_type, 
                        'charge': float(atom_charge), 
                        'mass': mass_dict[atom_type][0],
                        'desc': mass_dict[atom_type][1]
                    }
                }
                cur_atom_group[cur_group_i].update(cur_atom_dict)
            elif l.startswith('DONO'):
                if 'H_donors' not in cur_res:
                    cur_res['H_donors'] = []
                donors = tuple(l.split()[1:])
                cur_res['H_donors'].append(donors)
            elif l.startswith('ACCE'):
                if 'H_acceptors' not in cur_res:
                    cur_res['H_acceptors'] = []
                acceptors = tuple(l.split()[1:])
                cur_res['H_acceptors'].append(acceptors)
            elif l.startswith('BOND'):
                single_bonds = pairwise_parser(l)
                cur_res['bonds']['single'].extend(single_bonds)
            elif l.startswith('DOUBLE'):
                double_bonds = pairwise_parser(l)
                cur_res['bonds']['double'].extend(double_bonds)
            elif l.startswith('TRIPLE'):
                triple_bonds = pairwise_parser(l)
                cur_res['bonds']['triple'].extend(triple_bonds)
            elif l.startswith('AROMATIC'):
                aromatic_bonds = pairwise_parser(l)
                cur_res['bonds']['triple'].extend(aromatic_bonds)
            elif l.startswith('IMPR'):
                # Improper (branching structures)
                ##TODO: quad_parser is not correct for IMPR
                impropers = quad_parser(l)
                cur_res['impropers'].extend(impropers)
            elif l.startswith('CMAP'):
                # Dihedral crossterm energy correction map
                cmap = octa_parser(l)
                cur_res['cmap'].append(cmap)
            elif l.startswith('IC') or l.startswith('ic'):
                # Internal Coordinates
                ic_key, ic_param_dict = ic_parser(l)
                cur_res['ic'][ic_key] = ic_param_dict
            elif l.startswith('DELETE'):
                if 'delete' not in cur_res:
                    cur_res['delete'] = []
                delete_entry = delete_parser(l)
                cur_res['delete'].append(delete_entry)
            elif l.startswith('DIHE'):
                pass
            elif l.startswith('ANGLE'):
                pass
            elif l.startswith('PATCH') or l.startswith('patch'):
                pass
            elif l.startswith('END') or l.startswith('end'):
                break
            else:
                self.unparsed_lines.append(l)