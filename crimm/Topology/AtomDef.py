from copy import copy
import crimm.StructEntities as Entities
from crimm.Data.atom_types import atom_type_dict

class AtomTypeDefinition:
    """Atom type definition class. This class is used to define the atom type and
    other properties of an atom. It is also used to supplement the AtomDefinition 
    class."""
    def __init__(self, atom_type, mass, element, desc):
        self.atom_type = atom_type
        self.mass = mass
        self.desc = desc
        self.element = element

    def __repr__(self):
        repr_str = f"<Atom Type Definition {self.atom_type}>"
        return repr_str


class AtomDefinition:
    """Atom definition class. This class is used to define the specific atom for 
    a residue topology with atom type and other properties. It is used to create
    ResidueDefinition class."""
    def __init__(
            self, name, atom_type_def, parent_res_def, charge
        ):
        self.atom_type_def = atom_type_def
        self.parent_def = parent_res_def
        self.name = name
        self.is_donor = False
        self.is_acceptor = False
        self.charge = charge

    def __repr__(self):
        atom_type = self.atom_type_def.atom_type
        repr_str = f"<Atom Definition name={self.name}, type={atom_type}>"
        return repr_str
    
    def create_new_atom(self, coords = None, serial_number = 0):
        """Create a new atom instance from the Atom definition. The default coordinates
        will be none if not specified, and a MissingAtom instance will be created."""

        return Entities.Atom(
            name = self.name,
            coord=coords,
            bfactor=0.0,
            occupancy=1.0,
            altloc=' ',
            serial_number=serial_number,
            element=self.atom_type_def.element,
            fullname=self.name,
            topo_definition=self,
            mass=self.atom_type_def.mass,
        )


atom_type_definitions = {}
for entity_type, cur_atom_type_dict in atom_type_dict.items():
    cur_def_dict = {}
    atom_type_definitions[entity_type] = cur_def_dict
    for atom_type_code, atom_dict in cur_atom_type_dict.items():
        cur_def_dict[atom_type_code] = AtomTypeDefinition(
            atom_type = atom_type_code,
            mass = atom_dict['mass'],
            element = atom_dict['element'],
            desc = atom_dict['description']
        )