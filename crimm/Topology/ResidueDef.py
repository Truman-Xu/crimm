import warnings
from typing import List, Dict, Tuple
from collections import OrderedDict
import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
import crimm.StructEntities as Entities
from crimm.Topology.AtomDef import AtomDefinition, atom_type_definitions

RTF_VERSION = 'c36'
aa_3to1 = protein_letters_3to1_extended.copy()
aa_3to1.update({'HSE':'H', 'HSD':'H', 'HSP':'H'})
na_3to1 = {
    'GUA':'G', 'ADE':'A', 'CYT':'C', 'THY':'T', 'URA':'U'
}
na_1to3 = {
    'A': 'ADE', 'C': 'CYT', 'G': 'GUA', 'T': 'THY', 'U': 'URA',
}
bond_order_dict = {'single':1, 'double':2, 'triple':3, 'aromatic':2}

class ResidueDefinition:
    def __init__(self) -> None:
        self.resname = None
        self.restype = None
        self.atom_groups = []
        self.atom_dict = {}
        self.bonds = []
        self.impropers = []
        self.cmap = []
        self.H_donors = []
        self.H_acceptors = []
        self.description = None
        self.ic = {}

    def __repr__(self) -> str:
        return f'<ResidueDef {self.resname}>'
    
    def __str__(self) -> str:
        return str(self.resname)
    
    def __len__(self) -> int:
        return len(self.atom_dict)
    
    def __iter__(self):
        return iter(self.atom_dict.values())
    
    def __getitem__(self, key):
        return self.atom_dict[key]
    
    def __contains__(self, key):
        return key in self.atom_dict
    
    def get_atom_defs(self, atom_name, atom_type, charge):
        atom_type_def = atom_type_definitions[self.restype][atom_type]
        return AtomDefinition(atom_name, atom_type_def, self, charge)

class Alanine(ResidueDefinition):
    def __init__(self):
        super().__init__()
        self.resname = 'ALA'
        self.restype = 'protein'
        self.atom_groups = [
            ('N', 'HN', 'CA', 'HA'), ('CB', 'HB1', 'HB2', 'HB3'), ('C', 'O')
        ]
        self.atom_dict = {
            'N': self.get_atom_defs('N', 'NH1', -0.47),
            'HN': self.get_atom_defs('HN', 'H', 0.31),
            'CA': self.get_atom_defs('CA', 'CT1', 0.07),
            'HA': self.get_atom_defs('HA', 'HB1', 0.09),
            'CB': self.get_atom_defs('CB', 'CT3', -0.27),
            'HB1': self.get_atom_defs('HB1', 'HA', 0.09),
            'HB2': self.get_atom_defs('HB2', 'HA', 0.09),
            'HB3': self.get_atom_defs('HB3', 'HA', 0.09),
            'C': self.get_atom_defs('C', 'C', 0.51),
            'O': self.get_atom_defs('O', 'O', -0.51)
        }
        self.bonds = {
            'single':(
                ('CA', 'CB'), ('N', 'HN'), ('N', 'CA'), 
                ('CA', 'C'), ('C', '+N'), ('CA', 'HA'), ('CB', 'HB1'), 
                ('CB', 'HB2'), ('CB', 'HB3')
            ),
            'double':(('O','C'),),
            'triple':(),
            'aromatic':()
        }
        self.impropers = (('N', '-C', 'CA', 'HN'), ('C', 'CA', '+N', 'O'))
        self.cmap = (('-C', 'N', 'CA', 'C'), ('N', 'CA', 'C', '+N'))
        self.H_donors = (('HN', 'N'), )
        self.H_acceptors = (('O', 'C'), )
