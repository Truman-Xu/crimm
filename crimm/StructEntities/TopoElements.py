"""This module defines the TopoEntity class and its subclasses Bond, Angle, Dihedral, and Improper."""
from collections import namedtuple

class TopoEntity:
    """A TopoEntity is a base class for Topology entities such as Bonds, Angles, 
    Dihedrals, and Impropers."""

    def _create_full_id(self):
        """Create unique id by comparing the end atoms to decide if the sequence need to be flipped"""
        atom_ids = tuple(a.full_id for a in self)
        if atom_ids[0] > atom_ids[-1]:
            atom_ids = tuple(reversed(atom_ids))
        return atom_ids
    
    def __getnewargs__(self):
        "Support for pickle protocol 2: http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances"
        return *self, *self.__dict__.values()

    def __getstate__(self):
        """
        Additional support for pickle since parent class implements its own __getstate__
        so pickle does not store or restore the type and order, python 2 problem only
        https://www.python.org/dev/peps/pep-0307/#case-3-pickling-new-style-class-instances-using-protocol-2
        """
        return self.__dict__

    def __hash__(self) -> int:
        return hash(self.full_id)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return self.full_id == other.full_id
    
    def __deepcopy__(self, memo):
        return type(self)(*self, **self.__dict__)

BondTuple = namedtuple('Bond', ['atom1', 'atom2'])
class Bond(TopoEntity, BondTuple):
    """A Bond object represents a bond between two Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 2-element tuple of Atom objects.
    It also has fields that can optionally be used to describe the bond order and type of bond."""

    bond_order_dict = {'single':1, 'double':2, 'triple':3, 'aromatic':2}
    def __new__(cls, atom1, atom2, bond_type=None, param=None):
        """Create a new Bond. """
        bond = super(Bond, cls).__new__(cls, atom1, atom2)
        bond.type = bond_type
        bond.order = cls.bond_order_dict.get(type)
        bond.param = param
        bond.full_id = bond._create_full_id()
        return bond

    def __repr__(self):
        s = "Bond(%s, %s" % (self[0], self[1])
        if self.type is not None:
            s = "%s, type=%s" % (s, self.type)
        if self.order is not None:
            s = "%s, order=%d" % (s, self.order)
        if self.length is not None:
            s = "%s, length=%.2f" % (s, self.length)
        s += ")"
        return s

    def __deepcopy__(self, memo):
        return Bond(*self, self.type, self.param)
    
    @property
    def length(self):
        """return the current bond length"""
        if self[0].coord is None or self[1].coord is None:
            return None
        return (((self[0].coord - self[1].coord)**2).sum())**0.5
    
AngleTuple = namedtuple('Angle', ['atom1', 'atom2', 'atom3'])
class Angle(TopoEntity, AngleTuple):
    """An Angle object represents an angle between three Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 2-element tuple of Bond objects.
    """
    
    def __new__(cls, atom1, atom2, atom3, param=None):
        """Create a new Entity. """
        angle = super(Angle, cls).__new__(cls, atom1, atom2, atom3)
        angle.param = param
        angle.full_id = angle._create_full_id()
        return angle
    
    def __repr__(self):
        s = f"Angle({tuple(self)}, angle = {self.angle:.2f})"
        return s
    
    @property
    def angle(self):
        """return the current angle in degrees"""
        return 0.000
    
DiheTuple = namedtuple('Dihedral', ['i', 'j', 'k', 'l'])
class Dihedral(TopoEntity, DiheTuple):
    """A Dihedral object represents a dihedral angle between four Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 4-element tuple (i, j, k, l) of 
    Atom objects.
    """

    def __new__(cls, atom_i, atom_j, atom_k, atom_l, param=None):
        """Create a new Entity. """
        dihedral = super(Dihedral, cls).__new__(cls, atom_i, atom_j, atom_k, atom_l)
        dihedral.param = param
        dihedral.full_id = dihedral._create_full_id()
        return dihedral

    def __repr__(self):
        s = "Dihedral({}, {}, {}, {}, angle = {:.2f})".format(*self, self.dihe)
        return s

    @property
    def dihe(self):
        """return the current dihedral in degrees"""
        return 0.000
        

ImprTuple = namedtuple('Improper', ['i', 'j', 'k', 'l'])
class Improper(TopoEntity, ImprTuple):
    """A Dihedral object represents a dihedral angle between four Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 4-element tuple (i, j, k, l) of 
    Atom objects.
    """

    def __new__(cls, atom_i, atom_j, atom_k, atom_l, param=None):
        """Create a new Entity. """
        improper = super(Improper, cls).__new__(cls, atom_i, atom_j, atom_k, atom_l)
        improper.param = param
        improper.full_id = improper._create_full_id()
        return improper
    
    def __repr__(self):
        s = "Improper({}, {}, {}, {}, angle = {:.2f})".format(*self, self.impr)
        return s
    
    @property
    def impr(self):
        """return the current improper in degrees"""
        return 0.000