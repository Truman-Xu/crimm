
# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Consumer class that builds a Structure object.
This is used by the PDBParser and MMCIFparser classes.
"""

import warnings

# SMCRA hierarchy
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning, PDBConstructionException
from PeptideChain import PeptideChain, StandardChain

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Atom import Atom, DisorderedAtom
from collections import namedtuple
class ChmStructureBuilder(StructureBuilder):
    """Deals with constructing the Structure object.
    The StructureBuilder class is used by the PDBParser classes to
    translate a file to a Structure object.
    """

    def __init__(self):
        """Initialize the class."""
        super().__init__()

    # Public methods called by the Parser classes

    def init_structure(self, structure_id):
        """Initialize a new Structure object with given id.

        Arguments:
         - id - string

        """
        self.structure = Structure(structure_id)

    def init_chain(self, chain_id):
        """Create a new Chain object with given id.

        Arguments:
         - chain_id - string

        """
        if self.model.has_id(chain_id):
            self.chain = self.model[chain_id]
            warnings.warn(
                "WARNING: Chain %s is discontinuous at line %i."
                % (chain_id, self.line_counter),
                PDBConstructionWarning,
            )
        else:
            self.chain = PeptideChain(chain_id)
            self.model.add(self.chain)

    def init_schain(self, chain_id, chain_info: namedtuple):
        if self.model.has_id(chain_id):
            self.chain = self.model[chain_id]
            warnings.warn(
                "WARNING: Chain %s is discontinuous at line %i."
                % (chain_id, self.line_counter),
                PDBConstructionWarning,
            )
        else:
            self.chain = StandardChain(chain_id, **chain_info._asdict())
            self.model.add(self.chain)

    def reset_disordered_res(self,disordered_residue):
        for resname, child_residue in disordered_residue.child_dict.items():
            if child_residue.id == disordered_residue.id:
                disordered_residue.disordered_select(resname)

    def finish_chain_construction(self):
        if not hasattr(self, 'chain'):
            return
        for res in self.chain:
            if isinstance(res, DisorderedResidue):
                self.reset_disordered_res(res)
        self.chain.update()
            
    def init_seg(self, segid):
        """Flag a change in segid.

        Arguments:
         - segid - string

        """
        self.segid = segid

    def process_duplicate_res(self, res_id, resname):
        field, resseq, icode = res_id
        # There already is a residue with the id (field, resseq, icode).
        # This only makes sense in the case of a point mutation.
        warnings.warn(
            "WARNING: Residue ('%s', %i, '%s') redefined at line %i."
            % (field, resseq, icode, self.line_counter),
            PDBConstructionWarning,
        )
        duplicate_residue = self.chain[res_id]
        self.process_disordered_res(duplicate_residue, res_id, resname)

    def process_duplicate_het(self, res_id, resname):
        het_res_id = self.chain.find_het_by_seq(res_id[1])[0]
        duplicate_residue = self.chain[het_res_id]
        self.process_disordered_res(duplicate_residue, res_id, resname)

    def process_disordered_res(self, duplicate_residue, res_id, resname):
        field, resseq, icode = res_id
        if duplicate_residue.is_disordered() == 2:
            # The residue in the chain is a DisorderedResidue object.
            # So just add the last Residue object.
            if duplicate_residue.disordered_has_id(resname):
                # The residue was already made
                self.residue = duplicate_residue
                duplicate_residue.disordered_select(resname)
                return
            
            # Make a new residue and add it to the already
            # present DisorderedResidue
            new_residue = Residue(res_id, resname, self.segid)
            duplicate_residue.disordered_add(new_residue)
            self.residue = duplicate_residue
            return
        
        if resname == duplicate_residue.resname:
            # Not disordered but resname and id already exist
            warnings.warn(
                "WARNING: Residue ('%s', %i, '%s','%s') already defined "
                "with the same name at line  %i."
                % (field, resseq, icode, resname, self.line_counter),
                PDBConstructionWarning,
            )
            self.residue = duplicate_residue
            return
        
        # Make a new DisorderedResidue object and put all
        # the Residue objects with the id (field, resseq, icode) in it.
        # These residues each should have non-blank altlocs for all their atoms.
        # If not, the PDB file probably contains an error.
        if not self._is_completely_disordered(duplicate_residue):
            # if this exception is ignored, a residue will be missing
            self.residue = None
            raise PDBConstructionException(
                "Blank altlocs in duplicate residue %s ('%s', %i, '%s')"
                % (resname, field, resseq, icode)
            )
        
        self.chain.detach_child(duplicate_residue.id)
        new_residue = Residue(res_id, resname, self.segid)
        disordered_residue = DisorderedResidue(duplicate_residue.id)
        self.chain.add(disordered_residue)
        disordered_residue.disordered_add(duplicate_residue)
        disordered_residue.disordered_add(new_residue)
        self.residue = disordered_residue

    def init_residue(self, resname, field, resseq, icode):
        """Create a new Residue object.

        Arguments:
         - resname - string, e.g. "ASN"
         - field - hetero flag, "W" for waters, "H" for
           hetero residues, otherwise blank.
         - resseq - int, sequence identifier
         - icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)
        
        if self.chain.has_id(res_id):
            self.process_duplicate_res(res_id, resname)
        elif len(self.chain.find_het_by_seq(resseq)) > 0:
            self.process_duplicate_het(res_id, resname)
        else:
            self.residue = Residue(res_id, resname, self.segid)
            self.chain.add(self.residue)

    def init_atom(
        self,
        name,
        coord,
        b_factor,
        occupancy,
        altloc,
        fullname,
        serial_number=None,
        element=None,
        pqr_charge=None,
        radius=None,
        is_pqr=False,
    ):
        """Create a new Atom object.

        Arguments:
         - name - string, atom name, e.g. CA, spaces should be stripped
         - coord - Numeric array (Float0, size 3), atomic coordinates
         - b_factor - float, B factor
         - occupancy - float
         - altloc - string, alternative location specifier
         - fullname - string, atom name including spaces, e.g. " CA "
         - element - string, upper case, e.g. "HG" for mercury
         - pqr_charge - float, atom charge (PQR format)
         - radius - float, atom radius (PQR format)
         - is_pqr - boolean, flag to specify if a .pqr file is being parsed

        """
        residue = self.residue
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_id(name):
            duplicate_atom = residue[name]
            # atom name with spaces of duplicate atom
            duplicate_fullname = duplicate_atom.get_fullname()
            if duplicate_fullname != fullname:
                # name of current atom now includes spaces
                name = fullname
                warnings.warn(
                    "Atom names %r and %r differ only in spaces at line %i."
                    % (duplicate_fullname, fullname, self.line_counter),
                    PDBConstructionWarning,
                )
        if not is_pqr:
            self.atom = Atom(
                name,
                coord,
                b_factor,
                occupancy,
                altloc,
                fullname,
                serial_number,
                element,
            )
        elif is_pqr:
            self.atom = Atom(
                name,
                coord,
                None,
                None,
                altloc,
                fullname,
                serial_number,
                element,
                pqr_charge,
                radius,
            )
        if altloc != " ":
            # The atom is disordered
            if residue.has_id(name):
                # Residue already contains this atom
                duplicate_atom = residue[name]
                if duplicate_atom.is_disordered() == 2:
                    duplicate_atom.disordered_add(self.atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.detach_child(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add(disordered_atom)
                    disordered_atom.disordered_add(self.atom)
                    disordered_atom.disordered_add(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn(
                        "WARNING: disordered atom found with blank altloc before "
                        "line %i.\n" % self.line_counter,
                        PDBConstructionWarning,
                    )
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.disordered_add(self.atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add(self.atom)

    def set_symmetry(self, spacegroup, cell):
        """Set symmetry."""
        pass