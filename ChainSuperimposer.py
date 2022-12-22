
from Bio.Align import PairwiseAligner
from Bio.PDB import Superimposer
from Bio.PDB.Chain import Chain
from PeptideChain import PeptideChain, StandardChain
import warnings

class ChainSuperimposer(Superimposer):
    
    def __init__(self) -> None:
        super().__init__()

    def find_all_common_res(self, 
        ref_chain: PeptideChain | StandardChain, 
        mov_chain: PeptideChain | StandardChain):
        """ Find all common residues from two chains """

        # if any of the chains is not PeptideChain class, convert them
        ## FIXME: change to CanonicalPeptideChain for both
        is_ref_standard = isinstance(ref_chain, StandardChain)
        is_mov_standard = isinstance(mov_chain, StandardChain)

        if is_ref_standard and is_mov_standard:
            ref_aligned, mov_aligned = self._find_common_res_for_standard_chain(ref_chain, mov_chain)
        else:
            ref_aligned, mov_aligned = self._find_common_res_for_simple_chain(ref_chain, mov_chain)       

        assert len(ref_aligned) == len(mov_aligned)

        return ref_aligned, mov_aligned
         
    def _find_common_res_for_standard_chain(
            self, ref_chain: StandardChain, mov_chain: StandardChain
        ):
        if ref_chain.can_seq == mov_chain.can_seq:
            # If canonical sequences are present and identical
            # Use all present residues from both chains
            r_present = self._get_present_res_ids_from_sele(ref_chain)
            m_present = self._get_present_res_ids_from_sele(mov_chain)
            # Get the lists of residues handles that are present in both chains
            # for superposition
            ref_aligned, mov_aligned = self._get_common_res(
                ref_chain, mov_chain, r_present, m_present
            )
        else:
            ## TODO: TEST THIS!!
            # Canonical sequence is present for both but they are not identical
            # Align them first to find the common segments
            ref_range, mov_range = self._get_aligned_ranges(ref_chain.can_seq, 
                                                            mov_chain.can_seq)
            # Create a set of the common residue ids from alignment
            r_sele_res = self._get_aligned_res_ids(1, ref_range)
            # Get residues ids that are present in these common residues
            r_present = self._get_present_res_ids_from_sele(ref_chain, r_sele_res)
            m_sele_res = self._get_aligned_res_ids(1, mov_range)
            m_present = self._get_present_res_ids_from_sele(mov_chain, m_sele_res)
            # Get the lists of residue handles that are present in both chains 
            # for superposition
            ref_aligned, mov_aligned = self._get_common_res(
                ref_chain, mov_chain, r_present, m_present
            )
        return ref_aligned, mov_aligned

    def _find_common_res_for_simple_chain(
            self, ref_chain: PeptideChain, mov_chain: PeptideChain
        ):
        ref_seq = ref_chain.extract_all_seq()
        mov_seq = mov_chain.extract_all_seq()
        if ref_seq != mov_seq:
            # No canonical sequence. Fall back to present sequence.
            # Align to find the common residue first
            ref_range, mov_range = self._get_aligned_ranges(ref_seq, mov_seq)
            # Since no canonical sequence, all residues from aligned 
            # sequence is used as common residues for superposition
            ref_aligned = self._get_aligned_res(ref_chain, ref_range)
            mov_aligned = self._get_aligned_res(mov_chain, mov_range) 
        else:
            # If the sequences are identical, no alignment will be performed,
            # and superposition will be attempted on all residues
            ref_aligned = ref_chain.child_list
            mov_aligned = mov_chain.child_list
        return ref_aligned, mov_aligned

    def _get_aligned_res_ids(self, start_id, ranges):
        """
        Get a set of residue ids from the aligned ranges
        """
        aligned_res = []
        for start, end in ranges:
            # Res ids starts at different location than the aligned range ids
            # no hetflags included
            cur_id_range = [
                (' ', i, ' ') for i in range(start+start_id,end+start_id)
            ]
            aligned_res.extend(cur_id_range)
        return set(aligned_res)

    def _get_present_res_ids_from_sele(
            self, chain: PeptideChain, selected_res_id: set = set()
        ):
        """
        Get a set of residue ids present in chain from a set of 
        selected res ids. Default to select all present residues
        """
        present = set(chain.child_dict.keys())
        if selected_res_id:
            # TODO: hetflag not included now, include the ones with hetflags
            return selected_res_id.intersection(present)
        return present

    def _get_common_res(
            self, ref_chain, mov_chain, ref_present: set, mov_present: set
        ):

        common = sorted(ref_present.intersection(mov_present))
        ref_common = []
        mov_common = []
        for i in common:
            # i here corresponds to residue id now
            ref_common.append(ref_chain[i])
            mov_common.append(mov_chain[i])
        return ref_common, mov_common

    def _get_aligned_ranges(self, ref_seq, mov_seq):
        aligner = PairwiseAligner()
        alignments = aligner.align(ref_seq, mov_seq)
        top_alignment = alignments[0]
        ref_algn_ranges, mov_algn_ranges = top_alignment.aligned
        return ref_algn_ranges, mov_algn_ranges

    def _get_aligned_res(self, chain, ranges):
        aligned_res = []
        for start, end in ranges:
            for i in range(start, end):
                aligned_res.append(chain.child_list[i])
        return aligned_res
    
    def check_chain_type(self, chain):
        if (
            isinstance(chain, PeptideChain) or \
            isinstance(chain, StandardChain)
        ):
            return chain
        elif isinstance(chain, Chain):
            return PeptideChain(chain)
        else:
            raise TypeError(
                'Chain, PeptideChain, or StandardChain class is required to use'
                'ChainImposer for superposition.'
                )
        
    def set_chains(self, ref_chain, mov_chain, on_atoms = 'CA'):
        
        ref_chain = self.check_chain_type(ref_chain)
        mov_chain = self.check_chain_type(mov_chain)
        ref_common, mov_common = self.find_all_common_res(ref_chain, mov_chain)

        if on_atoms == 'CA':
            ref_align_atoms, mov_align_atoms = \
                    self._find_common_CA_atoms(ref_common, mov_common)
        elif on_atoms == 'backbone':
            ref_align_atoms, mov_align_atoms = \
                    self._find_common_backbone_atoms(ref_common, mov_common)
        elif on_atoms == 'all':
            ref_align_atoms, mov_align_atoms = \
                    self._find_all_common_atoms(ref_common, mov_common)
        else:
            raise ValueError('on_atoms has to be selected from '
                '{"backbone", "all", "CA"}') 
        self.set_atoms(ref_align_atoms, mov_align_atoms)

    def _find_all_common_atoms(self, ref_res_list, mov_res_list):
        """Find all common atoms from two lists of residues that are aligned"""
        ref_align_atoms = []
        mov_align_atoms = []
        for ref_res, mov_res in zip(ref_res_list, mov_res_list):
            for atom_id in ref_res.child_dict.keys():
                # Make sure move chain residues have the same atoms
                # Otherwise, the missing atom will be excluded
                if atom_id in mov_res:
                    ref_align_atoms.append(ref_res[atom_id])
                    mov_align_atoms.append(mov_res[atom_id])

        return ref_align_atoms, mov_align_atoms

    def _find_common_CA_atoms(self, ref_res_list, mov_res_list):
        ref_align_atoms = []
        mov_align_atoms = []
        for ref_res, mov_res in zip(ref_res_list, mov_res_list):
            # Make sure both residues have the CA atom
            # Otherwise, the missing atom will be excluded
            if 'CA' in ref_res and 'CA' in mov_res:
                ref_align_atoms.append(ref_res['CA'])
                mov_align_atoms.append(mov_res['CA'])

        return ref_align_atoms, mov_align_atoms
    
    def _find_common_backbone_atoms(self, ref_res_list, mov_res_list):
        ref_align_atoms = []
        mov_align_atoms = []
        for ref_res, mov_res in zip(ref_res_list, mov_res_list):
            for atom_id in ['N','CA','C','O']:
                # Make sure both residues have the backbone atom
                # Otherwise, the missing atom will be excluded
                if atom_id in ref_res and atom_id in mov_res:
                    ref_align_atoms.append(ref_res[atom_id])
                    mov_align_atoms.append(mov_res[atom_id])

        return ref_align_atoms, mov_align_atoms

    def set_around_gap(self, ref_chain, mov_chain, start, end, cutoff=10):
        """
        Set atoms for superposition around a gap/missing residues on the ref 
        chain. Start and end residue index needed. Alignment will be performed on
        backbone atoms.
        """
        ## TODO: add options for on_atoms, and find common res list first 
        ## then run the for-loop on those.

        if hasattr(ref_chain,'can_seq'):
            ref_len = len(ref_chain.can_seq)
        else:
            ref_len = len(ref_chain.child_list)

        if hasattr(mov_chain,'can_seq'):
            mov_len = len(mov_chain.can_seq)
        else:
            mov_len = len(mov_chain.child_list)

        seq_len = min(ref_len, mov_len)
        res_ids = []
        # select the residue ids around the gap by a cutoff
        if start < cutoff and end < (seq_len - 2*cutoff):
            # Case of missing N terminal and the cutoff ends before the 
            # last residue of the chain
            res_ids = list(range(end, end+2*cutoff))
        elif end > seq_len - cutoff and start > 2*cutoff:
            # Case of missing C terminal and the cutoff ends before the
            # first residue of the chain
            res_ids = list(range(start-2*cutoff))
        elif start > cutoff and end < (seq_len - cutoff):
            # The case of an actual gap in the chain and
            # the cutoff ends before the first and last residue of the chain
            res_ids = list(range(start-cutoff, start))+list(range(end, end+cutoff))
        else:
            # In other cases, the specified cutoff runs 
            warnings.warn('The cutoff value extends out of the available '
                    'residues in the chain. All backbone atoms on both chains '
                    'are selected for superposition')
            ref_common, mov_common = \
                    self.find_all_common_res(ref_chain, mov_chain)
            ref_atoms, mov_atoms = \
                    self._find_common_backbone_atoms(ref_common, mov_common) 
            self.set_atoms(ref_atoms, mov_atoms) 
            return

                                   
        ref_atoms = []
        mov_atoms = []
        for i in res_ids:
            if i in mov_chain and i in ref_chain:
                if mov_chain[i].resname == ref_chain[i].resname:
                    for atom in ['N','CA','C','O']:
                        if atom in mov_chain[i]:
                            ref_atoms.append(ref_chain[i][atom])
                            mov_atoms.append(mov_chain[i][atom])
        
        self.set_atoms(ref_atoms, mov_atoms)
    
    def apply_transform(self, chain: Chain):
        self.apply(chain.get_atoms())

    def show(self, ref_chain, mov_chain):

        from IPython.display import display
        try:
            import nglview as nv
        except ImportError:
            raise ImportError(
                "WARNING: nglview not found! Install nglview to show"
                "protein structures."
                "http://nglviewer.org/nglview/latest/index.html#installation"
            )
        # if any of the chains is not PeptideChain class, convert them
        if not isinstance(ref_chain, PeptideChain):
            ref_chain = PeptideChain(ref_chain)
        if not isinstance(mov_chain, PeptideChain):
            mov_chain = PeptideChain(mov_chain)

        view = nv.NGLWidget()
        ref_chain.load_nglview(view)
        mov_chain.load_nglview(view)
        
        display(view)