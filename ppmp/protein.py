from Bio import PDB
import json


class Protein:
    """Instance is one protein"""

    # Initialize modules length, read into memory.
    modules_path = './data/modules-length'
    with open(modules_path, 'r') as f:
        modules_length = dict()
        for line in f:
            words = line.split()
            modules_length[words[0]] = int(words[1])

    # Set up biopyhton PDB parser
    parser = PDB.PDBParser()
    # Set up biopython Superimposer
    sup = PDB.Superimposer()

    @staticmethod
    def kabsch(q, p, modules_range=None):
        """Calculate RMSD between parts of two proteins using the Kabsch method.

        :param q: First instance of Protein.
        :param p: Second instance of Protein.
        :param modules_range: Range of the modules that should be considered.
        If not given take the whole proteins.
        :return: RMSD between the proteins for the given range.
        """
        assert len(q.residue_chain) == len(p.residue_chain)

        if modules_range is None:
            Protein.sup.set_atoms(q.residue_chain, p.residue_chain)
            return Protein.sup.rms
        else:
            modules_rmsd = []
            for i in modules_range:
                start = sum(q.modules_sections[:i])
                end = start + q.modules_sections[i]
                Protein.sup.set_atoms(q.residue_chain[start:end],
                                      p.residue_chain[start:end])
                modules_rmsd.append(Protein.sup.rms)
            return modules_rmsd

    def __init__(self, structure_name, pdb_path, json_path, strict=True):
        self.name = structure_name
        self.path = pdb_path
        self.structure = Protein.parser.get_structure(self.name, self.path)

        with open(json_path, 'r') as f:
            self.modules_chain = json.load(f)['nodes']

        self.modules_sections = [Protein.modules_length[x + '.pdb'] for x in self.modules_chain]

        self.residue_chain = list()
        for residue in self.structure.get_residues():
            self.residue_chain.append(residue['CA'])

        if strict:
            assert sum(self.modules_sections) == len(self.residue_chain)
