import sys
from Bio.PDB import *

class SelectResidues(Select):
    """
    Classe pour la sélection de résidus spécifiques.
    """
    def __init__(self, selected_residues):
        self.selected_residues = selected_residues
    
    def accept_residue(self, residue):
        """
        Accepte les résidus sélectionnés.
        """
        if residue in self.selected_residues:
            return True
        else:
            return False

def split_chain(pdb_file, output_file, *args):
    """
    Sélectionne les résidus spécifiés dans un fichier PDB, en utilisant les identifiants de chaîne et les numéros de résidus de début et de fin.
    Enregistre les résidus sélectionnés dans un nouveau fichier PDB.
    """
    # Création du parser PDB
    parser = PDBParser()
    # Lecture du fichier PDB
    structure = parser.get_structure("pdb", pdb_file)
    # Récupération des chaînes spécifiées
    model = structure[0]
    print(model)
    chains = [model[c] for c in args[0::3]]
    print(chains)
    # Sélection des résidus spécifiés pour chaque chaîne
    selected_res = []
    for i, chain in enumerate(chains):
        res_list = Selection.unfold_entities(chain, 'R')
        start_res, end_res = int(args[i*3+1]), int(args[i*3+2])
        selected_res += [res for res in res_list if start_res <= int(res.get_id()[1]) <= end_res]
    # Écriture des résidus sélectionnés dans un nouveau fichier PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, select=SelectResidues(selected_res))

# Appel de la fonction avec les arguments en ligne de commande
if __name__ == "__main__":
    pdb_file = sys.argv[1]
    output_file = sys.argv[-1]
    args = sys.argv[2:-1]
    print(args)
    split_chain(pdb_file, output_file, *args)
