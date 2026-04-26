import json
from typing import TypeAlias, Dict, List
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.MMCIFParser import MMCIFParser


JSON: TypeAlias = dict[str, "JSON"] | list["JSON"] | str | int | float | bool | None
BindingSite: TypeAlias = List[str]
ApoPdbId: TypeAlias = str
HoloPdbId: TypeAlias = str
Chain: TypeAlias = str
Ligand: TypeAlias = tuple[str, str, str]


def load_dataset(path: str) -> JSON:
    """Loads dataset from JSON file.

    Args:
        path (str): Path to JSON file

    Returns:
        JSON: dataset
    """
    with open(path) as f:
        dataset = json.load(f)
    return dataset


def get_apo_binding_residues(dataset: JSON) -> Dict[ApoPdbId, BindingSite]:
    """Loads binding residues for each APO structure from the dataset. As the APO structure can be associated with more than one HOLO structure, 
    you need to loop over all structures to receive every binding residue.

    Returns:
        Dict[ApoPdbId, BindingSite]: Dictionary of all auth_seq_id indices
    """
    apo_residues = {}
    for apo_pdb_id, holo_structures in dataset.items():
        apo_residues[apo_pdb_id] = set()
        for holo_structure in holo_structures:
            apo_residues[apo_pdb_id].update(holo_structure['apo_pocket_selection'])
    return apo_residues


def get_apo_holo_pairs(dataset: JSON) -> set[tuple[ApoPdbId, HoloPdbId]]:
    """Retrieves every apo-holo pair from the dataset

    Args:
        dataset (JSON): dataset

    Returns:
        set[tuple[ApoPdbId, HoloPdbId]]: every apo-holo pair
    """
    apo_holo_pairs = set()
    for apo_pdb_id, holo_structures in dataset.items():
        for holo_structure in holo_structures:
            holo_pdb_id = holo_structure['holo_pdb_id'] 
            apo_holo_pairs.add((apo_pdb_id, holo_pdb_id))
    return apo_holo_pairs


def extract_apo_chains(dataset: JSON) -> Dict[ApoPdbId, Chain]:
    """Extracts all chains for each apo structure from the dataset. As the APO structure can be associated with more than one HOLO structure, 
    you need to loop over all structures to receive every chain.

    Args:
        dataset (JSON): dataset

    Returns:
        Dict[ApoPdbId, Chain]: Dictionary of apo pdb ids and their corresponding chains
    """
    apo_chains = {}
    for apo_pdb_id, holo_structures in dataset.items():
        for holo_structure in holo_structures:
            chain = holo_structure['apo_chain']
            apo_chains[apo_pdb_id] = chain
            break

    return apo_chains


def get_sequence_with_auth_and_label_indices(crypto_path: str, pdb_id: ApoPdbId, chain_id: Chain) -> tuple[str, list[str], list[str]]:
    '''
    Extracts the sequence of the specified chain from the cif file and returns it along with the auth_seq_id and label_seq_id indices.

    In PDB files, the author-assigned residue numbers (`auth_seq_id`)  may differ from the sequential numbering assigned by PDB (`label_seq_id`). 
    In CryptoBench, binding residues are identified using author-defined residue numbers.
    '''

    parser = MMCIFParser(QUIET=True)
    auth_structure = parser.get_structure(pdb_id, f"{crypto_path}/cryptobench-dataset/auxiliary-data/cif-files/{pdb_id.lower()}.cif")

    parser = MMCIFParser(auth_residues=False, QUIET=True)
    label_structure = parser.get_structure(pdb_id, f"{crypto_path}/cryptobench-dataset/auxiliary-data/cif-files/{pdb_id.lower()}.cif")

    sequence = ''
    label_indices = []
    auth_indices = []

    for auth_residue, label_residue in zip(auth_structure[0][chain_id].get_residues(), label_structure[0][chain_id].get_residues()):
        if auth_residue.get_id()[0][0] == ' ':
            sequence += protein_letters_3to1[auth_residue.get_resname().title()]
            label_indices.append(str(label_residue.get_id()[1]))
            auth_indices.append(str(auth_residue.get_id()[1]))

    return sequence, label_indices, auth_indices