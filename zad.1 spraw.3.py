from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt


def calc_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))


def load_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure


def generate_contact_map(structure, distance_threshold=8.0):
  
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue: 
                    ca_coords.append(residue['CA'].get_coord())
    
    n = len(ca_coords)
    contact_map = np.zeros((n, n), dtype=int)


    for i in range(n):
        for j in range(i + 1, n):  
            distance = calc_distance(ca_coords[i], ca_coords[j])
            if distance <= distance_threshold:
                contact_map[i, j] = 1
                contact_map[j, i] = 1  
    
    return contact_map


def plot_contact_map(contact_map):
    plt.figure(figsize=(8, 8))
    plt.imshow(contact_map, cmap='Greys', interpolation='nearest')
    plt.title("Mapa kontaktów")
    plt.xlabel("Reszta 1")
    plt.ylabel("Reszta 2")
    plt.colorbar(label="Kontakt (0/1)")
    plt.show()


pdb_file = "C:/politechnika/semestr 5 bioinformatyka/bioinformatyka strukturalna-lab/4ywo (2).pdb"  
structure = load_structure(pdb_file)

contact_map = generate_contact_map(structure, distance_threshold=8.0)


plot_contact_map(contact_map)
