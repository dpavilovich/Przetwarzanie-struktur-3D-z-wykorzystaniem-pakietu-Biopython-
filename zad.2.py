from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
import numpy as np
import matplotlib.pyplot as plt


def load_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure


def calculate_phi_psi(structure):
    phi_psi_values = []
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())  
            for i in range(1, len(residues) - 1):  
                try:
                  
                    n_prev = residues[i - 1]['N'].get_vector()
                    c_alpha = residues[i]['CA'].get_vector()
                    c = residues[i]['C'].get_vector()
                    n_next = residues[i + 1]['N'].get_vector()
                    
                  
                    phi = calc_dihedral(n_prev, c_alpha, c, n_next) * (180 / np.pi)  
                    
                
                    ca_next = residues[i + 1]['CA'].get_vector()
                    psi = calc_dihedral(c_alpha, c, n_next, ca_next) * (180 / np.pi)
                    
                    phi_psi_values.append((phi, psi))
                except KeyError:
                  
                    continue
    return phi_psi_values


def plot_ramachandran(phi_psi_values):
    phi, psi = zip(*phi_psi_values) 
    plt.figure(figsize=(8, 8))
    plt.scatter(phi, psi, s=10, alpha=0.7, color='blue')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.title("Wykres Ramachandrana")
    plt.xlabel("Kąt phi (°)")
    plt.ylabel("Kąt psi (°)")
    plt.grid(color='lightgray', linestyle='--', linewidth=0.5)
    plt.show()


pdb_file = "C:/politechnika/semestr 5 bioinformatyka/bioinformatyka strukturalna-lab/4ywo (2).pdb" 


structure = load_structure(pdb_file)
phi_psi_values = calculate_phi_psi(structure)
plot_ramachandran(phi_psi_values)
