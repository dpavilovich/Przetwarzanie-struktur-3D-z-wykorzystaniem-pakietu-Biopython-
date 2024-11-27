from Bio import PDB


def convert_to_coarse_grained(pdb_file, output_file):
   
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_file)
    
  
    io = PDB.PDBIO()

   
    coarse_atoms = set([
        'N9', 'C2', 'C6',  
        'N1', 'C2', 'C4',  
        'P', 'C4'          
    ])
    
    atoms_to_write = []

   
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    
                    if atom.get_name() in coarse_atoms:
                        atoms_to_write.append(atom)
    
   
    class CoarseGrainedSelect(PDB.Select):
        def accept_atom(self, atom):
            
            return atom in atoms_to_write
    
   
    io.set_structure(structure)

   
    with open(output_file, 'w') as output:
        for atom in atoms_to_write:
           
            residue = atom.get_parent()
            chain = residue.get_parent()
            
            
            output.write(
                "ATOM  {:5d}  {:<3}   {:<1} {:<1}{:>4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}           {:<2}\n".format(
                    atom.get_serial_number(),
                    atom.get_name(),
                    residue.get_resname(),
                    chain.get_id(),
                    residue.get_id()[1],
                    atom.get_coord()[0],
                    atom.get_coord()[1],
                    atom.get_coord()[2],
                    atom.bfactor,
                    atom.element
                )
            )


convert_to_coarse_grained('430D.pdb', '430D_coarse.pdb.pdb')
