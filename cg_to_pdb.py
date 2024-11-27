from Bio.PDB import PDBParser
import os
import argparse


def load_templates(templates_dir):
    
    parser = PDBParser(QUIET=True)
    templates = {}
    for template_file in os.listdir(templates_dir):
        if template_file.endswith(".pdb"):
            template_name = template_file.split(".")[0]
            template_structure = parser.get_structure(template_name, os.path.join(templates_dir, template_file))
            templates[template_name] = template_structure
    return templates


def combine_pdb_files(cg_file, templates, output_file):
  
    parser = PDBParser(QUIET=True)
    cg_structure = parser.get_structure("CG_Structure", cg_file)

    output_atoms = [] 
    atom_counter = 1  

    added_residues = set()  

    for model in cg_structure:
        for chain in model:
            for residue in chain.get_residues():
                res_name = residue.get_resname()
                res_id = residue.get_id()

           
                if (res_name, res_id) in added_residues:
                    continue

                if res_name in templates:
                    template = templates[res_name]
                    template_chain = next(template.get_chains())

                    for template_residue in template_chain.get_residues():
                        if template_residue.get_id() == res_id:
                            for atom in template_residue.get_atoms():
                                atom.serial_number = atom_counter
                                atom.bfactor = getattr(atom, 'bfactor', 1.0)  
                                atom.element = atom.element or atom.name[0].upper()  
                                atom_counter += 1
                                output_atoms.append(atom)
                else:
                   
                    for atom in residue.get_atoms():
                        atom.serial_number = atom_counter
                        atom.bfactor = getattr(atom, 'bfactor', 1.0)  
                        atom.element = atom.element or atom.name[0].upper() 
                        atom_counter += 1
                        output_atoms.append(atom)

                added_residues.add((res_name, res_id))

   
    save_pdb(output_atoms, output_file)


def save_pdb(atoms, output_file):
    
    with open(output_file, 'w') as output:
        for atom in atoms:
            residue = atom.get_parent()
            chain = residue.get_parent()

            element = atom.element if atom.element else atom.name[0].upper()
            bfactor = atom.bfactor if hasattr(atom, 'bfactor') else 0.0

            output.write(
                "ATOM  {:5d} {:>4} {:>3} {:>1}{:>4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}           {:>2}\n".format(
                    atom.serial_number,
                    atom.name,
                    residue.get_resname(),
                    chain.id,
                    residue.get_id()[1],
                    atom.coord[0],
                    atom.coord[1],
                    atom.coord[2],
                    atom.bfactor,
                    element
                )
            )
    print(f"Plik wyjściowy zapisany jako: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Łączenie plików PDB (gruboziarnisty + szablony).")
    parser.add_argument("--input", required=True, help="Gruboziarnisty plik PDB (np. 430_cg.pdb)")
    parser.add_argument("--output", required=True, help="Plik wyjściowy PDB (np. 430_rec.pdb)")
    parser.add_argument("--templates", required=True, help="Katalog z szablonami PDB (np. ./templates/)")
    args = parser.parse_args()

 
    templates = load_templates(args.templates)

    
    combine_pdb_files(args.input, templates, args.output)
