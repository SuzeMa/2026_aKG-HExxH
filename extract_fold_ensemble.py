import os
import glob
import pymol
from pymol import cmd

try:
    from PIL import Image
except ImportError:
    print("Error: The 'Pillow' library is required to save TIFF files.")
    exit(1)

TEMPLATE_PDB = "aKG-HExxH_core-fold.pdb"
TARGET_DIR = "./alphafold_TIGR04267"               
OUTPUT_DIR = "./extracted_corefolds"     
ENSEMBLE_PDB = "ensemble_corefolds.pdb"  
OUTPUT_FASTA = "ensemble_corefolds_sequences.fasta.txt" 
OUTPUT_TIFF = "ensemble_publishable.tiff"

def get_clean_resis(selection):
    resis = set()
    cmd.iterate(f"({selection}) and name CA", "resis.add(resi)", space={'resis': resis})
    valid_resis = []
    for r in resis:
        try:
            valid_resis.append(int(''.join(filter(str.isdigit, r))))
        except ValueError:
            pass
    return sorted(valid_resis)

def main():
    pymol.pymol_argv = ['pymol', '-qc'] 
    pymol.finish_launching()

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"--- Stage 1: Topology-Anchored Boundary Extension & Extraction ---")
    cmd.load(TEMPLATE_PDB, "template")
    
    temp_true_resis = get_clean_resis("template")
    if not temp_true_resis:
        print("Error: Could not read template residues.")
        return
    temp_true_min = temp_true_resis[0]
    temp_true_max = temp_true_resis[-1]
    
    target_files = glob.glob(os.path.join(TARGET_DIR, "*.pdb"))
    ensemble_models = []
    fasta_sequences = []

    print(f"Found {len(target_files)} target PDB files.")

    for target_file in target_files:
        input_name = os.path.splitext(os.path.basename(target_file))[0]
        cmd.load(target_file, "target_struct")
        chains = cmd.get_chains("target_struct")

        for chain in chains:
            chain_obj = f"target_chain_{chain}"
            cmd.create(chain_obj, f"target_struct and chain {chain}")

            try:
                cmd.super(chain_obj, "template", object="aln_obj")

                target_matched_resis = get_clean_resis(f"{chain_obj} and aln_obj")
                temp_matched_resis = get_clean_resis("template and aln_obj")
                
                if not target_matched_resis or not temp_matched_resis:
                    cmd.delete(chain_obj)
                    continue
                
                target_matched_min = target_matched_resis[0]
                target_matched_max = target_matched_resis[-1]
                
                temp_matched_min = temp_matched_resis[0]
                temp_matched_max = temp_matched_resis[-1]

                n_extension = temp_matched_min - temp_true_min
                c_extension = temp_true_max - temp_matched_max

                target_final_min = target_matched_min - n_extension
                target_final_max = target_matched_max + c_extension

                extracted_obj = f"{input_name}_chain_{chain}_corefold"
                cmd.create(extracted_obj, f"{chain_obj} and resi {target_final_min}-{target_final_max}")
                
                cmd.super(extracted_obj, "template")

                out_path = os.path.join(OUTPUT_DIR, f"{input_name}_chain_{chain}_corefold.pdb")
                cmd.save(out_path, extracted_obj)
                
                fasta_str = cmd.get_fastastr(extracted_obj)
                fasta_sequences.append(fasta_str)

                ensemble_models.append(extracted_obj)
                cmd.delete("aln_obj")

            except pymol.CmdException:
                pass 
            
            cmd.delete(chain_obj)

        cmd.delete("target_struct")

    if ensemble_models:
        print(f"\nGenerating ensemble structure: {ENSEMBLE_PDB}")
        for i, model in enumerate(ensemble_models):
            cmd.create("ensemble_all", model, 1, i + 1)
        cmd.save(ENSEMBLE_PDB, "ensemble_all", state=0)
        
        print(f"Generating Multi-FASTA sequence file: {OUTPUT_FASTA}")
        with open(OUTPUT_FASTA, "w") as f_out:
            for seq in fasta_sequences:
                f_out.write(seq.strip() + "\n")
                
        print(f"Successfully processed {len(ensemble_models)} core-fold structures.")
    else:
        print("No matching structural motifs were extracted. Exiting.")
        cmd.quit()
        return

    print(f"\n--- Stage 2: High-Resolution Rendering (600 DPI, TIFF, OpenGL Draw) ---")
    cmd.reinitialize()

    cmd.bg_color("white")
    cmd.set("antialias", 2)            
    
    cmd.set("cartoon_fancy_helices", 1) 
    cmd.set("cartoon_highlight_color", "white")
    cmd.set("light_count", 0) 

    cmd.load(ENSEMBLE_PDB, "ensemble_view")
    cmd.set("all_states", "on") 
    
    cmd.hide("everything")
    cmd.show("cartoon", "ensemble_view")
    cmd.color("wheat", "ensemble_view")

    cmd.orient("ensemble_view")             
    cmd.zoom("ensemble_view", buffer=5.0)   

    temp_png = "temp_render.png"
    render_width = 6000
    render_height = 6000

    print(f"Drawing high-resolution image using OpenGL (width={render_width}, height={render_height})...")
    
    try:
        cmd.draw(render_width, render_height, antialias=2)
        cmd.png(temp_png, dpi=600)
        
        img = Image.open(temp_png)
        img.save(OUTPUT_TIFF, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
        os.remove(temp_png)
        
        print(f"Success! Publishable TIFF image saved to: {OUTPUT_TIFF}")
    except Exception as e:
        print(f"\nDraw rendering failed: {e}")

    cmd.quit()

if __name__ == "__main__":
    main()
