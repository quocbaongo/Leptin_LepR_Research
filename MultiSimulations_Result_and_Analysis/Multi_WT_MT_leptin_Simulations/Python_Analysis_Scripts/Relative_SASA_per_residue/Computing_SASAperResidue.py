# Import libraries
import os
import argparse
from pymol import cmd
from pymol import stored
import MDAnalysis as mda

def validate_file(f):
        if not os.path.exists(f):
                # Argparse uses the ArgumentTypeError to give a rejection message like:
                # error: argument input: x does not exist
                raise argparse.ArgumentTypeError("{0} does not exist".format(f))
        return f


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Detecting solvent accessible surface area per input protein's amino acids")
	parser.add_argument("--tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)  
	parser.add_argument("--FileOut", help="Specifying the name of out put file (.txt format)", required=True)
	args = parser.parse_args()
	
	TPR = args.tpr
	file_output = args.FileOut
	
	# Extract SASA per residue
	cmd.load(TPR)
	StructureName=os.path.basename(TPR).split(".")[0]
	SASAProfile=cmd.get_sasa_relative(StructureName)
	
	# Write to file
	for key, value in SASAProfile.items():
		file_output_content=open(file_output,"a")
		file_output_content.write(f"{key[3]}{key[2]} {float(value)}\n")
		file_output_content.close()
