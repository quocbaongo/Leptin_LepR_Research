import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import json
import os
import sys
import MDAnalysis
	
def validate_file(f):

	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError('{0} does not exist'.format(f))
    
	return f	
	
if __name__ == "__main__":
    
	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Plotting distance distribution")
	parser.add_argument("--WTRMSF", type=validate_file, help="RMSF file of WT structure in nm (.xvg)", required=True)
	parser.add_argument("--MTRMSF", type=validate_file, help="RMSF file of MT structure in nm (.xvg)", required=True)
	parser.add_argument("--PDBFile", type=validate_file, help="Input structure file in pdb format", required=True)
	
	parser.add_argument("--FileOut", help="Specifying the name of out put file", required=True)
	args = parser.parse_args()
	
	WTRMSF = args.WTRMSF
	MTRMSF = args.MTRMSF
	PDBFile = args.PDBFile
	file_out = args.FileOut


	# WT protein
	with open(WTRMSF, "r") as f:
		WTRMSF_file=f.readlines() 
    
	WTRMSF_Value=[]	
    
	count=0
	for line in WTRMSF_file:
		if not (line.startswith("#") or line.startswith("@")):
			count+=1
			WTRMSF_Value.append([float(line.split()[1]),count])

	# MT protein
	with open(MTRMSF, "r") as f:
		MTRMSF_file=f.readlines() 
    
	MTRMSF_Value=[]	
    
	count=0
	for line in MTRMSF_file:
		if not (line.startswith("#") or line.startswith("@")):
			count+=1
			MTRMSF_Value.append([float(line.split()[1]),count])

	# PDB Structure
	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(PDBFile)

	# Select resid and resname
	PDB_resid_resname={}

	for atom in uPDB.atoms:
		if atom.resid not in PDB_resid_resname.keys():
			PDB_resid_resname.update({atom.resid: atom.resname.capitalize()})

	# Plotting
	fig, ax = plt.subplots(1, 1, figsize=(12, 7))

	# Manually generate xticks
	xticks=[]
	for i in range(list(PDB_resid_resname.keys())[0], list(PDB_resid_resname.keys())[-1] +1):
		if i == list(PDB_resid_resname.keys())[0]:
			xticks.append(str(i))
		elif i % 20 == 0:
			xticks.append(str(i))
		elif i == list(PDB_resid_resname.keys())[-1]:
			xticks.append(str(i))
		else:
			xticks.append("")

	
	# Leptin RMSF
	ax.plot([i[1] for i in WTRMSF_Value], [i[0] for i in WTRMSF_Value],"k-",color="cyan",label="WT")
	ax.plot([i[1] for i in MTRMSF_Value], [i[0] for i in MTRMSF_Value], "k-",color="red",label="MT")	

	ax.legend(loc="upper right")
	ax.set_xticks([i[1] for i in WTRMSF_Value],xticks)
	ax.set_xticklabels(xticks, rotation=30, ha="right",fontsize=7)
	ax.set_xlabel("Residue ID",labelpad=30,size=13)
	ax.set_ylabel("RMSF (nm)",labelpad=10,size=13)
	ax.legend(prop={"size": 8})

	fig.tight_layout()
	plt.savefig("RMSF.png",dpi=600,bbox_inches="tight")




	




