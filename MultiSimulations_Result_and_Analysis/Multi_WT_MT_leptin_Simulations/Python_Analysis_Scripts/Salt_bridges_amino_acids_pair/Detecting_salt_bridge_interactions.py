# imports
import numpy as np
import MDAnalysis
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Pool, Manager
from functools import partial
import argparse
import sys
import os
import json

def validate_file(f):
    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
       	raise argparse.ArgumentTypeError("{0} does not exist".format(f))

    return f

def Screening_salt_bridge(u, ts, InputAA1, InputAA2):

	# AA id extraction
	InputAA1=args.residue1
	ResID1=int(''.join(filter(str.isdigit, InputAA1)))
	ResName1=InputAA1[:3]
	ChainID1=InputAA1[-1]


	InputAA2=args.residue2
	ResID2=int(''.join(filter(str.isdigit, InputAA2)))
	ResName2=InputAA2[:3]
	ChainID2=InputAA2[-1]


	if (ResName1.upper() == "GLU" or ResName1.upper() == "ASP") and (ResName2.upper() == "ARG" or ResName2.upper() == "LYS"):
		u.universe.trajectory[ts]
			
		acidic_o = u.select_atoms(f"(resid {ResID1} and segid {ChainID1}) and (name OD* or name OE*)")
		basic_n = u.select_atoms(f"(resid {ResID2} and segid {ChainID2}) and (name NZ NH*)")
		
		salt_bridges_group=u.select_atoms("group basic and around 3.2 group acidic", acidic=acidic_o, basic=basic_n)
				
		Salt_bridges_dict.update({ts: [f"{ResName1}{ResID1}{ChainID1}_{basic_aa.resname}{basic_aa.resid}{basic_aa.segid}" for basic_aa in salt_bridges_group]})

	elif (ResName2.upper() == "GLU" or ResName2.upper() == "ASP") and (ResName1.upper() == "ARG" or ResName1.upper() == "LYS"):
		u.universe.trajectory[ts]
				
		acidic_o = u.select_atoms(f"(resid {ResID2} and segid {ChainID2}) and (name OD* or name OE*)")
		basic_n = u.select_atoms(f"(resid {ResID1} and segid {ChainID1}) and (name NZ NH*)")
			
		salt_bridges_group=u.select_atoms("group acidic and around 3.2 group basic", acidic=acidic_o, basic=basic_n)
			
		Salt_bridges_dict.update({ts: [f"{acidic_aa.resname}{acidic_aa.resid}{acidic_aa.segid}_{ResName1}{ResID1}{ChainID1}" for acidic_aa in salt_bridges_group]})

if __name__ == '__main__':

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Detecting all potential salt bridges involving selected two amino acids over the course of simulation")
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): pdb", required=True)	
	parser.add_argument("-residue1", help="Specifying the id of amino acid to track for salt bridge. The first amino acid id is based on the id specifid\
						in the flag -tpr (example: TYR140A -> amino acid Tyr at amino acid id of 140 and chain id of A).\
						Note: amino acid triple code must be written in capital", required=True)
	parser.add_argument("-residue2", help="Specifying the id of amino acid to track for salt bridge. The second amino acid id is based on the id specifid\
						in the flag -tpr (example: TYR140A -> amino acid Tyr at amino acid id of 140 and chain id of A).\
						Note: amino acid triple code must be written in capital", required=True)


	args = parser.parse_args()
	

	# Initiate universe
	Environ = MDAnalysis.Universe(args.tpr, args.traj)	
	# Initiate universe
	Salt_bridges_dict=Manager().dict()
	Salt_bridges_count={}	


	# Calculate volume multiprocessing					
	pool = mp.Pool(mp.cpu_count())  
	pool.starmap(Screening_salt_bridge, [(Environ,time_frame,args.residue1,args.residue2) for time_frame in range(len(Environ.universe.trajectory))])
	
	# Re-arrange salt bridges dict
	for key, value in Salt_bridges_dict.items():
		if value:
			contact = value[0]
			if contact not in Salt_bridges_count.keys():
				Salt_bridges_count.update({contact:1})
			else:
				Salt_bridges_count[contact]+=1
	
	# Write to file
	with open(f"{args.residue1}_{args.residue2}_raw.json", "w") as outfile:
		json.dump(Salt_bridges_dict.copy(), outfile)	

	with open(f"{args.residue1}_{args.residue2}_count.json", "w") as outfile:
		json.dump(Salt_bridges_count.copy(), outfile)	























