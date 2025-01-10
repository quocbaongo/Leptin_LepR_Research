import os
import json
import sys
import argparse
import MDAnalysis
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

if __name__ == "__main__":

	# Define
	InputFile=sys.argv[1]
	PDBFile=sys.argv[2]
	ResID=sys.argv[3]
	
	TotalFrame=100001

	# Atom information in universe u using .pdb file as topology
	# Load topology and trajectory to MDAnalysis
	uPDB = MDAnalysis.Universe(PDBFile)

	# Select resid and resname
	PDB_resid_resname={}

	for atom in uPDB.atoms:
		if atom.resid not in PDB_resid_resname.keys():
			PDB_resid_resname.update({f"{atom.resid}{atom.segid}": f"{atom.resname.capitalize()}{atom.resid}{atom.segid}"})

	HelicesResidues=[f"{i}A" for i in range(23,48)]
	HelicesResidues+=[f"{i}A" for i in range(72,89)]
	HelicesResidues+=[f"{i}A" for i in range(92,116)]
	HelicesResidues+=[f"{i}A" for i in range(141,165)]
	HelicesResidues=[PDB_resid_resname[i] for i in HelicesResidues]
	
	ABLoopResidues=[f"{i}A" for i in range(48,72)]
	ABLoopResidues=[PDB_resid_resname[i] for i in ABLoopResidues]

	CDLoopResidues=[f"{i}A" for i in range(116,141)]	
	CDLoopResidues=[PDB_resid_resname[i] for i in CDLoopResidues]


	# Dataframe to add
	df = pd.DataFrame(columns = ["Interaction amino acid pairs","Interaction atom pairs","# frames","Probability"])	

	# Open input file
	with open(InputFile) as f:
		HydrogenFileContent = f.read().splitlines()    


	# Read each hydrogen bond and filter out intramolecular hbonding
	HydrogenBondConnection={}
	HydrogenBondCountPerTime={}
	Count_BondsToHelices=0
	Count_BondsToABLoop=0
	Count_BondsToCDLoop=0
	
	for i in range(100001):
		HydrogenBondCountPerTime.update({i:0})

	for line in HydrogenFileContent:
		Time=int(line.split()[0])
		Temp1_ResID=f"{line.split()[1].split('_')[1].capitalize()}{line.split()[1].split('_')[0]}"
		Temp2_ResID=f"{line.split()[3].split('_')[1].capitalize()}{line.split()[3].split('_')[0]}"
		
		if Temp1_ResID[3:] == ResID:
			if Temp2_ResID in HelicesResidues:
				Count_BondsToHelices+=1
			elif Temp2_ResID in ABLoopResidues:
				Count_BondsToABLoop+=1
			elif Temp2_ResID in CDLoopResidues:
				Count_BondsToCDLoop+=1
		elif Temp2_ResID[3:] == ResID:
			if Temp1_ResID in HelicesResidues:
				Count_BondsToHelices+=1
			elif Temp1_ResID in ABLoopResidues:
				Count_BondsToABLoop+=1
			elif Temp1_ResID in CDLoopResidues:
				Count_BondsToCDLoop+=1
		
		HbondPattern=f"{line.split()[1].split('_')[2]} {Temp1_ResID}······ {line.split()[3].split('_')[2]} {Temp2_ResID}"
				
		HydrogenBondCountPerTime[Time]+=1
		
		key=f"{Temp1_ResID}-{Temp2_ResID}"
		
		if key not in HydrogenBondConnection.keys():
			HydrogenBondConnection.update({key:(HbondPattern, 0)})
		else:
			count=HydrogenBondConnection[key][1]
			count+=1
			HydrogenBondConnection[key]=(HydrogenBondConnection[key][0],count)

	# Number of hydrogen bond count per time
	HydrogenBondCountPerTimeKeys = list(HydrogenBondCountPerTime.keys())
	HydrogenBondCountPerTimeKeys.sort()
	
	# Sorted dictionary
	HydrogenBondCountPerTimeSorted = {i: HydrogenBondCountPerTime[i] for i in HydrogenBondCountPerTimeKeys}
	
	for key, value in HydrogenBondCountPerTimeSorted.items():
		f = open("HbondPerTime.txt", "a")
		f.write(f"{key} {value}\n")
		f.close()

	# Normalize
	for key in HydrogenBondConnection.keys():
		HydrogenBondConnection[key]=(HydrogenBondConnection[key][0], HydrogenBondConnection[key][1], HydrogenBondConnection[key][1] / TotalFrame)
		
	for key, value in HydrogenBondConnection.items():
		row_to_append = pd.DataFrame([{"Interaction amino acid pairs": key,
						"Interaction atom pairs": value[0],
						"# frames": value[1],
						"Probability": value[2]}])

		df = pd.concat([df, row_to_append])
	
	df.reset_index(drop=True, inplace=True)
	sorted_df = df.sort_values(by=['# frames'], ascending=False)	
	
	# Write to file
	with pd.ExcelWriter(f"{ResID}.xlsx") as writer:
		sorted_df.to_excel(writer)
		
	print(f"Total Hbonds to Helices: {Count_BondsToHelices}")
	print(f"Total Hbonds to AB loop: {Count_BondsToABLoop}")
	print(f"Total Hbonds to CD loop: {Count_BondsToCDLoop}")














