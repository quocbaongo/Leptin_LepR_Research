import argparse
import numpy as np
from pymol import cmd
from pymol import stored
import os
from itertools import groupby
import sys
import MDAnalysis
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Pool
from functools import partial

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f
	
def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName='interface'):
	'''
	interfaceResidues -- finds 'interface' residues between two chains in a complex.
 
	PARAMS
		cmpx
			The complex containing cA and cB
 
		cA
			The first chain in which we search for residues at an interface
			with cB
 
		cB
			The second chain in which we search for residues at an interface
			with cA
 
		cutoff
			The difference in area OVER which residues are considered
			interface residues.  Residues whose dASA from the complex to
			a single chain is greater than this cutoff are kept.  Zero
			keeps all residues.
 
		selName
			The name of the selection to return.
 
	RETURNS
		* A selection of interface residues is created and named
			depending on what you passed into selName
		* An array of values is returned where each value is:
			( modelName, residueNumber, dASA )
 
	NOTES
		If you have two chains that are not from the same PDB that you want
		to complex together, use the create command like:
			create myComplex, pdb1WithChainA or pdb2withChainX
		then pass myComplex to this script like:
			interfaceResidues myComlpex, c. A, c. X
 
		This script calculates the area of the complex as a whole.  Then,
		it separates the two chains that you pass in through the arguments
		cA and cB, alone.  Once it has this, it calculates the difference
		and any residues ABOVE the cutoff are called interface residues.
 
	AUTHOR:
		Jason Vertrees, 2009.		
	'''
	# Save user's settings, before setting dot_solvent
	oldDS = cmd.get('dot_solvent')
	cmd.set('dot_solvent', 1)
 
	# set some string names for temporary objects/selections
	tempC, selName1 = 'tempComplex', selName+'1'
	chA, chB = 'chA', 'chB'
 
	# operate on a new object & turn off the original
	cmd.create(tempC, cmpx)
	cmd.disable(cmpx)
 
	# remove cruft and inrrelevant chains
	cmd.remove(tempC + ' and not (polymer and (%s or %s))' % (cA, cB))
 
	# get the area of the complete complex
	cmd.get_area(tempC, load_b=1)
	# copy the areas from the loaded b to the q, field.
	cmd.alter(tempC, 'q=b')
 
	# extract the two chains and calc. the new area
	# note: the q fields are copied to the new objects
	# chA and chB
	cmd.extract(chA, tempC + ' and (' + cA + ')')
	cmd.extract(chB, tempC + ' and (' + cB + ')')
	cmd.get_area(chA, load_b=1)
	cmd.get_area(chB, load_b=1)
 
	# update the chain-only objects w/the difference
	cmd.alter( '%s or %s' % (chA,chB), 'b=b-q' )
 
	# The calculations are done.  Now, all we need to
	# do is to determine which residues are over the cutoff
	# and save them.
	stored.r, rVal, seen = [], [], []
	cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')
 
	cmd.enable(cmpx)
	cmd.select(selName1, None)
	for (model,resi,diff) in stored.r:
		key=resi+'-'+model
		if abs(diff)>=float(cutoff):
			if key in seen: continue
			else: seen.append(key)
			rVal.append( (model,resi,diff) )
			# expand the selection here; I chose to iterate over stored.r instead of
			# creating one large selection b/c if there are too many residues PyMOL
			# might crash on a very large selection.  This is pretty much guaranteed
			# not to kill PyMOL; but, it might take a little longer to run.
			cmd.select( selName1, selName1 + ' or (%s and i. %s)' % (model,resi))
 
	# this is how you transfer a selection to another object.
	cmd.select(selName, cmpx + ' in ' + selName1)
	# clean up after ourselves
	cmd.delete(selName1)
	cmd.delete(chA)
	cmd.delete(chB)
	cmd.delete(tempC)
	# show the selection
	cmd.enable(selName)
 
	# reset users settings
	cmd.set('dot_solvent', oldDS)
 
	return rVal

def DetectingInterfaceResidues(HA_chains, Abs_chains, MinimizedPDBFile):

	
	'''
	Function to detect interface residues 

	Arguments:
	
	HA_chains          -- Hemagglutinin chains 
	Abs_chains         -- Antibody chains
	MinimizedPDBFile   -- Energy minimized PDB structure
	
	Returns:
	
	List containing interface residue belonging to hemagglutinin and antibody	

	'''
	
	HAInterfaceRes = []
	cmd.load(MinimizedPDBFile)

	HA_chainsList = list(HA_chains)
	Abs_chainsList = list(Abs_chains)

	for i, j in zip(HA_chainsList, Abs_chainsList):
		if i == ' ':
			HA_chainsList.remove(i)
		elif j == ' ':
			Abs_chainsList.remove(j)


	for chain1 in HA_chainsList:
		for chain2 in Abs_chainsList:
			interface = interfaceResidues(os.path.basename(MinimizedPDBFile).split('.')[0], 
						      cA='c. ' + chain1, 
						      cB='c. ' + chain2, 
						      cutoff=1.0,		
						      selName='interface')
			
			for residue in interface:
				stored.list = []
				if residue[0] == 'chA':
					cmd.iterate('resi ' + str(residue[1]) + ' and chain ' + chain1, 
						'stored.list.append(resn)')
					HAInterfaceRes.append([residue[1], np.unique(stored.list)[0], 
							       chain1])
	
	output = [list(v) for i, v in groupby(HAInterfaceRes, lambda x: x[2])]
	
	# Delete repetition
	for i in range(len(output)):
		output[i] = [list(j) for j in set(map(tuple, output[i]))]	
			
	return output
	
def MeasuringDistance(Ref, Target):
	
	x1=Ref[0]
	y1=Ref[1]
	z1=Ref[2]
	
	x2=Target[0]
	y2=Target[1]
	z2=Target[2]
	
	dist=((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
	
	return dist

def MeasuringDistance_function(u,ts,ListRes1,ListRes2,FileOut):

	u.universe.trajectory[ts]
	# First selection
	Selection1=''
	for res in ListRes1:
	    Selection1+=f'(resid {int(res[0])} and segid {res[2]} and name CA) or '
	    
	Selection1=Selection1[:-4]
	Sele1_COG=u.select_atoms(Selection1).center_of_geometry()
	
	
	# Second selection
	Selection2=''
	for res in ListRes2:
	    Selection2+=f'(resid {int(res[0])} and segid {res[2]} and name CA) or '
	
	Selection2=Selection2[:-4]
	Sele2_COG=u.select_atoms(Selection2).center_of_geometry()
	
	# Measuring distance
	Dist=MeasuringDistance(Sele1_COG, Sele2_COG)
	
	# Write to file
	File=open(FileOut,"a")
	File.write(f"{ts} {float(Dist)/10.0}\n")
	File.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Detecting the interface residues between 2 molecules (specified by chainID) in the initial structure provided by -tpr flag. For each frame of the trajectory, the center of mass of each group of residues belonging to a specific molecule is detected, and their distance is measured.')
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)
	parser.add_argument("-ChainID1", help="Specifying chain ID of first molecule", required=True)
	parser.add_argument("-ChainID2", help="Specifying chain ID of second molecule", required=True)
	parser.add_argument("-File_Out", help="Name of output file, and measured distance displayed in nm", required=True)
	args = parser.parse_args()

	TRAJ = args.traj
	TPR = args.tpr
	FirstChainID = args.ChainID1
	SecondChainID = args.ChainID2
	file_out = args.File_Out
	
	FirstInterfaceRes=DetectingInterfaceResidues(FirstChainID, SecondChainID, TPR)
	SecondInterfaceRes=DetectingInterfaceResidues(SecondChainID, FirstChainID, TPR)
	
	# Initiate universe
	Environ = MDAnalysis.Universe(TPR, TRAJ)
	
	# Start measuring
	pool = mp.Pool(mp.cpu_count())
	pool.starmap(MeasuringDistance_function, [(Environ,frameID,FirstInterfaceRes[0],SecondInterfaceRes[0],file_out) for frameID in range(len(Environ.universe.trajectory))])








    
    
    
    
    
    
    
    	
	

