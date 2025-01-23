# imports
import os
import MDAnalysis
import argparse
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Pool

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f

def MeasuringDistanceFunc1(Environ, TargetResID, RefResID, frame_index):
    
	Environ.universe.trajectory[frame_index]
	
	# Center of geometry of target residue side chain
	if Environ.select_atoms(f"resid {TargetResID}").resnames[0] == 'GLY':
		TargetCoordGeometry = Environ.select_atoms(f"resid {TargetResID} and name CA").center_of_geometry()
		x1=TargetCoordGeometry[0]
		y1=TargetCoordGeometry[1]
		z1=TargetCoordGeometry[2]        
	else:
		TargetCoordGeometry = Environ.select_atoms(f"resid {TargetResID} and not backbone and not type H").center_of_geometry()
		x1=TargetCoordGeometry[0]
		y1=TargetCoordGeometry[1]
		z1=TargetCoordGeometry[2]    

	# Center of geometry of reference residue side chain
	if Environ.select_atoms(f"resid {RefResID}").resnames[0] == 'GLY':
		RefCoordGeometry = Environ.select_atoms(f"resid {RefResID} and name CA").center_of_geometry()
		x2=RefCoordGeometry[0]
		y2=RefCoordGeometry[1]
		z2=RefCoordGeometry[2]
	else:
		RefCoordGeometry = Environ.select_atoms(f"resid {RefResID} and not backbone and not type H").center_of_geometry()
		x2=RefCoordGeometry[0]
		y2=RefCoordGeometry[1]
		z2=RefCoordGeometry[2]    

	# Measuring distance
	dist=((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
	
	# Write in file
	Distance = open(f"file_{TargetResID}_{RefResID}.txt", 'a')
	Distance.write(f"{frame_index} {dist}\n")
	Distance.close()

def MeasuringDistanceFunc2(Environ, TargetResID, RefResID, TargetResChainID, RefResChainID, frame_index):
    
	Environ.universe.trajectory[frame_index]
	
	# Center of geometry of target residue side chain
	if Environ.select_atoms(f"resid {TargetResID}").resnames[0] == 'GLY':
		TargetCoordGeometry = Environ.select_atoms(f"resid {TargetResID} and segid {TargetResChainID} and name CA").center_of_geometry()
		x1=TargetCoordGeometry[0]
		y1=TargetCoordGeometry[1]
		z1=TargetCoordGeometry[2]        
	else:
		TargetCoordGeometry = Environ.select_atoms(f"(resid {TargetResID} and segid {TargetResChainID}) and (not backbone and not type H)").center_of_geometry()
		x1=TargetCoordGeometry[0]
		y1=TargetCoordGeometry[1]
		z1=TargetCoordGeometry[2]    

	# Center of geometry of reference residue side chain
	if Environ.select_atoms(f"resid {RefResID}").resnames[0] == 'GLY':
		RefCoordGeometry = Environ.select_atoms(f"resid {RefResID} and segid {RefResChainID} and name CA").center_of_geometry()
		x2=RefCoordGeometry[0]
		y2=RefCoordGeometry[1]
		z2=RefCoordGeometry[2]
	else:
		RefCoordGeometry = Environ.select_atoms(f"(resid {RefResID} and segid {RefResChainID}) and (not backbone and not type H)").center_of_geometry()
		x2=RefCoordGeometry[0]
		y2=RefCoordGeometry[1]
		z2=RefCoordGeometry[2]    

	# Measuring distance
	dist=((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
	
	# Write in file
	Distance = open(f"file_{TargetResID}_{RefResID}.txt", 'a')
	Distance.write(f"{frame_index} {dist}\n")
	Distance.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Measuring sampled side chain distances between two specified amino acids over the course of simulation')
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)
	parser.add_argument("-resi1", help="Specifying selected first residue ID", required=True)
	parser.add_argument("-resi2", help="Specifying selected second residue ID", required=True)
	parser.add_argument("-ChainID1", help="Specifying selected first residue's chainID", required=False)
	parser.add_argument("-ChainID2", help="Specifying selected second residue's chainID", required=False)
	args = parser.parse_args()


	################################################################
	################## Trajectory and topology #####################
	################################################################	
	
	TRAJ = args.traj
	TPR = args.tpr
	TargetResidueID = args.resi1
	TargetResidueChainID = args.ChainID1
	ReferenceResidueID = args.resi2
	ReferenceResidueChainID = args.ChainID2
	
	u = MDAnalysis.Universe(TPR, TRAJ)
	
	pool = mp.Pool(mp.cpu_count())
		
	# Start collecting sampled distance
	if TargetResidueChainID and ReferenceResidueChainID:
		pool.starmap(MeasuringDistanceFunc2, [(u,str(TargetResidueID),str(ReferenceResidueID),str(TargetResidueChainID),str(ReferenceResidueChainID),ts) for ts in range(len(u.universe.trajectory))])
		
	else:
		pool.starmap(MeasuringDistanceFunc1, [(u,str(TargetResidueID),str(ReferenceResidueID),ts) for ts in range(len(u.universe.trajectory))])



































