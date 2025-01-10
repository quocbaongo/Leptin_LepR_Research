# Import libraries
import os
import numpy as np
import argparse
from scipy.spatial import ConvexHull
import MDAnalysis
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Pool
from functools import partial
import sys

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f



def MeasuringVolume(Environ, frame_index, SelectedRes, FileOut):

	
	Environ.universe.trajectory[frame_index]

	CoordGeometry = []
			
	for element in SelectedRes:
		for i in range(element[0], element[1]+1):
			CoordGeometry.append(Environ.select_atoms(f'resid {i} and segid {element[2]} and name CA').positions[0])
			
	CoordGeometry = np.array(CoordGeometry)
	hull = ConvexHull(CoordGeometry)
			
	# Write in file
	Volume = open(FileOut, 'a')
	Volume.write(f'{frame_index} {hull.volume}\n')
	Volume.close()



if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Measuring sampled side chain distances between two specified amino acids over the course of simulation')
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)	
	parser.add_argument("-Region", type=validate_file, help="File in .txt format to define amino acids that belongs to the target region", required=True)
	parser.add_argument("-FileOut", help="Specifying the name of out put file", required=True)
	args = parser.parse_args()

	################################################################
	################## Trajectory and topology #####################
	################################################################	
	
	TRAJ = args.traj
	TPR = args.tpr
	ResIDRegion = args.Region
	file_output = args.FileOut
	
	# Initiate universe
	u = MDAnalysis.Universe(TPR, TRAJ)
	
	# Open input file
	DefinedRegion=[]
	with open(f'{ResIDRegion}', 'r') as f:
		FileContent=f.readlines()
		
		for line in FileContent:
			if line.strip():
				DefinedRegion.append([int(line.split()[0]), int(line.split()[1]), line.split()[2]])


	# Calculate volume multiprocessing					
	pool = mp.Pool(mp.cpu_count())  
	pool.starmap(MeasuringVolume, [(u,ts,DefinedRegion,file_output) for ts in range(len(u.universe.trajectory))])


	
	
	
	
	
