import os
import argparse
import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc
from biotite.application.dssp import DsspApp

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Detecting protein secondary structure change during the MD simulation')
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)	
	parser.add_argument("-FileOut", help="Specifying the name of out put file", required=True)
	args = parser.parse_args()

	################################################################
	################## Trajectory and topology #####################
	################################################################	
	
	TRAJ = args.traj
	TPR = args.tpr
	file_out = args.FileOut
	
	xtc_file = xtc.XTCFile.read(TRAJ)
	traj = xtc_file.get_structure(template=strucio.load_structure(TPR))
	time = xtc_file.get_time()
	traj = traj[:, struc.filter_amino_acids(traj)]
	
	# DSSP does not assign an SSE to the last residue -> -1
	sse = np.empty((traj.shape[0], struc.get_residue_count(traj)-1), dtype='U1')

	for idx, frame in enumerate(traj):
		app = DsspApp(traj[idx])
		app.start()
		app.join()
		sse[idx] = app.get_sse()


	# Save the array to a binary file
	np.save(f"{file_out}.npy", sse)

