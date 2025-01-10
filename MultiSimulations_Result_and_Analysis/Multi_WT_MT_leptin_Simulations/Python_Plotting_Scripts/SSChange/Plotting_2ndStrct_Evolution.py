import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors
import matplotlib as mpl
import MDAnalysis
import gromacs.formats

def validate_file(f):

	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    
	return f

def sse_to_num(sse):

	num = np.empty(sse.shape, dtype=int)
	num[sse == "Coil"] = 0
	num[sse == "B-Sheet"] = 1
	num[sse == "B-Bridge"] = 2
	num[sse == "Bend"] = 3
	num[sse == "Turn"] = 4
	num[sse == "A-Helix"] = 5
	num[sse == "5-Helix"] = 6
	num[sse == "3-Helix"] = 7

	return num


if __name__ == "__main__":

	# flag list
	parser = argparse.ArgumentParser()
    
	# Add an argument for the list of strings
	parser.add_argument("--InputFile", type=validate_file, help="Input file for protein's secondary structures over an MD simulation in .xpm format")
	parser.add_argument("--InputStructure", type=validate_file, help="Input simulated structure in .pdb format")
	parser.add_argument("--File_Out", help="Output image file (in .png form)")

	# Parse the command-line arguments
	args = parser.parse_args()
    
	# Reading secondary structure file
	xpm=gromacs.formats.XPM(f"{args.InputFile}")
	Proteinsse = sse_to_num(xpm.array)
	
	# Start Plotting
	# SSE colormap
	color_assign = {r"coil": "white",
			r"$\beta$-sheet": "red",
			r"$\beta$-bridge": "black",
			r"bend": "green",
			r"turn": "yellow",
			r"$\alpha$-helix": "blue",
			r"$3_{10}$-helix": "gray",
			r"$\pi$-helix": "purple",
			}

	cmap = colors.ListedColormap(color_assign.values())


	# Initialise the subplot function using number of rows and columns 
	fig, ax = plt.subplots(1, 1, figsize=(14,5)) 
	# WT
	ax.imshow(Proteinsse.T, cmap=cmap, origin='lower', aspect='auto')
	ax.set_xlabel("# Structural frames",labelpad=6, size=10)
	ax.set_ylabel("Residue ID",labelpad=6, size=10)

	# x axis tick label
	xticks = np.arange(1, len(Proteinsse))
	xticks = [str(i) for i in xticks]

	for i in range(len(xticks)):
		if int(xticks[i]) == 1:
			xticks[i] = str(1)
		elif not (int(xticks[i]) % 10000 == 0):
			xticks[i] = ""

	ax.set_xticks([i for i in range(len(xticks))],xticks)
	ax.set_xticklabels(xticks, rotation = 30)

	# y axis tick
	yticks = []
	u = MDAnalysis.Universe(f"{args.InputStructure}", f"{args.InputStructure}")

	for residue in u.select_atoms("all").residues:
		yticks.append(residue.resid)

	yticks = [str(i) for i in yticks]

	for i in range(len(yticks)):
		if (int(yticks[i]) % 10 == 0) or (int(yticks[i]) == 22) or (int(yticks[i]) == 166):
			yticks[i] = str(yticks[i])
		else:
			yticks[i] = ""
        
	ax.set_yticks([i for i in range(len(yticks))],yticks)

	# Custom legend below the DSSP plot
	custom_lines = [Line2D([0], [0], color=cmap(i), lw=4) for i in range(len(color_assign))]

	ax.legend(custom_lines, color_assign.keys(), ncol=len(color_assign), fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.3))

	# Adjust the layout to make room for the legend
	plt.tight_layout()
	# Save the plot
	plt.savefig(f"{args.File_Out}.png", dpi=600)

