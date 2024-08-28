#!/bin/bash -l


# Coarse-grained molecular dynamics simulations analysis
# Note: in the 3:3 leptin-LR complex, leptin molecules are designated as chain A, C, E
#LR molecules are designated as chain B, D, F


# Required input
WorkDir=/path/to/analysis/directory
WTTrajDir=/path/to/generated/trajectories/WT/Leptin/LR
MTTrajDir=/path/to/generated/trajectories/MT/Lepin/LR
GMX=$GMX


# The simulation of WT/MT leptin-LR were repeated 5 times each
#Thus, we have 5 generated trajectories for each
#1. First step of the anaylsis workflow is to correct for periodicity

# WT leptin simulations

for i in 1 2 3 4 5
do
        cd $WTTrajDir/Rep$i
        mkdir $WTTrajDir/Rep$i/FixPBC
        
        cd $WTTrajDir/Rep$i/FixPBC

        echo "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/MD/traj_comp.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/trjout.whole.xtc -pbc whole 
        
        echo "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/trjout.whole.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -pbc nojump -o $WTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc 
        
        echo "Protein" "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans 

        echo "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -dump 0 -o $WTTrajDir/Rep$i/FixPBC/frame0.pdb 

        rm $WTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc
        rm $WTTrajDir/Rep$i/FixPBC/trjout.whole.xtc

done


# L72S leptin simulations

for i in 1 2 3 4 5
do
        cd $MTTrajDir/Rep$i
        mkdir $MTTrajDir/Rep$i/FixPBC
        
        cd $MTTrajDir/Rep$i/FixPBC

        echo "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/MD/traj_comp.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/trjout.whole.xtc -pbc whole 
        
        echo "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/trjout.whole.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -pbc nojump -o $MTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc 
        
        echo "Protein" "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans 

        echo "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -dump 0 -o $MTTrajDir/Rep$i/FixPBC/frame0.pdb 

        rm $MTTrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc
        rm $MTTrajDir/Rep$i/FixPBC/trjout.whole.xtc

done



#2. Measuring the RMSD evolution of the complex leptin-LR CRH2 domain within the 3:3 leptin-LR complexes

# WT leptin simulations

mkdir $WorkDir/WT
mkdir $WorkDir/WT/RMSDEvolution
cd $WorkDir/WT/RMSDEvolution

# Generate index file containing the backbone (BB) beads of the leptin molecule and LR CRH2 molecule of each monomeric structure(~ from amino acid ILE428 - VAL633)

echo -e 'a BB & chain A \n a BB & chain B & r 428-633 \n "BB_&_chA" | "BB_&_chB_&_r_428-633" \n q' | $GMX make_ndx -f $WTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/WT/RMSDEvolution/indexAB.ndx

echo -e 'a BB & chain C \n a BB & chain D & r 428-633 \n "BB_&_chC" | "BB_&_chD_&_r_428-633" \n q' | $GMX make_ndx -f $WTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/WT/RMSDEvolution/indexCD.ndx

echo -e 'a BB & chain E \n a BB & chain F & r 428-633 \n "BB_&_chE" | "BB_&_chF_&_r_428-633" \n q' | $GMX make_ndx -f $WTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/WT/RMSDEvolution/indexEF.ndx


# Measuring RMSD
for i in 1 2 3 4 5
do
	mkdir $WorkDir/WT/RMSDEvolution/Rep$i
	cd $WorkDir/WT/RMSDEvolution/Rep$i
	
	echo -e "BB_&_chA_BB_&_chB_&_r_428-633" "BB_&_chA_BB_&_chB_&_r_428-633" | $GMX rms -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/WT/RMSDEvolution/Rep$i/fileAB.xvg -n $WorkDir/WT/RMSDEvolution/indexAB.ndx
	
	echo -e "BB_&_chC_BB_&_chD_&_r_428-633" "BB_&_chC_BB_&_chD_&_r_428-633" | $GMX rms -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/WT/RMSDEvolution/Rep$i/fileCD.xvg -n $WorkDir/WT/RMSDEvolution/indexCD.ndx

	echo -e "BB_&_chE_BB_&_chF_&_r_428-633" "BB_&_chE_BB_&_chF_&_r_428-633" | $GMX rms -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/WT/RMSDEvolution/Rep$i/fileEF.xvg -n $WorkDir/WT/RMSDEvolution/indexEF.ndx

done





# L72S leptin simulations

mkdir $WorkDir/L72S
mkdir $WorkDir/L72S/RMSDEvolution
cd $WorkDir/L72S/RMSDEvolution

# Generate index file containing the backbone (BB) beads of the leptin molecule and LR CRH2 molecule of each monomeric structure(~ from amino acid ILE428 - VAL633)

echo -e 'a BB & chain A \n a BB & chain B & r 428-633 \n "BB_&_chA" | "BB_&_chB_&_r_428-633" \n q' | $GMX make_ndx -f $MTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/L72S/RMSDEvolution/indexAB.ndx

echo -e 'a BB & chain C \n a BB & chain D & r 428-633 \n "BB_&_chC" | "BB_&_chD_&_r_428-633" \n q' | $GMX make_ndx -f $MTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/L72S/RMSDEvolution/indexCD.ndx

echo -e 'a BB & chain E \n a BB & chain F & r 428-633 \n "BB_&_chE" | "BB_&_chF_&_r_428-633" \n q' | $GMX make_ndx -f $MTTrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/L72S/RMSDEvolution/indexEF.ndx


# Measuring RMSD
for i in 1 2 3 4 5
do
	mkdir $WorkDir/L72S/RMSDEvolution/Rep$i
	cd $WorkDir/L72S/RMSDEvolution/Rep$i
	
	echo -e "BB_&_chA_BB_&_chB_&_r_428-633" "BB_&_chA_BB_&_chB_&_r_428-633" | $GMX rms -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/L72S/RMSDEvolution/Rep$i/fileAB.xvg -n $WorkDir/L72S/RMSDEvolution/indexAB.ndx
	
	echo -e "BB_&_chC_BB_&_chD_&_r_428-633" "BB_&_chC_BB_&_chD_&_r_428-633" | $GMX rms -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/L72S/RMSDEvolution/Rep$i/fileCD.xvg -n $WorkDir/L72S/RMSDEvolution/indexCD.ndx

	echo -e "BB_&_chE_BB_&_chF_&_r_428-633" "BB_&_chE_BB_&_chF_&_r_428-633" | $GMX rms -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $WorkDir/L72S/RMSDEvolution/Rep$i/fileEF.xvg -n $WorkDir/L72S/RMSDEvolution/indexEF.ndx

done




#3. Tracking the distance between leptin and the LR CRH2 domain over the course of simulation.

# The python script /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py works in such a way that
#it first detects the interface amino acids residues from each molecule using pymol (check https://pymolwiki.org/index.php/InterfaceResidues) in the crystal structure 
#of 3:3 leptin-LR (PDB ID: 8AVF) in coarse-grained level. For each structural frame in the trjectory, 
#the center of geometry of each molecule's group of the mentioned amino acids is then detected, and
#their distance is calculated and saved


# WT leptin simulations

mkdir $WorkDir/WT
mkdir $WorkDir/WT/CenterOfGeometry_distance
cd $WorkDir/WT/CenterOfGeometry_distance



					################################## ATTENTION !!!! ##############################################
					# Modify /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py accordingly


for i in 1 2 3 4 5
do
	mkdir $WorkDir/WT/CenterOfGeometry_distance/Rep$i
	cd $WorkDir/WT/CenterOfGeometry_distance/Rep$i
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 A -ChainID2 B -File_Out $WorkDir/WT/CenterOfGeometry_distance/Rep$i/FileAB.txt
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 C -ChainID2 D -File_Out $WorkDir/WT/CenterOfGeometry_distance/Rep$i/FileCD.txt
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 E -ChainID2 F -File_Out $WorkDir/WT/CenterOfGeometry_distance/Rep$i/FileEF.txt


done





# L72S leptin simulations

mkdir $WorkDir/L72S
mkdir $WorkDir/L72S/CenterOfGeometry_distance
cd $WorkDir/L72S/CenterOfGeometry_distance



					################################## ATTENTION !!!! ##############################################
					# Modify /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py accordingly


for i in 1 2 3 4 5
do
	mkdir $WorkDir/L72S/CenterOfGeometry_distance/Rep$i
	cd $WorkDir/L72S/CenterOfGeometry_distance/Rep$i
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 A -ChainID2 B -File_Out $WorkDir/L72S/CenterOfGeometry_distance/Rep$i/FileAB.txt
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 C -ChainID2 D -File_Out $WorkDir/L72S/CenterOfGeometry_distance/Rep$i/FileCD.txt
	
	python3 /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/VaccumEM/Complex_box.pdb -ChainID1 E -ChainID2 F -File_Out $WorkDir/L72S/CenterOfGeometry_distance/Rep$i/FileEF.txt


done


















