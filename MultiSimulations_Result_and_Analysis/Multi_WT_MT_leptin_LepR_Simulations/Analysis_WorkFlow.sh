#!/bin/bash -l

# Atomistic molecular dynamics simulations analysis
# Required input
WorkDir=/path/to/analysis/directory
AnalysisScripts=/path/to/generated/Python_Analysis_Scripts
TrajDir=/path/to/generated/trajectories
GMX=gmx_2021.2

# Coarse-grained molecular dynamics simulations analysis
# Note: in the 3:3 leptin-LR complex, leptin molecules are designated as chain A, C, E
#LR molecules are designated as chain B, D, F

# The simulation of WT/MT leptin-LepR were repeated 10 times
#Thus, we have 10 generated trajectories for each
#1. First step of the anaylsis workflow is to correct for periodicity

for i in 1 2 3 4 5 6 7 8 9 10
do
        cd $TrajDir/Rep$i
        mkdir $TrajDir/Rep$i/FixPBC

        cd $TrajDir/Rep$i/FixPBC

        echo "Protein" | $GMX trjconv -f $TrajDir/Rep$i/MD/traj_comp.xtc -s $TrajDir/Rep$i/EM/em.tpr -o $TrajDir/Rep$i/FixPBC/trjout.whole.xtc -pbc whole

        echo "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/trjout.whole.xtc -s $TrajDir/Rep$i/EM/em.tpr -pbc nojump -o $TrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc

        echo "Protein" "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc -s $TrajDir/Rep$i/EM/em.tpr -o $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans

        echo "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/Rep$i/EM/em.tpr -dump 0 -o $TrajDir/Rep$i/FixPBC/frame0.pdb

        rm $TrajDir/Rep$i/FixPBC/trjout.whole.nojump.xtc
        rm $TrajDir/Rep$i/FixPBC/trjout.whole.xtc

done


# Concatenate trjectories
mkdir $WorkDir/ConcatTraj
cd $WorkDir/ConcatTraj

$GMX trjcat -f $TrajDir/Rep*/FixPBC/md_nojump_whole_fit.xtc -cat -o $WorkDir/ConcatTraj/trajout.xtc


#2. Measuring the RMSD evolution of the complex leptin-LR CRH2 domain and leptin-LR IgD domain within the 3:3 leptin-LR complexes

                                ############################################# Leptin-LR CRH2 #######################################################


mkdir $WorkDir/RMSDEvolution_leptin_CRH2
cd $WorkDir/RMSDEvolution_leptin_CRH2

# Generate index file containing the CA atoms of the leptin molecule and LR CRH2 molecule of each monomeric structure(~ from amino acid ILE428 - VAL633)

echo -e 'a CA & chain A \n a CA & chain B & r 428-633 \n "CA_&_chA" | "CA_&_chB_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexAB.ndx

echo -e 'a CA & chain C \n a CA & chain D & r 428-633 \n "CA_&_chC" | "CA_&_chD_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexCD.ndx

echo -e 'a CA & chain E \n a CA & chain F & r 428-633 \n "CA_&_chE" | "CA_&_chF_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexEF.ndx


echo -e '"Protein" & chain A \n "Protein" & chain B & r 428-633 \n "Protein_&_chA" | "Protein_&_chB_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexAB_all_atoms.ndx

echo -e '"Protein" & chain C \n "Protein" & chain D & r 428-633 \n "Protein_&_chC" | "Protein_&_chD_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexCD_all_atoms.ndx

echo -e '"Protein" & chain E \n "Protein" & chain F & r 428-633 \n "Protein_&_chE" | "Protein_&_chF_&_r_428-633" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/indexEF_all_atoms.ndx


# Measuring RMSD
for i in 1 2 3 4 5 6 7 8 9 10
do
        mkdir $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i
        cd $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i

        # CA atoms

        echo -e "CA_&_chA_CA_&_chB_&_r_428-633" "CA_&_chA_CA_&_chB_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileAB.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexAB.ndx

        echo -e "CA_&_chC_CA_&_chD_&_r_428-633" "CA_&_chC_CA_&_chD_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileCD.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexCD.ndx

        echo -e "CA_&_chE_CA_&_chF_&_r_428-633" "CA_&_chE_CA_&_chF_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileEF.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexEF.ndx


        # All atoms
        echo -e "Protein_&_chA_Protein_&_chB_&_r_428-633" "Protein_&_chA_Protein_&_chB_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileAB_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexAB_all_atoms.ndx

        echo -e "Protein_&_chC_Protein_&_chD_&_r_428-633" "Protein_&_chC_Protein_&_chD_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileCD_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexCD_all_atoms.ndx

        echo -e "Protein_&_chE_Protein_&_chF_&_r_428-633" "Protein_&_chE_Protein_&_chF_&_r_428-633" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_CRH2/Rep$i/fileEF_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_CRH2/indexEF_all_atoms.ndx

done

# Concat files
cat $WorkDir/RMSDEvolution_leptin_CRH2/Rep*/fileAB.xvg $WorkDir/RMSDEvolution_leptin_CRH2/Rep*/fileCD.xvg $WorkDir/RMSDEvolution_leptin_CRH2/Rep*/fileEF.xvg > $WorkDir/RMSDEvolution_leptin_CRH2/RMSDLeptin_CRH2.xvg


                                ############################################# Leptin-LR IgD #######################################################

mkdir $WorkDir/RMSDEvolution_leptin_IgD
cd $WorkDir/RMSDEvolution_leptin_IgD

# Generate index file containing the CA atoms of the leptin molecule and LR IgD molecule of each monomeric structure (~ from amino acid ILE636 - THR829)

echo -e 'a CA & chain A \n a CA & chain F & r 333-426 \n "CA_&_chA" | "CA_&_chF_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexAF.ndx

echo -e 'a CA & chain C \n a CA & chain B & r 333-426 \n "CA_&_chC" | "CA_&_chB_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexCB.ndx

echo -e 'a CA & chain E \n a CA & chain D & r 333-426 \n "CA_&_chE" | "CA_&_chD_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexED.ndx


echo -e '"Protein" & chain A \n "Protein" & chain F & r 333-426 \n "Protein_&_chA" | "Protein_&_chF_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexAF_all_atoms.ndx

echo -e '"Protein" & chain C \n "Protein" & chain B & r 333-426 \n "Protein_&_chC" | "Protein_&_chB_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexCB_all_atoms.ndx

echo -e '"Protein" & chain E \n "Protein" & chain D & r 333-426 \n "Protein_&_chE" | "Protein_&_chD_&_r_333-426" \n q' | $GMX make_ndx -f $TrajDir/Rep1/FixPBC/frame0.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/indexED_all_atoms.ndx


# Measuring RMSD
for i in 1 2 3 4 5 6 7 8 9 10
do
        mkdir $WorkDir/RMSDEvolution_leptin_IgD/Rep$i
        cd $WorkDir/RMSDEvolution_leptin_IgD/Rep$i

        # CA atoms
        echo -e "CA_&_chA_CA_&_chF_&_r_333-426" "CA_&_chA_CA_&_chF_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileAF.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexAF.ndx

        echo -e "CA_&_chC_CA_&_chB_&_r_333-426" "CA_&_chC_CA_&_chB_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileCB.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexCB.ndx

        echo -e "CA_&_chE_CA_&_chD_&_r_333-426" "CA_&_chE_CA_&_chD_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileED.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexED.ndx


        # All atoms
        echo -e "Protein_&_chA_Protein_&_chF_&_r_333-426" "Protein_&_chA_Protein_&_chF_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileAF_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexAF_all_atoms.ndx

        echo -e "Protein_&_chC_Protein_&_chB_&_r_333-426" "Protein_&_chC_Protein_&_chB_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileCB_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexCB_all_atoms.ndx

        echo -e "Protein_&_chE_Protein_&_chD_&_r_333-426" "Protein_&_chE_Protein_&_chD_&_r_333-426" | $GMX rms -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/VaccumEM/box.pdb -o $WorkDir/RMSDEvolution_leptin_IgD/Rep$i/fileED_all_atoms.xvg -n $WorkDir/RMSDEvolution_leptin_IgD/indexED_all_atoms.ndx

done

# Concat files
cat $WorkDir/RMSDEvolution_leptin_IgD/Rep*/fileAF.xvg $WorkDir/RMSDEvolution_leptin_IgD/Rep*/fileCB.xvg $WorkDir/RMSDEvolution_leptin_IgD/Rep*/fileED.xvg > $WorkDir/RMSDEvolution_leptin_IgD/RMSDLeptin_IgD.xvg


#3. Tracking the distance between leptin and the LR CRH2 domain over the course of simulation.

# The python script /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py works in such a way that
#it first detects the interface amino acids residues from each molecule using pymol (check https://pymolwiki.org/index.php/InterfaceResidues) using
#the complex structure provided through the flag '-tpr'. We used the crystal structure of 3:3 leptin-LR (PDB ID: 8AVF) for initial interface residues detection.
#Afterwards, for each structural frame in the trjectory, the center of geometry of each molecule's group of the mentioned amino acids is then detected, and
#their distance is calculated and saved in the text file with name defined through the flag '-File_Out'

mkdir $WorkDir/CenterOfGeometry_distance
cd $WorkDir/CenterOfGeometry_distance



                                        ################################## ATTENTION !!!! ##############################################
                                        # Modify /path/to/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py accordingly


for i in 1 2 3 4 5 6 7 8 9 10
do
        mkdir $WorkDir/CenterOfGeometry_distance/Rep$i
        cd $WorkDir/CenterOfGeometry_distance/Rep$i

        # leptin - LR CRH2
        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 A -ChainID2 B -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileAB.txt

        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 C -ChainID2 D -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileCD.txt

        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 E -ChainID2 F -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileEF.txt


        # leptin - LR IgD
        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 A -ChainID2 F -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileAF.txt

        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 C -ChainID2 B -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileCB.txt

        python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 E -ChainID2 D -File_Out $WorkDir/CenterOfGeometry_distance/Rep$i/FileED.txt

done



# Measuring initial distance

# leptin - LR CRH2

python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 A -ChainID2 B -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileAB.txt

python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 C -ChainID2 D -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileCD.txt

python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 E -ChainID2 F -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileEF.txt

# leptin - LR IgD
python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 A -ChainID2 F -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileAF.txt

python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 C -ChainID2 B -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileCB.txt

python3 $AnalysisScripts/CenterOfMass_InterfaceResidues_distance/Measuring_COM_distance.py -traj $TrajDir/VaccumEM/box.pdb -tpr $TrajDir/VaccumEM/box.pdb -ChainID1 E -ChainID2 D -File_Out $WorkDir/CenterOfGeometry_distance/InitialFileED.txt


# Concat files
cat $WorkDir/CenterOfGeometry_distance/Rep*/FileAB.txt $WorkDir/CenterOfGeometry_distance/Rep*/FileCD.txt $WorkDir/CenterOfGeometry_distance/Rep*/FileEF.txt > $WorkDir/CenterOfGeometry_distance/DistanceLeptin_CRH2.txt

cat $WorkDir/CenterOfGeometry_distance/Rep*/FileAF.txt $WorkDir/CenterOfGeometry_distance/Rep*/FileCB.txt $WorkDir/CenterOfGeometry_distance/Rep*/FileED.txt > $WorkDir/CenterOfGeometry_distance/DistanceLeptin_IgD.txt


cat $WorkDir/CenterOfGeometry_distance/InitialFileAB.txt $WorkDir/CenterOfGeometry_distance/InitialFileCD.txt $WorkDir/CenterOfGeometry_distance/InitialFileEF.txt > $WorkDir/CenterOfGeometry_distance/InitialFileLeptin_CRH2_dist.txt

cat $WorkDir/CenterOfGeometry_distance/InitialFileAF.txt $WorkDir/CenterOfGeometry_distance/InitialFileCB.txt $WorkDir/CenterOfGeometry_distance/InitialFileED.txt > $WorkDir/CenterOfGeometry_distance/InitialFileLeptin_IgD_dist.txt













