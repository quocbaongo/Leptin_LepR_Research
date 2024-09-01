#!/bin/bash -l


# Atomistic molecular dynamics simulations analysis


# Required input
WorkDir=/path/to/analysis/directory
WTTrajDir=/path/to/generated/trajectories/WT
MTTrajDir=/path/to/generated/trajectories/MT
GMX=$GMX


# The simulation of WT/MT leptin were repeated 10 times
#Thus, we have 10 generated trajectories for each
#1. First step of the anaylsis workflow is to correct for periodicity

# WT leptin simulations
for i in 1 2 3 4 5 6 7 8 9 10
do 
	cd $WTTrajDir/Rep$i
	mkdir $WTTrajDir/Rep$i/FixPBC
	
	cd $WTTrajDir/Rep$i/FixPBC

	echo "Protein" "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/MD/traj_comp.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/md_nojump.xtc -center -pbc nojump 
	
	echo "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/md_nojump.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -pbc whole

	echo "Protein" "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans 
  	
	echo "Protein" | $GMX trjconv -f $WTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $WTTrajDir/Rep$i/MD/md.tpr -o $WTTrajDir/Rep$i/FixPBC/frame0.pdb -dump 0
	
	rm $WTTrajDir/Rep$i/FixPBC/md_nojump.xtc
	rm $WTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc  	
done	


# L72S leptin simulations
for i in 1 2 3 4 5 6 7 8 9 10
do 
	cd $MTTrajDir/Rep$i
	mkdir $MTTrajDir/Rep$i/FixPBC
	
	cd $MTTrajDir/Rep$i/FixPBC

	echo "Protein" "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/MD/traj_comp.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/md_nojump.xtc -center -pbc nojump 
	
	echo "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/md_nojump.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -pbc whole

	echo "Protein" "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans 
  	
	echo "Protein" | $GMX trjconv -f $MTTrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $MTTrajDir/Rep$i/MD/md.tpr -o $MTTrajDir/Rep$i/FixPBC/frame0.pdb -dump 0
	
	rm $MTTrajDir/Rep$i/FixPBC/md_nojump.xtc
	rm $MTTrajDir/Rep$i/FixPBC/md_nojump_whole.xtc  	
done




# 2. Concatenate generated trajectories to generate an ensemble of possible protein conformation of each WT/MT simulations

# WT leptin simulations
cd $WorkDir
mkdir $WorkDir/WT
mkdir $WorkDir/L72S

cd $WorkDir/WT
mkdir $WorkDir/WT/ConcatenatedTraj
cd $WorkDir/WT/ConcatenatedTraj

$GMX trjcat -f $WTTrajDir/Rep*/FixPBC/md_nojump_whole_fit.xtc -o $WorkDir/WT/ConcatenatedTraj/trajout.xtc -cat -settime <<EOF
l
c
c
c
c
c
c
c
c
c
EOF


cp $WTTrajDir/Rep1/MD/md.tpr $WorkDir/WT/ConcatenatedTraj
cp $WTTrajDir/VaccumEM/VaccumEM.tpr $WorkDir/WT/ConcatenatedTraj
cp $WTTrajDir/PDB2GMX/Structure_pdb2gmx.pdb $WorkDir/WT/ConcatenatedTraj
cp $WTTrajDir/Rep1/FixPBC/frame0.pdb $WorkDir/WT/ConcatenatedTraj



# L72S leptin simulations
cd $WorkDir/L72S
mkdir $WorkDir/L72S/ConcatenatedTraj
cd $WorkDir/L72S/ConcatenatedTraj

$GMX trjcat -f $MTTrajDir/Rep*/FixPBC/md_nojump_whole_fit.xtc -o $WorkDir/L72S/ConcatenatedTraj/trajout.xtc -cat -settime <<EOF
l
c
c
c
c
c
c
c
c
c
EOF


cp $MTTrajDir/Rep1/MD/md.tpr $WorkDir/L72S/ConcatenatedTraj
cp $MTTrajDir/VaccumEM/VaccumEM.tpr $WorkDir/L72S/ConcatenatedTraj
cp $MTTrajDir/PDB2GMX/Structure_pdb2gmx.pdb $WorkDir/L72S/ConcatenatedTraj
cp $MTTrajDir/Rep1/FixPBC/frame0.pdb $WorkDir/L72S/ConcatenatedTraj



# 3. Collecting sampled side chain - side chain distance for Amino Acid Interaction Network analysis

# WT leptin simulations
cd $WorkDir/WT
mkdir $WorkDir/WT/SC_SC_distance
cd $WorkDir/WT/SC_SC_distance


# Collecting sampled distances along the trajectories
mkdir $WorkDir/WT/SC_SC_distance/SampledDistances
cd $WorkDir/WT/SC_SC_distance/SampledDistances


				##################################### ATTENTION ! ####################################################				
				# Modify '/path/to/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py' to correct file location
				
# Solvent accessible surface area (SASA): the following amino acids including  THR31, ILE35, ILE38, LEU72, MET75, THR78, LEU79, TYR82, ILE97, 
#LEU101, LEU104, LEU108, ALA146, LEU150, SER153, LEU154, and MET157 were found to be fully buried inside the WT leptin core. 
#Among these amino acids, seven (THR31, MET75, THR78, LEU79, LEU104, LEU108, and MET157) exhibit close proximity to the mutated residue LEU72.
#The sampled distance between interaction center of amino acid 72 and those of the listed amino acids over the course of simulations were 
#collected using the below analysis.


for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do 

	python3 /path/to/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $WorkDir/WT/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/WT/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i

done


# Collecting initial distances in the initial protein structure prior to the MD simulations
mkdir $WorkDir/WT/SC_SC_distance/InitialDistances
cd $WorkDir/WT/SC_SC_distance/InitialDistances

for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do 
	for file in $WTTrajDir/Rep*/FixPBC/frame0.pdb
	do
	
		python3 /path/to/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $file -tpr $WorkDir/WT/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i
		
	done

done



# L72S leptin simulations
cd $WorkDir/L72S
mkdir $WorkDir/L72S/SC_SC_distance
cd $WorkDir/L72S/SC_SC_distance

# Collecting sampled distances along the trajectories
mkdir $WorkDir/L72S/SC_SC_distance/SampledDistances
cd $WorkDir/L72S/SC_SC_distance/SampledDistances


for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do 

	python3 /path/to/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $WorkDir/L72S/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/L72S/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i

done


# Collecting initial distances in the initial protein structure prior to the MD simulations
mkdir $WorkDir/L72S/SC_SC_distance/InitialDistances
cd $WorkDir/L72S/SC_SC_distance/InitialDistances

for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do 
	for file in $MTTrajDir/Rep*/FixPBC/frame0.pdb
	do

		python3 /path/to/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $file -tpr $WorkDir/L72S/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i

	done
done


# 4. Collecting sampled volume of WT and MT leptin protein

# WT leptin simulations
cd $WorkDir/WT
mkdir $WorkDir/WT/ProteinVolumeMeasurement
cd $WorkDir/WT/ProteinVolumeMeasurement

# Create a .txt file that define interested region to compute volume
echo '22 167 A
' > $WorkDir/WT/ProteinVolumeMeasurement/InputFile.txt

				##################################### ATTENTION ! ####################################################				
				# Modify '/path/to/Measuring_protein_volume/Volume_computations.py' to correct file location

# Collecting sampled WT leptin structure volume values over the course of simulation
python3 /path/to/Measuring_protein_volume/Volume_computations.py -traj $WorkDir/WT/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/WT/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/WT/ProteinVolumeMeasurement/InputFile.txt -FileOut $WorkDir/WT/ProteinVolumeMeasurement/SampledVolumeValues.txt

# Collecting initial WT leptin structure volume value
for file in $WTTrajDir/Rep*/FixPBC/frame0.pdb
do
	python3 /path/to/Measuring_protein_volume/Volume_computations.py -traj $file -tpr $WorkDir/WT/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/WT/ProteinVolumeMeasurement/InputFile.txt -FileOut $WorkDir/WT/ProteinVolumeMeasurement/InitialVolumeValue.txt
	
done




# L72S leptin simulations
cd $WorkDir/L72S
mkdir $WorkDir/L72S/ProteinVolumeMeasurement
cd $WorkDir/L72S/ProteinVolumeMeasurement

# Create a .txt file that define interested region to compute volume
echo '22 167 A
' > $WorkDir/L72S/ProteinVolumeMeasurement/InputFile.txt


# Collecting sampled L72S leptin structure volume values over the course of simulation
python3 /path/to/Measuring_protein_volume/Volume_computations.py -traj $WorkDir/L72S/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/L72S/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/L72S/ProteinVolumeMeasurement/InputFile.txt -FileOut $WorkDir/L72S/ProteinVolumeMeasurement/SampledVolumeValues.txt

# Collecting initial L72S leptin structure volume value
for file in $MTTrajDir/Rep*/FixPBC/frame0.pdb
do
	python3 /path/to/Measuring_protein_volume/Volume_computations.py -traj $file -tpr $WorkDir/L72S/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/L72S/ProteinVolumeMeasurement/InputFile.txt -FileOut $WorkDir/L72S/ProteinVolumeMeasurement/InitialVolumeValue.txt
	
done



# 5. Measuring change in secondary structure of WT and MT leptin protein

# WT leptin simulations
cd $WorkDir/WT
mkdir $WorkDir/WT/SSChange
cd $WorkDir/WT/SSChange


				##################################### ATTENTION ! ####################################################				
				# Modify '/path/to/SecondaryStructure_Detection/SSDetection.py' to correct file location

# Collecting sampled WT leptin structure volume values over the course of simulation
python3 /path/to/SecondaryStructure_Detection/SSDetection.py -traj $WorkDir/WT/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/WT/ConcatenatedTraj/Structure_pdb2gmx.pdb -FileOut $WorkDir/WT/SSChange/SSOut



# L72S leptin simulations
cd $WorkDir/L72S
mkdir $WorkDir/L72S/SSChange
cd $WorkDir/L72S/SSChange


				##################################### ATTENTION ! ####################################################				
				# Modify '/path/to/SecondaryStructure_Detection/SSDetection.py' to correct file location

# Collecting sampled WT leptin structure volume values over the course of simulation
python3 /path/to/SecondaryStructure_Detection/SSDetection.py -traj $WorkDir/L72S/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/L72S/ConcatenatedTraj/Structure_pdb2gmx.pdb -FileOut $WorkDir/L72S/SSChange/SSOut

