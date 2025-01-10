#!/bin/bash -l

# Atomistic molecular dynamics simulations analysis
# Required input
WorkDir=/path/to/analysis/directory
AnalysisScripts=/path/to/generated/Python_Analysis_Scripts
TrajDir=/path/to/generated/trajectories
GMX=gmx_2021.2


# The simulation of WT/MT leptin were repeated 10 times
#Thus, we have 10 generated trajectories for each
#1. First step of the anaylsis workflow is to correct for periodicity

for i in 1 2 3 4 5 6 7 8 9 10
do
        cd $TrajDir/Rep$i
        mkdir $TrajDir/Rep$i/FixPBC

        cd $TrajDir/Rep$i/FixPBC

        echo "Protein" "Protein" | $GMX trjconv -f $TrajDir/Rep$i/MD/traj_comp.xtc -s $TrajDir/Rep$i/MD/md.tpr -o $TrajDir/Rep$i/FixPBC/md_nojump.xtc -center -pbc nojump

        echo "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/md_nojump.xtc -s $TrajDir/Rep$i/MD/md.tpr -o $TrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -pbc whole

        echo "Protein" "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/md_nojump_whole.xtc -s $TrajDir/Rep$i/MD/md.tpr -o $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -fit rot+trans

        echo "Protein" | $GMX trjconv -f $TrajDir/Rep$i/FixPBC/md_nojump_whole_fit.xtc -s $TrajDir/Rep$i/MD/md.tpr -o $TrajDir/Rep$i/FixPBC/frame0.pdb -dump 0

        rm $TrajDir/Rep$i/FixPBC/md_nojump.xtc
        rm $TrajDir/Rep$i/FixPBC/md_nojump_whole.xtc
done



# 2. Concatenate generated trajectories to generate an ensemble of possible protein conformation of each WT/MT simulations
cd $WorkDir
mkdir $WorkDir/ConcatenatedTraj
cd $WorkDir/ConcatenatedTraj

$GMX trjcat -f $TrajDir/Rep*/FixPBC/md_nojump_whole_fit.xtc -o $WorkDir/ConcatenatedTraj/trajout.xtc -cat -settime <<EOF
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


cp $TrajDir/Rep1/MD/md.tpr $WorkDir/ConcatenatedTraj
cp $TrajDir/VaccumEM/VaccumEM.tpr $WorkDir/ConcatenatedTraj
cp $TrajDir/PDB2GMX/Structure_pdb2gmx.pdb $WorkDir/ConcatenatedTraj


# 3. Collecting sampled volume of WT and MT leptin protein

mkdir $WorkDir/ProteinVolumeMeasurement
cd $WorkDir/ProteinVolumeMeasurement

# Create a .txt file that define interested region to compute volume, which is the area enclosed by the four major α-helices within leptin structure
#The main four anti-parallel α-helices include helix A (from residue Pro23 to His47), B (from residue Leu72 to Ser88), C (from residue Arg92 to Lys115), and D (from Ser141 to Ser164).

echo '23 47 A
72 88 A
92 115 A
141 164 A
' > $WorkDir/ProteinVolumeMeasurement/Helices.txt

# Collecting sampled WT/MT leptin helices structure volume values over the course of simulation
python3 $AnalysisScripts/Measuring_protein_volume/Volume_computations.py -traj $WorkDir/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/ProteinVolumeMeasurement/Helices.txt -FileOut $WorkDir/ProteinVolumeMeasurement/SampledHelicesVolume.txt

# Collecting initial WT/MT leptin structure volume value
for file in $TrajDir/Rep*/FixPBC/frame0.pdb
do

        python3 $AnalysisScripts/Measuring_protein_volume/Volume_computations.py -traj $file -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -Region $WorkDir/ProteinVolumeMeasurement/Helices.txt -FileOut $WorkDir/ProteinVolumeMeasurement/InitialHelicesVolume.txt

done

# 4. RMSF

# WT leptin simulations
mkdir $WorkDir/RMSFCA
cd $WorkDir/RMSFCA

echo "C-alpha" | $GMX rmsf -f $WorkDir/ConcatenatedTraj/trajout.xtc -s $WorkDir/ConcatenatedTraj/VaccumEM.tpr -o $WorkDir/RMSFCA/rmsfCA.xvg

# 5. Hbonding analysis
# Hydrogen bond analysis was conducted on the mutation site (Pro64) and on each amino acid at the N-terminal of the AB loop (residues Ile45 to Thr58), a region exhibiting increased flexibility in the MT leptin

cd $WorkDir/HbondNetwork

# From residues Ile45 to Thr58
for ((id=45;id<59;id++))
do
        for chain in A
        do
                mkdir $WorkDir/HbondNetwork/${id}${chain}
                cd $WorkDir/HbondNetwork/${id}${chain}

                python3 $AnalysisScripts/HbondNetworkAnalysis/DetectingHydrogenBond.py $WorkDir/ConcatenatedTraj/trajout.xtc $WorkDir/ConcatenatedTraj/VaccumEM.tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb $id $chain

                python3 $AnalysisScripts/HbondNetworkAnalysis/ReadJsonFile.py $WorkDir/HbondNetwork/${id}${chain}/HydrogenBondDetection.json $WorkDir/ConcatenatedTraj/VaccumEM.tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb

                python3 $AnalysisScripts/HbondNetworkAnalysis/DataForPlotting.py $WorkDir/HbondNetwork/${id}${chain}/HydrogenBondPerTimeStep.txt $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb $id $chain

        done

done

# Pro64/Ser64
for id in 64
do
        for chain in A
        do
                mkdir $WorkDir/HbondNetwork/${id}${chain}
                cd $WorkDir/HbondNetwork/${id}${chain}

                python3 $AnalysisScripts/HbondNetworkAnalysis/DetectingHydrogenBond.py $WorkDir/ConcatenatedTraj/trajout.xtc $WorkDir/ConcatenatedTraj/VaccumEM.tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb $id $chain

                python3 $AnalysisScripts/HbondNetworkAnalysis/ReadJsonFile.py $WorkDir/HbondNetwork/${id}${chain}/HydrogenBondDetection.json $WorkDir/ConcatenatedTraj/VaccumEM.tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb

                python3 $AnalysisScripts/HbondNetworkAnalysis/DataForPlotting.py $WorkDir/HbondNetwork/${id}${chain}/HydrogenBondPerTimeStep.txt $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb $id $chain

        done

done


# 6. Collecting sampled side chain - side chain distance for Amino Acid Interaction Network analysis

# WT/MT leptin simulations
cd $WorkDir
mkdir $WorkDir/SC_SC_distance
cd $WorkDir/SC_SC_distance


# Collecting sampled distances along the trajectories
mkdir $WorkDir/SC_SC_distance/SampledDistances
cd $WorkDir/SC_SC_distance/SampledDistances

# Solvent accessible surface area (SASA): the following amino acids including  THR31, ILE35, ILE38, LEU72, MET75, THR78, LEU79, TYR82, ILE97,
#LEU101, LEU104, LEU108, ALA146, LEU150, SER153, LEU154, and MET157 were found to be fully buried inside the WT leptin core.
#Among these amino acids, seven (THR31, MET75, THR78, LEU79, LEU104, LEU108, and MET157) exhibit close proximity to the mutated residue LEU72.
#The sampled distance between interaction center of amino acid 72 and those of the listed amino acids over the course of simulations were
#collected using the below analysis.

for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do
        mkdir $WorkDir/SC_SC_distance/SampledDistances/72_$i
        cd $WorkDir/SC_SC_distance/SampledDistances/72_$i
        python3 $AnalysisScripts/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $WorkDir/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i -ChainID1 A -ChainID2 A

done


# Collecting initial distances in the initial protein structure prior to the MD simulations
mkdir $WorkDir/SC_SC_distance/InitialDistances
cd $WorkDir/SC_SC_distance/InitialDistances

for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do
        mkdir $WorkDir/SC_SC_distance/InitialDistances/72_$i
        for file in $TrajDir/Rep*/FixPBC/frame0.pdb
        do
                cd $WorkDir/SC_SC_distance/InitialDistances/72_$i

                python3 $AnalysisScripts/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $file -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 72 -resi2 $i -ChainID1 A -ChainID2 A

        done

done



# The center-of-geometry distances between the side chains of Leu66 and each of hydrophobic residues including Ile63, Leu66, Leu79, Val81, Ile85, Leu125, Leu131, Val134, Leu135, Val145, Ala146, and Leu150) were monitored throughout the WT and MT concatenated trajectories

#The sampled distance between interaction center of Leu66 and those of the hydrophobic amino acids over the course of simulations were collected using the below analysis.
cd $WorkDir/SC_SC_distance/SampledDistances

for i in 63 64 79 81 85 125 131 134 135 145 146 150
do
        mkdir $WorkDir/SC_SC_distance/SampledDistances/${i}_66
        cd $WorkDir/SC_SC_distance/SampledDistances/${i}_66
        
        python3 $AnalysisScripts/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $WorkDir/ConcatenatedTraj/trajout.xtc -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 66 -resi2 $i -ChainID1 A -ChainID2 A

done

# Collecting initial distances in the initial protein structure prior to the MD simulations
cd $WorkDir/SC_SC_distance/InitialDistances

for i in 63 64 79 81 85 125 131 134 135 145 146 150
do
        mkdir $WorkDir/SC_SC_distance/InitialDistance/${i}_66
        cd $WorkDir/SC_SC_distance/InitialDistance/${i}_66

        for file in $TrajDir/Rep*/FixPBC/frame0.pdb
        do

                python3 $AnalysisScripts/AminoAcid_Interaction_Network/Measuring_SCToSC_Distance.py -traj $file -tpr $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -resi1 $id -resi2 66 -ChainID1 A -ChainID2 A

        done

done


# 7. Measuring change in secondary structure of WT and MT leptin protein

cd $WorkDir
mkdir $WorkDir/SSChange
cd $WorkDir/SSChange

echo "Protein" | $GMX do_dssp -f $WorkDir/ConcatenatedTraj/trajout.xtc -s $WorkDir/ConcatenatedTraj/Structure_pdb2gmx.pdb -o $WorkDir/SSChange/ss.xpm

















