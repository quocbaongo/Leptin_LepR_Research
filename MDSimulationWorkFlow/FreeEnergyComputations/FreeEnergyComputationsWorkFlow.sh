#!/bin/bash

# A general workflow that can be used for
#computing both leptin folding or binding free energy change upon L72S mutation

# Required input
WorkDir=/path/to/directory/for/computing/free/energy
MDPFile=/path/to/mdp/files/directory
Structure=/path/to/wild-type/bound/or/unbound/leptin/structure
GMX=gmx_2021.2

# Generating leptin hybrid structure containing both the wild-type and mutant amino acids simultaneously
mkdir $WorkDir/Mut_coordinate_topology
mkdir $WorkDir/Mut_coordinate_topology/PDB2GMX
cd $WorkDir/Mut_coordinate_topology/PDB2GMX

$GMX pdb2gmx -f $Structure -o $WorkDir/Mut_coordinate_topology/PDB2GMX/protein_pdb2gmx.pdb -p $WorkDir/Mut_coordinate_topology/PDB2GMX/topol.top -ff amber99sb-star-ildn-mut -water none -ignh


# Hybrid topology
mkdir $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol
cd $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol


			####################################################################
			########################### ATTENTION !!! ##########################
			####################################################################
			
# Change /path/to/mutate.py and /path/to/generate_hybrid_topology.py accordingly			
# Residue id is shifted to 1 => mutation residue id of 72 corresponds to id of 51 
#when using file mutate.py
		
python3 /path/to/mutate.py -f $WorkDir/Mut_coordinate_topology/PDB2GMX/protein_pdb2gmx.pdb -o $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/mut.pdb -ff amber99sb-star-ildn-mut <<EOF
51
Ser
n
EOF


$GMX pdb2gmx -f $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/mut.pdb -o $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/mut_pdb2gmx.pdb -ff amber99sb-star-ildn-mut -water tip3p -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/topol.top -i $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/posre.itp

python3 /path/to/generate_hybrid_topology.py -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/topol.top -o $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -ff amber99sb-star-ildn-mut



# Generate simulation box using protein hybrid structure
mkdir $WorkDir/Mut_coordinate_topology/solvation
cd $WorkDir/Mut_coordinate_topology/solvation

$GMX editconf -f $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/mut_pdb2gmx.pdb -o $WorkDir/Mut_coordinate_topology/solvation/box.pdb -bt dodecahedron -d 2.0 -c

$GMX solvate -cp $WorkDir/Mut_coordinate_topology/solvation/box.pdb -cs spc216.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Mut_coordinate_topology/solvation/water.pdb

$GMX grompp -f $MDPFile/genion.mdp -c $WorkDir/Mut_coordinate_topology/solvation/water.pdb -r $WorkDir/Mut_coordinate_topology/solvation/water.pdb -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Mut_coordinate_topology/solvation/tpr.tpr 

echo "SOL" | $GMX genion -s $WorkDir/Mut_coordinate_topology/solvation/tpr.tpr -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Mut_coordinate_topology/solvation/ions.pdb -conc 0.15 -pname NA -nname CL -neutral




			#############################################################################
			########################### Simulation replication ##########################
			#############################################################################	

# The whole procedure including energy minimization, equilibrium and non-equilibrium simulations were repeated 3 times

			
			# Equilibrium simulations
			

# Generating directory for each repetition

for i in 1 2 3 4 5
do
	mkdir $WorkDir/Rep$i
done
	
# Simulated molecules + solvent energy minimization
# StateA (hybrid amino acid is set at wild-type state)

for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/EM
	mkdir $WorkDir/Rep$i/EM/StateA
	cd $WorkDir/Rep$i/EM/StateA

	$GMX grompp -f $MDPFile/forward/equil_md/f_enmin.mdp -c $WorkDir/Mut_coordinate_topology/solvation/ions.pdb -r $WorkDir/Mut_coordinate_topology/solvation/ions.pdb -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/EM/StateA/tpr.tpr
	
done	

# StateB (hybrid amino acid is set at mutant state)
for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/EM/StateB
	cd $WorkDir/Rep$i/EM/StateB

	$GMX grompp -f $MDPFile/reverse/equil_md/r_enmin.mdp -c $WorkDir/Mut_coordinate_topology/solvation/ions.pdb -r $WorkDir/Mut_coordinate_topology/solvation/ions.pdb -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/EM/StateB/tpr.tpr
	
done	

# Running EM simulations on parallel
cd $WorkDir 
srun $GMX mdrun -s tpr.tpr -multidir $WorkDir/Rep*/EM/State*




# Equilibrating the simulated system
# StateA (hybrid amino acid is set at wild-type state)

for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/Equilibration
	mkdir $WorkDir/Rep$i/Equilibration/StateA
	cd $WorkDir/Rep$i/Equilibration/StateA

	$GMX grompp -f $MDPFile/forward/equil_md/f_npt.mdp -c $WorkDir/Rep$i/EM/StateA/confout.gro -r $WorkDir/Rep$i/EM/StateA/confout.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/Equilibration/StateA/tpr.tpr -maxwarn 1
	
done


# StateB (hybrid amino acid is set at mutant state)
for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/Equilibration/StateB
	cd $WorkDir/Rep$i/Equilibration/StateB

	$GMX grompp -f $MDPFile/reverse/equil_md/r_npt.mdp -c $WorkDir/Rep$i/EM/StateB/confout.gro -r $WorkDir/Rep$i/EM/StateB/confout.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/Equilibration/StateB/tpr.tpr -maxwarn 1
	
done	


# Running equilibrations on parallel
cd $WorkDir 
srun $GMX mdrun -s tpr.tpr -multidir $WorkDir/Rep*/Equilibration/State*


# Equilibrium simulation
# StateA (hybrid amino acid is set at wild-type state)
for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/MD
	mkdir $WorkDir/Rep$i/MD/StateA
	cd $WorkDir/Rep$i/MD/StateA
	
	$GMX grompp -f $MDPFile/forward/equil_md/f_equil.mdp -c $WorkDir/Rep$i/Equilibration/StateA/confout.gro -r $WorkDir/Rep$i/Equilibration/StateA/confout.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/MD/StateA/tpr.tpr

done


# StateB (hybrid amino acid is set at mutant state)
for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/MD/StateB
	cd $WorkDir/Rep$i/MD/StateB

	$GMX grompp -f $MDPFile/reverse/equil_md/r_equil.mdp -c $WorkDir/Rep$i/Equilibration/StateB/confout.gro -r $WorkDir/Rep$i/Equilibration/StateB/confout.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/MD/StateB/tpr.tpr
	
done


# Running Equilibrium simulations on parallel
cd $WorkDir 
srun $GMX mdrun -s tpr.tpr -multidir $WorkDir/Rep*/MD/State*


			# Non-Equilibrium simulations

# Prepare for non-equilibrium simulations
# State A (hybrid amino acid is set at wild-type state)

for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/FrameExtraction
	mkdir $WorkDir/Rep$i/FrameExtraction/StateA
	cd $WorkDir/Rep$i/FrameExtraction/StateA
	
	echo "System" | $GMX trjconv -f $WorkDir/Rep$i/MD/StateA/traj.trr -s $WorkDir/Rep$i/MD/StateA/tpr.tpr -o $WorkDir/Rep$i/FrameExtraction/StateA/frame_.gro -ur compact -pbc mol -b 20000 -skip 3 -sep 
	
done


# State B (hybrid amino acid is set at mutant state)

for i in 1 2 3 4 5
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/FrameExtraction/StateB
	cd $WorkDir/Rep$i/FrameExtraction/StateB
	
	echo "System" | $GMX trjconv -f $WorkDir/Rep$i/MD/StateB/traj.trr -s $WorkDir/Rep$i/MD/StateB/tpr.tpr -o $WorkDir/Rep$i/FrameExtraction/StateB/frame_.gro -ur compact -pbc mol -b 20000 -skip 3 -sep 
	
done


# Non-equilibrium transitions

for i in 1 2 3 4 5
do
	mkdir $WorkDir/Rep$i/NonEquilibriumTransition
	mkdir $WorkDir/Rep$i/NonEquilibriumTransition/StateA
	mkdir $WorkDir/Rep$i/NonEquilibriumTransition/StateB	

done	


# State A (hybrid amino acid is set at wild-type state)
for i in 1 2 3 4 5
do
	for ((j=0;j<100;j++))
	do
		mkdir $WorkDir/Rep$i/NonEquilibriumTransition/StateA/frame$j
		cd $WorkDir/Rep$i/NonEquilibriumTransition/StateA/frame$j

		$GMX grompp -f $MDPFile/forward/nonequil_md/f_nonequil.mdp -c $WorkDir/Rep$i/FrameExtraction/StateA/frame_$j.gro -r $WorkDir/Rep$i/FrameExtraction/StateA/frame_$j.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/NonEquilibriumTransition/StateA/frame$j/tpr.tpr -maxwarn 1
			
	done
done


# State B (hybrid amino acid is set at mutant state)
for i in 1 2 3 4 5
do
	for ((j=0;j<100;j++))
	do
		mkdir $WorkDir/Rep$i/NonEquilibriumTransition/StateB/frame$j
		cd $WorkDir/Rep$i/NonEquilibriumTransition/StateB/frame$j

		$GMX grompp -f $MDPFile/reverse/nonequil_md/r_nonequil.mdp -c $WorkDir/Rep$i/FrameExtraction/StateB/frame_$j.gro -r $WorkDir/Rep$i/FrameExtraction/StateB/frame_$j.gro -p $WorkDir/Mut_coordinate_topology/Mut_hybrid_topol/hybrid.top -o $WorkDir/Rep$i/NonEquilibriumTransition/StateB/frame$j/tpr.tpr -maxwarn 1
	
	done
done

# Running Non-Equilibrium transitions on parallel
cd $WorkDir 
srun $GMX mdrun -s tpr.tpr -multidir $WorkDir/Rep*/NonEquilibriumTransition/State*/frame*




































