# Required input
WorkDir=/path/to/working/directory
MDPFile=/path/to/mdp/files/directory
Structure=/path/to/simulated/structure
GMX=gmx_2021.2

# Generates system topology

mkdir $WorkDir/PDB2GMX
cd $WorkDir/PDB2GMX

$GMX pdb2gmx -f $Structure -o $WorkDir/PDB2GMX/Structure_pdb2gmx.pdb -p $WorkDir/PDB2GMX/topol.top -ignh -water spce <<EOF
6	
EOF

# The AMBER99SB-ILDN protein force-field was chosen


#Vaccum energy minimization
mkdir $WorkDir/VaccumEM
cd $WorkDir/VaccumEM

$GMX editconf -f $WorkDir/PDB2GMX/Structure_pdb2gmx.pdb -o $WorkDir/VaccumEM/box.pdb -bt cubic -d 1.0 -c

$GMX grompp -f $MDPFile/em.mdp -c $WorkDir/VaccumEM/box.pdb -r $WorkDir/VaccumEM/box.pdb -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/VaccumEM/VaccumEM.tpr -maxwarn 1

srun $GMX mdrun -s $WorkDir/VaccumEM/VaccumEM.tpr


mkdir $WorkDir/Solvation
cd $WorkDir/Solvation

$GMX solvate -cp $WorkDir/VaccumEM/confout.gro -cs spc216.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Solvation/solv.gro

$GMX grompp -f $MDPFile/ions.mdp -c $WorkDir/Solvation/solv.gro -r $WorkDir/Solvation/solv.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Solvation/ions.tpr -maxwarn 1

echo "SOL" | $GMX genion -s $WorkDir/Solvation/ions.tpr -o $WorkDir/Solvation/solv_ions.gro -p $WorkDir/PDB2GMX/topol.top -pname NA -nname CL -neutral -conc 0.15 


# Create directory for each independent simulation
for((i=1;i<11;i++))
do
	mkdir $WorkDir/Rep$i
done

# Simulated molecules + solvent energy minimization

for((i=1;i<11;i++))
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/EM
	cd $WorkDir/Rep$i/EM
	
	$GMX grompp -f $MDPFile/em.mdp -c $WorkDir/Solvation/solv_ions.gro -r $WorkDir/Solvation/solv_ions.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Rep$i/EM/em.tpr
	
done


# Running EM on parallel
cd $WorkDir 
srun $GMX mdrun -s em.tpr -multidir $WorkDir/Rep*/EM



# First phase of equilibration conducted under an NVT ensemble (constant Number of particles, Volume, and Temperature): Each simulation was initiated with a random seed for initial velocity.


for((i=1;i<11;i++))
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/NVT
	cd $WorkDir/Rep$i/NVT
	
	$GMX grompp -f $MDPFile/nvt.mdp -c $WorkDir/Rep$i/EM/confout.gro -r $WorkDir/Rep$i/EM/confout.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Rep$i/NVT/nvt.tpr
		
done


# Running NVT equilibration on parallel
cd $WorkDir 
srun $GMX mdrun -s nvt.tpr -multidir $WorkDir/Rep*/NVT


# Second phase of equilibration conducted under an NPT ensemble (constant Number of particles, Pressure, and Temperature)
for((i=1;i<11;i++))
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/NPT
	cd $WorkDir/Rep$i/NPT
	
	$GMX grompp -f $MDPFile/npt.mdp -c $WorkDir/Rep$i/NVT/confout.gro -r $WorkDir/Rep$i/NVT/confout.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Rep$i/NPT/npt.tpr
		
done

# Running NPT equilibration on parallel
cd $WorkDir 
srun $GMX mdrun -s npt.tpr -multidir $WorkDir/Rep*/NPT


# Production
for((i=1;i<11;i++))
do
	cd $WorkDir/Rep$i
	mkdir $WorkDir/Rep$i/MD
	cd $WorkDir/Rep$i/MD

	$GMX grompp -f $MDPFile/md.mdp -c $WorkDir/Rep$i/NPT/confout.gro -r $WorkDir/Rep$i/NPT/confout.gro -p $WorkDir/PDB2GMX/topol.top -o $WorkDir/Rep$i/MD/md.tpr

done

# Running production on parallel
cd $WorkDir 
srun $GMX mdrun -s md.tpr -multidir $WorkDir/Rep*/MD








