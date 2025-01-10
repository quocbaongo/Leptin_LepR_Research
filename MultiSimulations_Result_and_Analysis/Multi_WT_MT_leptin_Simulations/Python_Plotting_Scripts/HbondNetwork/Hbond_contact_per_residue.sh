# WT

mkdir WTRes_HbondNetwork
cd WTRes_HbondNetwork

for id in 45 46 47 48 49 50 51 52 53 54 55 56 57 58 64
do

	python3 ../Atom_atom_hydrogenBondLifeTime.py ../../../WTLeptin/HbondNetwork/${id}A/HydrogenBondPerTimeStep.txt ../../../WTLeptin/ConcatenatedTraj/Structure_pdb2gmx.pdb ${id}A	> Output_${id}A.txt

done

cd ..

# P64S

mkdir P64SRes_HbondNetwork
cd P64SRes_HbondNetwork

for id in 45 46 47 48 49 50 51 52 53 54 55 56 57 58 64
do

	python3 ../Atom_atom_hydrogenBondLifeTime.py ../../../P64SLeptin/HbondNetwork/${id}A/HydrogenBondPerTimeStep.txt ../../../P64SLeptin/ConcatenatedTraj/Structure_pdb2gmx.pdb ${id}A > Output_${id}A.txt

done

cd ..
