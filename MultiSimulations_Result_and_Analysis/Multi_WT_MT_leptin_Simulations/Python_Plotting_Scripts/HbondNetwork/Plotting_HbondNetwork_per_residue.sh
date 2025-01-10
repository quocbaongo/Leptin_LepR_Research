mkdir WT_P64S_HbondNetwork_Compare
cd WT_P64S_HbondNetwork_Compare

for id in 45 46 47 48 49 50 51 52 53 54 55 56 57 58 64
do

	python3 ../PlottingHbondNetwork.py --WTInput ../../../WTLeptin/HbondNetwork/${id}A/HydrogenBondCountDict.json --WTTotalFrame 100001 --MTInput ../../../P64SLeptin/HbondNetwork/${id}A/HydrogenBondCountDict.json --MTTotalFrame 100001 --WTPDB ../Structure_pdb2gmx.pdb --Output ${id}A --y_lower_limit -0.25 --y_upper_limit 0.25

done
