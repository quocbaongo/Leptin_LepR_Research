#!/bin/sh

for i in 63 64 79 81 85 125 131 134 135 145 146 150
do

	echo $i

	cd WT_P64S

	python3 Plotting_Distance_Distribution.py --WTSampledValue ../../WTLeptin/HydrophobicNetwork/SampledDistance/${i}_66/file_${i}_66.txt --WTInitialValue ../../WTLeptin/HydrophobicNetwork/InitialDistance/${i}_66/file_${i}_66.txt --MTSampledValue ../../P64SLeptin/HydrophobicNetwork/SampledDistance/${i}_66/file_${i}_66.txt --MTInitialValue ../../P64SLeptin/HydrophobicNetwork/InitialDistance/${i}_66/file_${i}_66.txt --FileOut Distance_Distribution_Plot/Dist${i}_66

	python3 DistributionStatistics.py --SampledDistanceFile ../../WTLeptin/HydrophobicNetwork/SampledDistance/${i}_66/file_${i}_66.txt --InitialDistance ../../WTLeptin/HydrophobicNetwork/InitialDistance/${i}_66/file_${i}_66.txt --FileOut Distribution_Statistics/WT_${i}_66

	python3 DistributionStatistics.py --SampledDistanceFile ../../P64SLeptin/HydrophobicNetwork/SampledDistance/${i}_66/file_${i}_66.txt --InitialDistance ../../P64SLeptin/HydrophobicNetwork/InitialDistance/${i}_66/file_${i}_66.txt --FileOut Distribution_Statistics/P64S_${i}_66

done
