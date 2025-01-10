#!/bin/sh

mkdir WT_L72S
mkdir WT_L72S/Distance_Distribution_Plot
mkdir WT_L72S/Distribution_Statistics

mkdir WT_P64S
mkdir WT_P64S/Distance_Distribution_Plot
mkdir WT_P64S/Distribution_Statistics


for i in 31 35 38 75 78 79 82 97 101 104 108 146 150 153 154 157
do

	echo $i
	cd WT_L72S
	
	
	python3 Plotting_Distance_Distribution.py --WTSampledValue ../../WTLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --WTInitialValue ../../WTLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --MTSampledValue ../../L72SLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --MTInitialValue ../../L72SLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distance_Distribution_Plot/Dist72_$i

	python3 DistributionStatistics.py --SampledDistanceFile ../../WTLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --InitialDistance ../../WTLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distribution_Statistics/WT_72_$i

	python3 DistributionStatistics.py --SampledDistanceFile ../../L72SLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --InitialDistance ../../L72SLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distribution_Statistics/L72S_72_$i



	cd WT_P64S

	python3 Plotting_Distance_Distribution.py --WTSampledValue ../../WTLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --WTInitialValue ../../WTLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --MTSampledValue ../../P64SLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --MTInitialValue ../../P64SLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distance_Distribution_Plot/Dist72_$i

	python3 DistributionStatistics.py --SampledDistanceFile ../../WTLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --InitialDistance ../../WTLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distribution_Statistics/WT_72_$i

	python3 DistributionStatistics.py --SampledDistanceFile ../../P64SLeptin/SC_SC_distance/SampledDistances/72_$i/file_72_$i.txt --InitialDistance ../../P64SLeptin/SC_SC_distance/InitialDistances/72_$i/file_72_$i.txt --FileOut Distribution_Statistics/P64S_72_$i

done
