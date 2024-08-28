for((i=22;i<168;i++))
do
	echo $i
	
	python3 Plotting_Distance_Distribution.py --WTSampledValue ../../WTLeptin/SC_SC_distance/SampledDistances/file_72_$i.txt --WTInitialValue ../../WTLeptin/SC_SC_distance/InitialDistances/File_72_$i.txt --MTSampledValue ../../L72SLeptin/SC_SC_distance/SampledDistances/file_72_$i.txt --MTInitialValue ../../L72SLeptin/SC_SC_distance/InitialDistances/File_72_$i.txt --FileOut Distance_Distribution_Plot/Dist72_$i

	python3 DistributionStatistics.py --SampledDistanceFile ../../WTLeptin/SC_SC_distance/SampledDistances/file_72_$i.txt --InitialDistance ../../WTLeptin/SC_SC_distance/InitialDistances/File_72_$i.txt --FileOut Distribution_Statistics/WT_72_$i
	
	python3 DistributionStatistics.py --SampledDistanceFile ../../L72SLeptin/SC_SC_distance/SampledDistances/file_72_$i.txt --InitialDistance ../../L72SLeptin/SC_SC_distance/InitialDistances/File_72_$i.txt --FileOut Distribution_Statistics/L72S_72_$i

done

