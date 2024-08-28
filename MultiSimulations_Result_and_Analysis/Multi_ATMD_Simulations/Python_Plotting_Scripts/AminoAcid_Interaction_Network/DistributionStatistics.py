import os
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import statistics

def validate_file(f):

	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    
	return f


if __name__ == "__main__":

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Analyzing distance distribution curve')
	parser.add_argument("--SampledDistanceFile", type=validate_file, help="File containing sampled distance (.txt)", required=True)
	parser.add_argument("--InitialDistance", type=validate_file, help="File containing initial distance (.txt)", required=True)
	parser.add_argument("--FileOut", help="Specifying the name of out put file", required=True)
	args = parser.parse_args()
	
	SampledDistanceValues = args.SampledDistanceFile
	InitialDistanceValues = args.InitialDistance
	file_out = args.FileOut

	# File input processing
	SampledDistanceFileContent=open(SampledDistanceValues, "r")
	SampledDistanceFileContent=SampledDistanceFileContent.readlines()
	SampledDistanceFileContent=[float(i.strip().split()[1]) for i in SampledDistanceFileContent]


	InitialDistanceFileContent=open(InitialDistanceValues, "r")
	InitialDistanceFileContent=InitialDistanceFileContent.readlines()
	InitialDistanceFileContent=[float(i.strip().split()[1]) for i in InitialDistanceFileContent]

	MeanInitialValue=np.mean(InitialDistanceFileContent)


	###################################################################################################################################################################
	
	# Calculate mean value
	MeanValue=0.0
	TotalValue=0.0
	
	n,bins,patches=plt.hist(SampledDistanceFileContent,bins=100)
	
	for i in range(len(n)):
		TotalValue+=n[i]
		
	for i in range(len(bins)-1):
		binCenter=(bins[i]+bins[i+1])/2
		MeanValue+=(binCenter*(n[i]/TotalValue))
		
	Var=statistics.variance(SampledDistanceFileContent)
	
	# Statistics
	# Count population within cutoff range
	# Acceptance range of 5%
	AcceptanceRange=0.05
	
	UpperLimit=MeanInitialValue+(MeanInitialValue*AcceptanceRange)
	LowerLimit=MeanInitialValue-(MeanInitialValue*AcceptanceRange)	
		
	UpIdx=np.digitize(UpperLimit, bins)-1
	LowIdx=np.digitize(LowerLimit, bins)-1
	
	# Population within cutoff range
	TotalPopulation=0
	
	for i in range(LowIdx,UpIdx+1):
		TotalPopulation+=n[i]
		
	NormalizedTotalPopulation=TotalPopulation/TotalValue
		
	# Data to be written
	dictionary = {
			"MeanDistance": round(MeanValue,2),
			"InitialDistance": round(MeanInitialValue, 2),
			"Variance": round(Var, 2),
			"InitialDistancePopulation": round(NormalizedTotalPopulation,2)
			}
			
	#print(dictionary)
	# Serializing json
	json_object = json.dumps(dictionary, indent=4)
 
	# Writing to sample.json
	with open(f"{file_out}.json", "w") as outfile:
		outfile.write(json_object)
    
    

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
