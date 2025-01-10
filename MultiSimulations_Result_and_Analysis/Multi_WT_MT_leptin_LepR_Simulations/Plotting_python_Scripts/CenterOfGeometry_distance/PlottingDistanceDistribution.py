import os
import argparse
import numpy as np
import pandas as pd
import statistics
import matplotlib.pyplot as plt

def validate_file(f):

	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    
	return f


if __name__ == "__main__":

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Plotting distance values distribution')
	parser.add_argument("--WTSampledValue", type=validate_file, help="File containing sampled distance values of WT structure, unit of nm (.txt)", required=True)
	parser.add_argument("--WTInitialValue", type=validate_file, help="File containing initial distance values of WT structure, unit of nm (.txt)", required=True)
	parser.add_argument("--MTSampledValue", type=validate_file, help="File containing sampled distance values of MT structure, unit of nm (.txt)", required=True)
	parser.add_argument("--MTInitialValue", type=validate_file, help="File containing initial distance values of MT structure, unit of nm (.txt)", required=True)	
	parser.add_argument("--FileOut", help="Specifying the name of out put file", required=True)
	args = parser.parse_args()
	
	WTVolumeValues = args.WTSampledValue
	WTInitialValues = args.WTInitialValue
	MTVolumeValues = args.MTSampledValue
	MTInitialValues = args.MTInitialValue	
	file_out = args.FileOut


	# WT Sampled volume values	
	WTVolumeValuesFile=open(WTVolumeValues,"r")
	WTVolumeValuesFileContent=WTVolumeValuesFile.readlines()
	WTVolumeValuesFileContent=[float(i.strip().split()[1]) * 10.0 for i in WTVolumeValuesFileContent]
	
	
	# WT Initial values
	WTAvgInitialValue=0.0
	Total=0
	
	WTInitialValuesFile=open(WTInitialValues, "r")
	WTInitialValuesFileContent=WTInitialValuesFile.readlines()
	
	for line in WTInitialValuesFileContent:
		Total+=1
		WTAvgInitialValue+=(float(line.strip().split()[1]) * 10.0)
	
	WTAvgInitialValue/=Total
	
	
	# MT Sampled volume values	
	MTVolumeValuesFile=open(MTVolumeValues,"r")
	MTVolumeValuesFileContent=MTVolumeValuesFile.readlines()
	MTVolumeValuesFileContent=[float(i.strip().split()[1]) * 10.0 for i in MTVolumeValuesFileContent]
	
	# MT Initial values
	MTAvgInitialValue=0.0
	Total=0
	
	MTInitialValuesFile=open(MTInitialValues, "r")
	MTInitialValuesFileContent=MTInitialValuesFile.readlines()
	
	for line in MTInitialValuesFileContent:
		Total+=1
		MTAvgInitialValue+=(float(line.strip().split()[1]) * 10.0)
	
	MTAvgInitialValue/=Total	
	
	# Merge data to find min max value
	merged_list = WTVolumeValuesFileContent + MTVolumeValuesFileContent
	MinValue=min(merged_list)
	MaxValue=max(merged_list)	

	# Plotting
	fig, ax = plt.subplots(figsize=(8, 4))
	WTweights = np.ones_like(np.array(WTVolumeValuesFileContent)) / len(WTVolumeValuesFileContent)
	MTweights = np.ones_like(np.array(MTVolumeValuesFileContent)) / len(MTVolumeValuesFileContent)	
	
	nWT,binsWT,patchesWT=plt.hist(np.array(WTVolumeValuesFileContent),bins=100,
                                alpha=0.3,color="blue", weights=WTweights, label="WT",
                                range=[MinValue-(MinValue*0.05), MaxValue+(MaxValue*0.05)])
                                
	nMT,binsMT,patchesMT=plt.hist(np.array(MTVolumeValuesFileContent),bins=100,
                                alpha=0.7,color="red", weights=MTweights, label="MT",
                                range=[MinValue-(MinValue*0.05), MaxValue+(MaxValue*0.05)])                                
                                
	# Compute the center of each bin for plotting
	WTbin_centers = (binsWT[1:] + binsWT[:-1]) / 2                                
	MTbin_centers = (binsMT[1:] + binsMT[:-1]) / 2
	
	# Initial value bar height
	InitialValueHeight=0.0
	
	if max(nWT) > max(nMT):
		InitialValueHeight=max(nWT)
		
	else:
		InitialValueHeight=max(nMT)	
	
	InitialValueHeight=InitialValueHeight+ (InitialValueHeight*0.2)
                                

	# Initial value bar width
	InitialValueWidth=0.0
    
	widthWT=binsWT[1]-binsWT[0]
	widthMT=binsMT[1]-binsMT[0]		
    
	if widthWT > widthMT:
		InitialValueWidth=widthWT
	else:
		InitialValueWidth=widthMT	

	InitialValueWidth/=2                                
                                
	plt.plot(WTbin_centers, nWT, linewidth=2, color='blue')
	plt.plot(MTbin_centers, nMT, linewidth=2, color='red')                                
                                

	plt.bar([WTAvgInitialValue],[InitialValueHeight],color='green',width=InitialValueWidth,
                                            label='WT starting volume value')
	plt.bar([MTAvgInitialValue],[InitialValueHeight],color='orange',width=InitialValueWidth,
                                            label='MT starting volume value')

	
	plt.tight_layout()
	plt.xlabel('Distance (Å)',labelpad=10,fontsize=15)
	plt.ylabel('Probability density',labelpad=10,fontsize=15)
	plt.legend(loc='upper right')
	plt.savefig(f'{file_out}.png',dpi=600,bbox_inches='tight')
	#plt.savefig(f'{file_out}.eps',dpi=600,bbox_inches='tight')
	
	
	# Generate distribution statistics
	# Calculate mean value
	# WT
	WTMeanValue=0.0
	WTTotalValue=0.0
	
	
	for i in range(len(nWT)):
		WTTotalValue+=nWT[i]
		
	for i in range(len(binsWT)-1):
		WTbinCenter=(binsWT[i]+binsWT[i+1])/2
		WTMeanValue+=(WTbinCenter*(nWT[i]/WTTotalValue))
		
	WTVar=statistics.variance(WTVolumeValuesFileContent)
	

	
	# MT
	MTMeanValue=0.0
	MTTotalValue=0.0
	
	
	for i in range(len(nMT)):
		MTTotalValue+=nMT[i]
		
	for i in range(len(binsMT)-1):
		MTbinCenter=(binsMT[i]+binsMT[i+1])/2
		MTMeanValue+=(MTbinCenter*(nMT[i]/MTTotalValue))
		
	MTVar=statistics.variance(MTVolumeValuesFileContent)
	
	print(len(WTVolumeValuesFileContent), len(MTVolumeValuesFileContent))	
	
	print(f"WTMeanDistance: {round(WTMeanValue,3)}")
	print(f"WTInitialDistance: {round(WTAvgInitialValue, 3)}")
	print(f"WTVariance: {round(WTVar, 3)}")	
	
	
	print(f"MTMeanDistance: {round(MTMeanValue,3)}")
	print(f"MTInitialDistance: {round(MTAvgInitialValue, 3)}")
	print(f"MTVariance: {round(MTVar, 3)}")	
	
	
	
	
	
	
