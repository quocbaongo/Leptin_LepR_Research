import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def validate_file(f):

	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError('{0} does not exist'.format(f))
    
	return f


if __name__ == '__main__':

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Plotting distance distribution')
	parser.add_argument('--WTSampledValue', type=validate_file, help='File containing sampled values of WT structure in Å (.txt)', required=True)
	parser.add_argument('--WTInitialValue', type=validate_file, help='File containing initial values of WT structure in Å (.txt)', required=True)
	parser.add_argument('--MTSampledValue', type=validate_file, help='File containing sampled values of MT structure in Å (.txt)', required=True)
	parser.add_argument('--MTInitialValue', type=validate_file, help='File containing initial values of MT structure in Å (.txt)', required=True)	
	parser.add_argument('--FileOut', help='Specifying the name of out put file', required=True)
	args = parser.parse_args()
	
	WTVolumeValues = args.WTSampledValue
	WTInitialValues = args.WTInitialValue
	MTVolumeValues = args.MTSampledValue
	MTInitialValues = args.MTInitialValue	
	file_out = args.FileOut


	# WT Sampled volume values	
	WTVolumeValuesFile=open(WTVolumeValues,'r')
	WTVolumeValuesFileContent=WTVolumeValuesFile.readlines()
	WTVolumeValuesFileContent=[float(i.strip().split()[1]) for i in WTVolumeValuesFileContent]
	
	# WT Initial values
	WTAvgInitialValue=0.0
	Total=0
	
	WTInitialValuesFile=open(WTInitialValues, 'r')
	WTInitialValuesFileContent=WTInitialValuesFile.readlines()
	
	for line in WTInitialValuesFileContent:
		Total+=1
		WTAvgInitialValue+=(float(line.strip().split()[1]))
	
	WTAvgInitialValue/=Total
	
	
	# MT Sampled volume values	
	MTVolumeValuesFile=open(MTVolumeValues,'r')
	MTVolumeValuesFileContent=MTVolumeValuesFile.readlines()
	MTVolumeValuesFileContent=[float(i.strip().split()[1]) for i in MTVolumeValuesFileContent]
	
	# MT Initial values
	MTAvgInitialValue=0.0
	Total=0
	
	MTInitialValuesFile=open(MTInitialValues, 'r')
	MTInitialValuesFileContent=MTInitialValuesFile.readlines()
	
	for line in MTInitialValuesFileContent:
		Total+=1
		MTAvgInitialValue+=(float(line.strip().split()[1]))
	
	MTAvgInitialValue/=Total	
	
	

	# Plotting
	fig, ax = plt.subplots(figsize=(8, 4))
	WTweights = np.ones_like(np.array(WTVolumeValuesFileContent)) / len(WTVolumeValuesFileContent)
	MTweights = np.ones_like(np.array(MTVolumeValuesFileContent)) / len(MTVolumeValuesFileContent)

	# Merge data to find min max value
	merged_list = WTVolumeValuesFileContent + MTVolumeValuesFileContent
	MinValue=min(merged_list)
	MaxValue=max(merged_list)	
	
	
	nWT,binsWT,patchesWT=plt.hist(np.array(WTVolumeValuesFileContent),bins=100,
                                alpha=0.3,color='blue', weights=WTweights, label='WT',
                                range=[MinValue-(MinValue*0.05), MaxValue+(MaxValue*0.05)])
                                
	nMT,binsMT,patchesMT=plt.hist(np.array(MTVolumeValuesFileContent),bins=100,
                                alpha=0.7,color='red', weights=MTweights, label='MT',
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
	plt.xlabel('Volume ($\mathrm{\AA^3}$)',labelpad=10,fontsize=15)
	plt.ylabel('Probability density',labelpad=10,fontsize=15)
	plt.legend(loc='upper right')
	plt.savefig(f'{file_out}.png',dpi=600,bbox_inches='tight')























