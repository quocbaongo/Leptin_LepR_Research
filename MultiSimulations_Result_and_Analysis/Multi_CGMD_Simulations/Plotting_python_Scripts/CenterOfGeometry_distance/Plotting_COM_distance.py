import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


if __name__ == "__main__":

	WTWorkDir='../../WTLeptin_LR/CenterOfGeometry_distance'
	MTWorkDir='../../L72SLeptin_LR/CenterOfGeometry_distance'

	# WT
	WTFileABContents=[]
	WTFileCDContents=[]
	WTFileEFContents=[]
	
	for i in range(1,6):
	
		# FileAB
		with open(f"{WTWorkDir}/Rep{i}/FileAB.txt", 'r') as f:
			FileABContent = f.readlines()
            
		FileABContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileABContent]
		FileABContent=sorted(FileABContent, key=lambda x: x[0])
		
		WTFileABContents.append(FileABContent)
		
		# FileCD
		with open(f"{WTWorkDir}/Rep{i}/FileCD.txt", 'r') as f:
			FileCDContent = f.readlines()
            
		FileCDContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileCDContent]
		FileCDContent=sorted(FileCDContent, key=lambda x: x[0])
		
		WTFileCDContents.append(FileCDContent)      
    
		# FileEF
		with open(f"{WTWorkDir}/Rep{i}/FileEF.txt", 'r') as f:
			FileEFContent = f.readlines()
            
		FileEFContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileEFContent]
		FileEFContent=sorted(FileEFContent, key=lambda x: x[0])
		
		WTFileEFContents.append(FileEFContent)


	# MT
	MTFileABContents=[]
	MTFileCDContents=[]
	MTFileEFContents=[]

	for i in range(1,6):
	
		# FileAB
		with open(f"{MTWorkDir}/Rep{i}/FileAB.txt", 'r') as f:
			FileABContent = f.readlines()
            
		FileABContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileABContent]
		FileABContent=sorted(FileABContent, key=lambda x: x[0])
        
		MTFileABContents.append(FileABContent)
            
		# FileCD
		with open(f"{MTWorkDir}/Rep{i}/FileCD.txt", 'r') as f:
			FileCDContent = f.readlines()
            
		FileCDContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileCDContent]
		FileCDContent=sorted(FileCDContent, key=lambda x: x[0])
        
		MTFileCDContents.append(FileCDContent)      
    
		# FileEF
		with open(f"{MTWorkDir}/Rep{i}/FileEF.txt", 'r') as f:
			FileEFContent = f.readlines()
            
		FileEFContent=[[int(i.strip().split()[0])/100,float(i.strip().split()[1])/10] for i in FileEFContent]
		FileEFContent=sorted(FileEFContent, key=lambda x: x[0])
        
		MTFileEFContents.append(FileEFContent)

	########################################## Start Plotting ######################################################
	# Plotting
	# Define size of figure
	fig = plt.figure(figsize=(14, 7))
	gs = gridspec.GridSpec(42, 10)
	
	# Define the positions of the subplots.
	ax0=fig.add_subplot(gs[:8, :5])
	ax1=fig.add_subplot(gs[9:17, :5])
	ax2=fig.add_subplot(gs[18:26, :5])
	ax3=fig.add_subplot(gs[27:34, :5])
	ax4=fig.add_subplot(gs[35:42, :5])

	ax0x=fig.add_subplot(gs[:8, 5:])
	ax1x=fig.add_subplot(gs[9:17, 5:])
	ax2x=fig.add_subplot(gs[18:26, 5:])
	ax3x=fig.add_subplot(gs[27:34, 5:])
	ax4x=fig.add_subplot(gs[35:42, 5:])

	# Plotting data
	# WT
	for j in range(5):
		globals()['ax%s' % j].plot([i[0] for i in WTFileABContents[j]], [i[1] for i in WTFileABContents[j]], color='green', label='Leptin1-LepR1')
		globals()['ax%s' % j].plot([i[0] for i in WTFileCDContents[j]], [i[1] for i in WTFileCDContents[j]], color='red', label='Leptin2-LepR2')
		globals()['ax%s' % j].plot([i[0] for i in WTFileEFContents[j]], [i[1] for i in WTFileEFContents[j]], color='blue', label='Leptin3-LepR3')
		
		# Set text to image
		globals()['ax%s' % j].text(0.5, 0.80, f"Rep{j+1}", transform = globals()['ax%s' % j].transAxes, fontsize=11)
		# Set y axis limit
		globals()['ax%s' % j].set_ylim([0.0, 2.1])
        
		# Set legend
		globals()['ax%s' % j].legend(prop={'size': 6})
		# y ticks rotation
		globals()['ax%s' % j].set_yticklabels(np.array([round(i, 2) for i in globals()['ax%s' % j].get_yticks()]),rotation=15)

	# MT
	for j in range(5):
		globals()['ax%sx' % j].plot([i[0] for i in MTFileABContents[j]], [i[1] for i in MTFileABContents[j]], color='green', label='Leptin1-LepR1')
		globals()['ax%sx' % j].plot([i[0] for i in MTFileCDContents[j]], [i[1] for i in MTFileCDContents[j]], color='red', label='Leptin2-LepR2')
		globals()['ax%sx' % j].plot([i[0] for i in MTFileEFContents[j]], [i[1] for i in MTFileEFContents[j]], color='blue', label='Leptin3-LepR3')
		# Set text to image
		globals()['ax%sx' % j].text(0.5, 0.8, f"Rep{j+1}", transform = globals()['ax%sx' % j].transAxes, fontsize=11)
		# Set y axis limit
		globals()['ax%sx' % j].set_ylim([0.0, 2.1])
		# Set legend
		globals()['ax%sx' % j].legend(prop={'size': 6})
        
	# x,yticks
	ax0.xaxis.set_tick_params(labelbottom=False)	
	ax1.xaxis.set_tick_params(labelbottom=False)
	ax2.xaxis.set_tick_params(labelbottom=False)
	ax3.xaxis.set_tick_params(labelbottom=False)

	ax0x.xaxis.set_tick_params(labelbottom=False)	
	ax1x.xaxis.set_tick_params(labelbottom=False)
	ax2x.xaxis.set_tick_params(labelbottom=False)
	ax3x.xaxis.set_tick_params(labelbottom=False)

	ax0x.yaxis.set_tick_params(labelbottom=False)	
	ax1x.yaxis.set_tick_params(labelbottom=False)
	ax2x.yaxis.set_tick_params(labelbottom=False)
	ax3x.yaxis.set_tick_params(labelbottom=False)
	ax4x.yaxis.set_tick_params(labelbottom=False)
    
	# x,y axis label
	ax4.set_xlabel('Time (ns)',size=10, labelpad=6)
	ax4x.set_xlabel('Time (ns)',size=10, labelpad=6)

	ax0.set_ylabel('Distance (nm)',size=10, labelpad=20)
	ax1.set_ylabel('Distance (nm)',size=10, labelpad=20)
	ax2.set_ylabel('Distance (nm)',size=10, labelpad=20)
	ax3.set_ylabel('Distance (nm)',size=10, labelpad=20)
	ax4.set_ylabel('Distance (nm)',size=10, labelpad=20)
	
	# Set fig label
	ax0.set_title('A', loc='left', fontdict={'fontsize': 20})
	ax0x.set_title('B', loc='left', fontdict={'fontsize': 20})
	
	# Save figure
	fig.tight_layout()
	plt.savefig(f"DistanceEvolution.png",dpi=600,bbox_inches='tight')
	plt.savefig(f"DistanceEvolution.eps",dpi=600,bbox_inches='tight')

	
	
	
	
	
	
	
	
	
	
	
	
