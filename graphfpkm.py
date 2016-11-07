from __future__ import division
import numpy as np
import glob
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import keyfiles
from keyfiles import * #this is where i bring in all the locations of the files and folders I use



def get_active_inactive_all_SRR(TF, TYPE=0):
	all, active, notactive = [], [], []
	"""TYPE=0 is non-promoter. TYPE=1 is promoter, TYPE=2 is both."""
	for infile in MDscorefiles:
	        f = open(infile)
                SRR = infile.split("/")[-1]
                SRR = SRR.split("_")[0]
                for line in f:
                        line = line.strip("\n")
                        line = line.split("\t")
                        if len(line)==6:
                                name = line[0].split("_")[0]
                                if name==TF:
##motif identifier      Total(NON,TSS,COMB)     MD score(NON,TSS,COMB)  Stationary p-value(NON,TSS,COMB)        Non-stationary p-value(NON,TSS,COMB)    Null Expectation, Variance
#NRF1_HUMAN.H10MO.A     288,3422,3710   0.489583,0.667738,0.653908      1.000000,1.000000,1.000000      0.000000,0.000000,0.000000      0.424508,0.008080 
      

                                        spval = line[3]
                                        nspval = line[4]
#                                        nonspval, tssspval, combspval = map(float, spval.split(","))
#                                        nonnspval, tssnspval, combnspval = map(float, nspval.split(","))
					spval = map(float, spval.split(","))
					typespval = spval[TYPE]
					nspval = map(float, nspval.split(","))
					typenspval = nspval[TYPE]
					all.append(SRR)
                                        if typespval>0.95:
                                                if typenspval<0.05:
                                                        active.append(SRR)
                                        if typenspval>0.95:
                                                if typespval<0.05:
                                                        notactive.append(SRR)
                f.close()
        return [sorted(all), sorted(active), sorted(notactive)]



def readfpkms():
	f = open(maxfpkmfile)
	titleline = f.readline()
	titleline = titleline.strip("\n")
	titleline = titleline.split("\t")
	SRRindex = {}
	TFfpkmdic = {}
	for i, SRR in enumerate(titleline):
		if not SRR in SRRindex.keys():
			SRRindex[SRR] = i-1
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		TF = line[0]
		fpkms = line[1:]
		TFfpkmdic[TF]=fpkms
	return SRRindex, TFfpkmdic

SRRindex, TFfpkmdic = readfpkms()


def getfpkmsfromSRR(TF, SRRlists):
	theindexes = []
	fpkmlists = []
	if TF in TFfpkmdic.keys():
		allfpkmsthisTF = TFfpkmdic[TF]
		print "number of SRR", [len(li) for li in SRRlists]
		for SRRlist in SRRlists:
			indexes = []
			for SRR in SRRlist:
				if SRR in SRRindex.keys():
					if SRRindex[SRR] not in indexes:
						indexes.append(SRRindex[SRR])
				else:
					print "missing SRR", SRR
		#some SRR are missing... I need to check which
		
			theindexes.append(indexes)
		print "number of fpkms", [len(li) for li in theindexes]
		for indexlist in theindexes:
			fpkmlist = []
			for i in indexlist:
				fpkm = float(allfpkmsthisTF[i])
				fpkmlist.append(fpkm)
			fpkmlists.append(np.array(fpkmlist))
		return np.array(fpkmlists)
	else:
		print "No fpkm data for this TF. See the findfpkm file for more info"
		a = np.array([])
		print np.array([a for SRRlist in SRRlists])
		return False


def boxwhisoneTF_all_active_inactive(TF):
	SRRlabels = ["all", "active", "inactive"]
	SRRlists = get_active_inactive_all_SRR(TF)
	data_to_plot = getfpkmsfromSRR(TF,SRRlists)
	if data_to_plot!=False:
		outfilename = TF+"fpkmbox.png"
		boxwiskersofpkms(outfilename, data_to_plot, SRRlabels)

def boxwhisallTF_all_active_inactive():
	SRRlabels = ["all", "active", "inactive"]
	allTFs = getTFs()
	finaldata_to_plot = 0
	for TF in allTFs:
		print TF
		SRRlists = get_active_inactive_all_SRR(TF)
		data_to_plot = getfpkmsfromSRR(TF,SRRlists)
		if data_to_plot!=False:
			if finaldata_to_plot==0:
				finaldata_to_plot = data_to_plot
			else:
				finaldata_to_plot = [np.concatenate([finaldata_to_plot[i], data_to_plot[i]]) for i in range(len(finaldata_to_plot))]
	outfilename = "all_fpkmbox.png"
	boxwiskersofpkms(outfilename, finaldata_to_plot, SRRlabels)

	
def meds_all_file(outfile="medTFfpkms.txt"):
	medfoldchange = []
	allTFs = getTFs()
	wf = open(cwd+outfile, "w")
        titleline = ["TF", "median fpkm all", "median fpkm active", "median fpkm inactive", "fpkm fold change", "n = acitve", "n = inactive"]
	wf.write("\t".join(map(str, titleline))+"\n")
	for TF in allTFs:
		print TF
		SRRlists = get_active_inactive_all_SRR(TF)
                data_to_plot = getfpkmsfromSRR(TF,SRRlists)
		if data_to_plot!=False:
			all, activefpkms, notactivefpkms = data_to_plot
			if len(activefpkms)>0 and len(activefpkms)>0:
				fpkmfoldchange = np.median(activefpkms)/np.median(notactivefpkms)
				line = [TF, np.median(all), np.median(activefpkms), np.median(notactivefpkms), fpkmfoldchange, len(activefpkms), len(notactivefpkms), min(activefpkms), min(notactivefpkms)]
			elif len(activefpkms)>0:
				line = [TF, np.median(all), np.median(activefpkms), "nan","nan", len(activefpkms), len(notactivefpkms), min(activefpkms), "nan"]
			elif len(notactivefpkms)>0:
                                line = [TF, np.median(all), "nan", np.median(notactivefpkms),"nan", len(activefpkms), len(notactivefpkms), "nan", min(notactivefpkms)]
			wf.write("\t".join(map(str, line))+"\n")
	wf.close()




def requireminfpkm_graphsub(fpkmcutoff=0, minsamples=0):
        f = open(outfile="medTFfpkms.txt")
        f.readline()#skip title
        medfoldchange = []
        for line in f:
                line = line.strip("\n")
                line = line.split("\t")
                if line[4]!="nan":
                        if line[4]!='0.0':
                                if line[4]!="inf":
					if float(line[2])>fpkmcutoff:
						if float(line[3])>fpkmcutoff:
							if int(line[5])>minsamples:
								if int(line[6])>minsamples:
		                                        #num = np.log(float(line[3]))
									num = np.log(float(line[1]))-np.log(float(line[2]))	
									medfoldchange.append(num)
        fig = plt.figure(1, figsize=(9, 6))
        # Create an axes instance        
        ax = fig.add_subplot(111)
        bins = range(-50, 50, 1)
        bins = [a/10 for a in bins]
        print bins
        ax.hist(medfoldchange, bins=bins)
        fig.savefig(cwd+"graphs/submedfoldchangeminfpkm"+str(fpkmcutoff)+"minsamples"+str(minsamples)+".png", bbox_inches='tight')
	print "graphs/submedfoldchangeminfpkm"+str(fpkmcutoff)+"minsamples"+str(minsamples)+".png"
        plt.clf()




def boxwiskersofpkms(outfilename, data_to_plot, catlabels):
	fig = plt.figure(1, figsize=(9, 6))
	# Create an axes instance
	ax = fig.add_subplot(111)
	# Create the boxplot
	bp = ax.boxplot(data_to_plot)
	SRRlabels = [label+" n="+str(len(data_to_plot[i])) for i, label in enumerate(catlabels)]
	ax.set_xticklabels(SRRlabels)
	#plt.ylim([0,3])
	plt.ylabel("fpkm")
	plt.yscale('log')
	print "saving fig ", "graphs/"+outfilename
	fig.savefig("graphs/"+outfilename, bbox_inches='tight')
	plt.clf()


		
	



#boxwhisoneTF_all_active_inactive("SOX2")
#boxwhisoneTF_all_active_inactive("P53")
#boxwhisoneTF_all_active_inactive("ESR1")
#boxwhisoneTF_all_active_inactive("NF2L1")
#boxwhisallTF_all_active_inactive()
#boxwhisoneTF_all_active_inactive("TF2LX")
#boxwhisallTF_all_active_inactive()
#meds_all_file(outfile="medTFfpkms.txt")
#requireminfpkm_graphsub(fpkmcutoff=0, minsamples=0)

