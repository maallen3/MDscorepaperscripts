from __future__ import division
import numpy as np
import glob
import sys
import os


#files and folders needed
cwd = os.getcwd()+"/"
othergenenamelist = cwd+"GeneName_manual2.txt"
translator = cwd+"TF.bed"
maxfpkmfile = cwd + "allTF_maxfpkm.txt"
isoformoffpkmfile = cwd+"allTF_maxisoform.txt"


dir_of_MDscores = "/scratch/Shares/dowell/MD_score_files/"
fileoffpkm="/scratch/Shares/dowell/counts/human20160929-114122/countsandfpkm.txt"
MDscorefiles = [infile for infile in glob.glob(os.path.join(dir_of_MDscores, '*.tsv'))]
aMDscorefile = MDscorefiles[0]
refgene = "/scratch/Shares/dowell/pubgro/genomefiles/human/refGene.bed"



def getTFs():
	f = open(aMDscorefile)
	names = []
	for line in f:
		line = line.strip()
		line = line.split("\t")
		if len(line)==6:
			nameline = line[0]
			name = nameline.split("_")[0]
			if not name.startswith("#"):
				if name not in names:
					names.append(name)
	return sorted(names)


