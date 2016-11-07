from __future__ import division
import glob
import sys
import os

previousconditionsfile = "/scratch/Shares/dowell/pubgro/conditions_20160613-1.txt"
#previousconditionsfile = "/scratch/Shares/dowell/pubgro/conditions20151211-121018.txt"
excludedir=["/scratch/Shares/dowell/pubgro/Meng2014b/mapped/bowtie2/sortedbam/"] #Meng2014b sort did not work

def millions_dic_from_flagstat(list_of_directory_of_flagstat):
        mil_dic = {}
	for directory_of_flagstat in list_of_directory_of_flagstat:
	        for flagstat in glob.glob(os.path.join(directory_of_flagstat, '*sorted.bam.flagstat')):
			print flagstat
        	        f = open(flagstat, "U")
	                line_with_mapped_reads = 2
        	        lines = f.readlines()
                	f.close()
	                mapped_reads = lines[line_with_mapped_reads].split(" ")[0]
        	        mapped_reads = int(mapped_reads)
                	bamfileroot = flagstat.split("/")[-1]
			if bamfileroot.startswith("ERR769523"):
				print "ERR769523", bamfileroot
			#problem is you only want undersore 1
	                bamfileroot = bamfileroot.strip(".flagstat")
        	        print bamfileroot, mapped_reads
                	mil_dic[bamfileroot] =   mapped_reads
			if bamfileroot.find("flip")>0:
				blist = bamfileroot.split(".flip")
				newbamfileroot = blist[0]+blist[1]
				mil_dic[newbamfileroot] =   mapped_reads
				
#			if bamfileroot.startswith("ERR"):
#				bamfilerootb = bamfileroot.split("_")[1]
#				if bamfilerootb.startswith("1"):
#					mil_dic[bamfileroot+"_1"] =   mapped_reads
        return mil_dic


def get_list_of_directory_of_flagstat(previousconditionsfile):
        pf = open(previousconditionsfile, "U")
        bamdirtouse = []
        for oline in pf:
                line = oline.strip("\n")
                line = line.split("\t")
                if line[0]=="":
                        break
                checkforprocressing = line[1]
                if not checkforprocressing.startswith("this"):
                        #keyword, bamfiledir,bamfile,posbedgraph,negbedgrah,tdffile,miniSRAnum, SRAtitle, organism, celltype,treatment_code,repnumber,treatmentnotes,exptype,whichprotocol,numberofreads, numberofmappedreads, intronexonratio,othernotes = line
#
#                        if organism=="human":
#				if bamfiledir not in bamdirtouse:
#					if bamfiledir not in excludedir:
#						bamdirtouse.append(bamfiledir)			

			organism = line[8]
                        bamfiledir = line[1]
                        bamfile = line[2]
                        if organism=="human":
				if bamfiledir not in bamdirtouse:
					if bamfiledir not in excludedir:
                                               bamdirtouse.append(bamfiledir) 

	return bamdirtouse

list_of_directory_of_flagstat = get_list_of_directory_of_flagstat(previousconditionsfile)
mil_dic = millions_dic_from_flagstat(list_of_directory_of_flagstat)
print mil_dic
for key in sorted(mil_dic.keys()):
	if key.startswith("SRR296"):
		print "this the key you were looking for", key
#dirofcouts = "/scratch/Shares/dowell/counts/human20151211-182056/"
dirofcouts = "/scratch/Shares/dowell/counts/human20160929-114122/"

def createfpkmfile():
	res = {}
	lineinfo = {}
	millmapped,labellist = [], []
	wf = open(dirofcouts+"countsandfpkm.txt", "w")
	col = 7#if counting starting at 1
	print "num files", len(glob.glob(os.path.join(dirofcouts, '*.bed')))
	for file_name in glob.glob(os.path.join(dirofcouts, '*.bed')):
		milmappedval = file_name.split("/")[-1]
		milmappedval = milmappedval.split("_")[0]
		print milmappedval
		if milmappedval.startswith("ERR"):
			millmapped.append(mil_dic[milmappedval+"_1.fastqbowtie2.sorted.bam"])
		else:
			millmapped.append(mil_dic[milmappedval])
	
		labellist.append(milmappedval)
	    	if len(lineinfo)==0:
			for line_nr, line in enumerate(open(file_name)):
				lineinfo.setdefault(line_nr, []).append(line.strip("\n").split('\t')[0:6])
		for line_nr, line in enumerate(open(file_name)):
        		res.setdefault(line_nr, []).append(line.strip("\n").split('\t')[col-1])
	titleline = ["chr", "start", "stop", "name", "score", "strand"]+[label+"_reads" for label in labellist]+[label+"_fpkm "+str(millmapped[i]) for i, label in enumerate(labellist)]
	print "len title line", len(titleline)
	print "len labellist", len(labellist)
	print "len info", 6
	wf.write("\t".join(titleline)+"\n")
	for line_nr in sorted(res):
		info = lineinfo[line_nr]
		info = info[0]
		chr, start, stop, name, score, strand = info 
		fpkms = []
		print "the counts", len(res[line_nr])
		for i, num_of_reads in enumerate(res[line_nr]):
			total_reads = millmapped[i]
			fragment_size = int(stop)-int(start)
			fpkm = calfpkm(num_of_reads, total_reads, fragment_size)	
			fpkms.append(fpkm)
#		print lineinfo[line_nr], res[line_nr], millmapped, fpkms
		wline = info+res[line_nr]+fpkms
		print name, len(wline)
		if len(wline)!=len(titleline):
			print "something is wrong with one of your bed files. Check the number of lines in the bed files you are using."
			break
		wf.write("\t".join(map(str,wline))+"\n")
		#wf.write(wline)
		#just need to write it to file now....

#	for line_nr in sorted(res):
#    		wf.write('\t'.join(res[line_nr])+"\n")


def calfpkm(num_of_reads, total_reads, fragment_size):
        fpkm = int(num_of_reads)/((int(fragment_size)/1000)*(int(total_reads)/1000000))
        return fpkm


createfpkmfile()
