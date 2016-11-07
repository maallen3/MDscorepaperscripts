import glob
import sys
import os

#previousconditionsfile = "/scratch/Shares/dowell/pubgro/conditions20151211-121018.txt"
previousconditionsfile = "/scratch/Shares/dowell/pubgro/conditions_20160613-1.txt"
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
	                bamfileroot = bamfileroot.strip(".flagstat")
        	        print bamfileroot, mapped_reads
                	mil_dic[bamfileroot] =   mapped_reads
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
                        keyword, bamfiledir,bamfile,posbedgraph,negbedgrah,tdffile,miniSRAnum, SRAtitle, organism, celltype,treatment_code,repnumber,treatmentnotes,exptype,whichprotocol,numberofreads, numberofmappedreads, intronexonratio,othernotes = line
                        if organism=="human":
				if bamfiledir not in bamdirtouse:
					if bamfiledir not in excludedir:
						bamdirtouse.append(bamfiledir)			


	return bamdirtouse

#list_of_directory_of_flagstat = get_list_of_directory_of_flagstat(previousconditionsfile)
#mil_dic = millions_dic_from_flagstat(list_of_directory_of_flagstat)
#dirofcouts = "/scratch/Shares/dowell/counts/human20151211-182056/"
dirofcouts = "/scratch/Shares/dowell/counts/human20160929-114122/" 

def createfpkmfile():
	res = {}
	wf = open("countsandfpkm.txt", "a")
	col = 7#if counting starting at 1
	for file_name in glob.glob(os.path.join(dirofcouts, '*.bed')):
	    for line_nr, line in enumerate(open(file_name)):
        	res.setdefault(line_nr, []).append(line.split('\t')[col-1])

	for line_nr in sorted(res):
    		wf.write('\t'.join(res[line_nr])+"\n")

createfpkmfile()
