import time


conditionsfile="/Users/allenma/Google Drive/pubgro/conditions_051216.txt"

def maketsv(conditionsfile):
	newlineordermatchednames = ["SRAnumber", "organism", "tissue", "general_celltype", "specific_celltype", "treatment_code","treated_or_like_treated" , "repnumber", "keyword", "exptype"]
	titleline=newlineordermatchednames+["mapped_reads", "total_reads", "percent_mapped"]
	f = open(conditionsfile, "U")
	orititleline = f.readline()
	orititlelinelist = orititleline.split("\t")
	indexlist = []
	for label in newlineordermatchednames:
		i = orititlelinelist.index(label)
		indexlist.append(i)
	NOW = time.strftime("%Y%m%d-%H%M%S")
	wf = open(conditionsfile+"_"+NOW+".csv", "w")
	wf.write(",".join(map(str,titleline))+"\n")
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		newline = []
		for i in indexlist:
			newline.append(line[i])
		numberofmappedreadsindex = orititlelinelist.index("numberofmappedreads")
		numberofreadsindex = orititlelinelist.index("numberofreads")
		mapped_reads = line[numberofmappedreadsindex].split(" ")[0]
		total_reads = line[numberofreadsindex].split(" ")[0]
		percent_mapped=line[numberofmappedreadsindex].split(" ")[-1]
		percent_mapped=percent_mapped.split(":")[0]
		percent_mapped=percent_mapped[1:-1]
		newline.append(mapped_reads)
		newline.append(total_reads)
		newline.append(percent_mapped)
		if mapped_reads!="":
			wf.write(",".join(map(str,newline))+"\n")
	f.close()
	wf.close()

maketsv(conditionsfile)
