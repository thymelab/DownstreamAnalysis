import glob,shutil

for lmfile in glob.glob("*/lin*"):
	#print(lmfile)
	lmdir = lmfile.split("/")[0]
	filename0 = lmfile.split("/")[1]
	filename = lmfile.split("/")[1]
	#print(lmdir)
	#print(filename)
	genofile = open(lmdir + "/" + "genotyping",'r')
	for line in genofile:
		geno = line.strip('*').split(":")[0]
		#print(geno)
		if len(line.strip('*').split(":")) > 1:
			if len(line.strip('*').split(":")[1].split(',')) > 0:
				counts = len(line.strip('*').split(":")[1].strip().split(','))
				#print(counts)
				if "_" + geno + "_" in filename:
					#print(geno,counts)
					filename = filename.replace("_" + geno + "_", "_" + str(counts) + "-" + geno + "_")
	resultname = lmdir + "_" + filename
	#print(resultname)
	shutil.copy(lmdir + "/" + filename0, "./linearmodelfiles/" + resultname)
	#break
