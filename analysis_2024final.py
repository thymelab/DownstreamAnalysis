#!/usr/bin/python
import glob,argparse,os,itertools
from itertools import chain

parser = argparse.ArgumentParser(description='options for sorting p-values')
parser.add_argument('-statsfile', type=str, action="store", dest="statsfile", default="*linearmodel*out")
#parser.add_argument('-statsfile', type=str, action="store", dest="statsfile", default="fam222a-ab*linearmodel*out")
parser.add_argument('-ofile', type=str, action="store", dest="ofile", default="")
parser.add_argument('-include', type=str, action="store", dest="includetxt", default="ribgraph")
parser.add_argument('-exclude', type=str, action="store", dest="excludetxt", default="weirdword")
parser.add_argument('-cutoff', type=float, action="store", dest="cutoff", default=0.05)
parser.add_argument('-overwriteanova', action="store_true", dest="overwriteanova", default=False) # Use the pval from linear mixed model instead of from anova # I am using anova data for everything here, since it goes lower than the lmm
parser.add_argument('-showfulldata', action="store_true", dest="showfulldata", default=True) # would almost always want to get the files back so they are all ready for analysis
parser.add_argument('-sortalpha', action="store_true", dest="sortalpha", default=False)

args = parser.parse_args()
statsfile = args.statsfile
ofile = args.ofile
includetxt = args.includetxt
excludetxt = args.excludetxt
cutoff = args.cutoff
overwriteanova = args.overwriteanova
showfulldata = args.showfulldata
sortalpha = args.sortalpha # default will be to sort by lowest pval, but this is another option

gene_pval = {}
gene_pval2 = {}
gene_pval3 = {}
# Example name
# ap1g1-ME_Box6_08_16_2022_JM_me_linearmodel_41-ap1g1Box6_het_vs_22-ap1g1Box6_hom_15870640_4294967294
for sfilename in glob.glob(statsfile):
	pvallist = []
	insiglist = []
	signames = []
	gene = sfilename.split("_")[0] + "_" + sfilename.split("linearmodel")[1].split("_")[2] + "_" + sfilename.split("linearmodel")[1].split("_")[5]
	sfile = open(sfilename, 'r')
	if ofile == "":
		ofile = "pvals_" + sfilename
	outfile = open(ofile, 'w')
	lines = sfile.readlines()
	for line in lines:
		# Example string
		#anova:  ribgraph_mean_time_day2dfall_numberofbouts_3600_controlgroup-het.data : N of control, test, Mean of array control, test, Variance of array control, test, SSMD, H-stat, P-value:  43 35 1196.7162790697676 1193.72 272040.68322336394 316867.31017142854 0.0039044380458149457 0.10660666975064714 0.7440409692821573
		if line.startswith("anova:"): # automatically skipping all the failed ones
			pval = line.strip().split()[-1] # I always leave pval at the end
			ssmd = line.strip().split()[-3] # I always leave pval at the end
			graph = "_".join([line.strip().split(":")[1].strip().split("_")[0]] + line.strip().split(":")[1].strip().split("_")[3:])
			if "a1%84%92" in graph: # deals with comparing old date that had slightly different timing for the ppi
				graph = graph.replace("a1%84%92","a1%89%97")
			fulldata = line.strip().split(":")[3].strip()
			pvallist.append([float(pval), graph, fulldata, float(ssmd)])
		else:
			if not overwriteanova: # don't even look at the linear mixed model info
				continue
			if line.startswith("ribgraph"):
				ribgraph = line.strip() # saving and then it will match up with the next instance of hitting the "mutornot"
				if "a1%84%92" in ribgraph:
					ribgraph = ribgraph.replace("a1%84%92","a1%89%97")
			if line.startswith("linear model failed"):
				ribgraph = line.split(":")[1].strip() # saving and then it will match up with the next instance of hitting the "mutornot"
				ribgraph = "_".join([line.strip().split(":")[1].strip().split("_")[0]] + line.strip().split(":")[1].strip().split("_")[3:])
				pvallist.append([float(1000.0), ribgraph, ""])
			if line.startswith("mutornot[T.wt] "):
				if len(line.split()) > 3:
					lmmpval = line.split()[4]
					coef = line.split()[1]
					if float(lmmpval) == 0:
						lmmpval = 0.001
					#if float(lmmpval) < cutoff:
					ribgraph = "_".join([line.strip().split(":")[1].strip().split("_")[0]] + line.strip().split(":")[1].strip().split("_")[3:])
					pvallist.append([float(lmmpval), ribgraph, ""])
	for stat in pvallist:
		if stat[2] == "": # if it's lmm data
			for s in range(0, len(pvallist)):
				if stat[1] == pvallist[s][1]: # found the lmm data
					pvallist[s][0] = stat[0] # replacing pval in the non-lmm data with lmm pval
	filteredpvallist0 = [x for x in pvallist if not x[2]==""] # eliminate all the original lmms
	filteredpvallistf = filteredpvallist0
	filteredpvallist = [x for x in filteredpvallist0 if not x[0]>cutoff] # filter by the cutoff
	if sortalpha:
		filteredpvallist.sort(key=lambda x: x[1])
		filteredpvallistf.sort(key=lambda x: x[1])
	else:
		filteredpvallist.sort(key=lambda x: x[0])
		filteredpvallistf.sort(key=lambda x: x[0])
	if showfulldata:
		for final in filteredpvallistf: # using the full pvalue list here for show full data
			if includetxt in str(final[1]):
				if excludetxt not in str(final[1]):
					# BASICALLY REWRITING ORIGINAL FILE, BUT CLEANED UP AND SORTED BY PVAL
					outfile.write("anova: " + str(final[1]) + " " + str(final[2]) + '\n')
					signames.append(str(final[1]))
					if final[0]>cutoff:
						insiglist.append(str(final[1]))
	else:
		for final in filteredpvallist:
			if includetxt in str(final[1]):
				if excludetxt not in str(final[1]):
					outfile.write(str(final[0]) + " " + str(final[1]) + '\n')
					signames.append(str(final[1]))
	ofile = ""
	if gene in gene_pval:
		gene_pval[gene].append(filteredpvallist)
		gene_pval2[gene].append(insiglist)
		gene_pval3[gene].append(signames)
	else:
		gene_pval[gene] = []
		gene_pval2[gene] = []
		gene_pval3[gene] = []
		gene_pval[gene].append(filteredpvallist)
		gene_pval2[gene].append(insiglist)
		gene_pval3[gene].append(signames)

for k,v in gene_pval.items():
	shared0 = {}
	print(k) # gene + the genetic comparison
	print(len(v)) # how many runs done
	if len(v) == 1:
		continue
	outfile2 = open("replicates_" + k + ".out", 'w')
	grouped_lists = {} #pulls all of the same behavior out from all of the trials
	for outer_list in v:
		for inner_list in outer_list:
			#print(inner_list)
			newkey = inner_list[1]
			#print(newkey)
			if newkey in grouped_lists:
				grouped_lists[newkey].append(inner_list)
			else:
				grouped_lists[newkey] = []
				grouped_lists[newkey].append(inner_list)
	for newkey,group in grouped_lists.items():
		if len(group) > 1:
			sorted_group = sorted(group, key=lambda x: x[0]) # this puts the smallest p-value first
			first_entry_sign = 1 if sorted_group[0][3] >= 0 else -1
			result = None
			signmatch = False
			for i in range(1, len(sorted_group)):
				entry_sign = 1 if sorted_group[i][3] >= 0 else -1
				if entry_sign == first_entry_sign:
					result = sorted_group[0][0]
					signmatch = True
					if newkey not in shared0.keys():
						shared0[newkey] = [(result,sorted_group[0][3])]
					else:
						print("THERE IS A PROBLEM")
					break
			if len(group) > 2 and signmatch == False:
				result = sorted_group[1][0] # necessarily the second one has to match signs with any later ones, since it did not match with the first
				if newkey not in shared0.keys():
					shared0[newkey] = [(result,sorted_group[1][3])]
				else:
					print("SECOND THERE IS A PROBLEM")
			else: # nothing matched, so don't keep
				continue
	for k2,v2 in shared0.items():
		outfile2.write("anova: " + k2 + " " +  str(v2[0][1]) + " sigspace " + str(v2[0][0]) + "\n")
	fulllist = list(set(itertools.chain.from_iterable(gene_pval2[k])))
	for l in fulllist:
		if l not in shared0.keys():
			outfile2.write("anova: " + l  + " 0 insignificant 2\n")
	outfile2.close()
