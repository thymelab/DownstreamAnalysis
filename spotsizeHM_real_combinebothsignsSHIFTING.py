import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import glob
import pandas
import re

def is_single_line_file(file_path):
	with open(file_path, 'r') as file:
		lines = file.readlines()
		return len(lines) == 1


def custom_sort(name):
	#print(name)
	alphanumeric_part = re.sub(r'\d+-', '', name)
	#print(alphanumeric_part)
	return (alphanumeric_part, name)


#means2HM_ssmd_stimuli_Kiaa0232_Box12_4_21BC_linearmodel_kiaa0232Box12_wt_vs_kiaa0232Box12_hom_11693279_4294967294
#,Assay Type,SSMD_x,Significance,Percents,SSMD_y
prefixlist = []
filelist = []
#posmeansepsign2HM_ssmd_baseline_replicates_atp1a3a-ko_wt_hom.csv
#posmeansepsign2HM_ssmd_stimuli_pvals_atp1a3a-ko_Box2_07_03_20_GK_cosscz_linearmodel_18-atp1a3aBox2_wt_vs_20-atp1a3aBox2_het_7357888_4294967294.csv
#old way above and below, new above
#for mfile0 in glob.glob("means2HM_ssmd_stimuli_Kiaa0232_Box12_4_21BC_linearmodel_kiaa0232Box12_hetandwt_vs_kiaa0232Box12_hom_11693280_4294967294*"):
for mfile0 in glob.glob("*sepsign*2HM_ssmd_*replicates*csv"):
	if is_single_line_file(mfile0):
		continue
	if "_linearmodel_" in mfile0:
		prefix = mfile0.split("_linearmodel_")[0].split("_ssmd_")[1]
	else:
		prefix = "_".join(mfile0.split(".")[0].split("_ssmd_")[1].split("_")[:-2])
	filelist.append(mfile0)
	prefixlist.append(prefix)
filelist = sorted(filelist, key=custom_sort)
prefixlist.sort()
for prefix in prefixlist:
	#print("PREFIX:",prefix)
	s = []
	ns = []
	c = []
	nc = []
	ylabels = []
	N = 0
	for mfile in filelist:
		#print(mfile)
		if "_linearmodel_" in mfile:
			premfile = mfile.split("_linearmodel_")[0].split("_ssmd_")[1]
		else:
			premfile = "_".join(mfile.split(".")[0].split("_ssmd_")[1].split("_")[:-2])
		if prefix != premfile:
			continue
		csv = pandas.read_csv(mfile)
		csv = csv.loc[:, ~csv.columns.str.contains('^Unnamed')]
		if "_linearmodel_" in mfile:
#posmeansepsign2HM_ssmd_stimuli_pvals_atp1a3a-ko_Box2_07_03_20_GK_cosscz_linearmodel_18-atp1a3aBox2_wt_vs_20-atp1a3aBox2_het_7357888_4294967294.csv
			gene = mfile.split("_linearmodel_")[0].split("_pvals_")[1].split("_")[0] + "_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[1] + "_vs_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[4]
		else:
#posmeansepsign2HM_ssmd_baseline_replicates_atp1a3a-ko_wt_hom.csv
			gene = mfile.split("_replicates_")[1].split(".")[0].split("_")[0] + "_" + mfile.split("_replicates_")[1].split(".")[0].split("_")[1] + "_vs_" + mfile.split("_replicates_")[1].split(".")[0].split("_")[2]
		if("neg" in mfile):
			ns.append(np.array(csv["Percents"].tolist())) # Significance is the % of the assays with significance, not sure why it ended up being named that)
			nc.append(np.array(csv["SSMD_y"].tolist()))
			ylabels.append(gene)
		if("pos" in mfile):
			N = N + 1
			s.append(np.array(csv["Percents"].tolist())) # Significance is the % of the assays with significance, not sure why it ended up being named that)
			c.append(np.array(csv["SSMD_y"].tolist()))
		print(gene)
		assays = (len(csv.index))
		M = assays
		xlabels = (csv["Assay Type"].tolist())

# THIS IS THE POINT WHERE WE STOP ADDING
	s = np.array(s)
	#print("TEST3")
	ns = np.array(ns)
	#print("TEST4")
	#print(c)
	c = np.array(c)
	#print("TEST5")
	#print(c)
	#c = [arr.flatten() for arr in c]
	#c = np.concatenate(c)
	#print("TEST6")
	fc = np.concatenate(c).flatten()
	#print(fc)
	#ac = np.reshape(ac, (1, len(c)))
	#print("TEST6")
	#print(ac)
	#zeros_array = np.reshape(zeros_array, (1, len(zeros_input)))
	#print("TEST5")
	nc = np.array(nc)
	fnc = np.concatenate(nc).flatten()
	#print("TEST6")
	x, y = np.meshgrid(np.arange(M), np.arange(N))
	#print("TEST7")
	fig, ax = plt.subplots()
	R = s/1/2
	fs = np.concatenate(s).flatten()
	fns = np.concatenate(ns).flatten()
	R = fs/1/2
	R2 = fns/1/2
#R = s/s.max()/2
	#print("TEST7")
	#print(R)
	#print("TEST8")
	#print(list(R.flat))
	#print("TEST9")
	#print(Rtest)
	#print("TEST10")
	#print(list(Rtest.flat))
	circles = [plt.Circle((j-0.2,i-0.2), radius=r*0.8) for r, j, i in zip(R.flat, x.flat, y.flat)]
	#print("TEST8")
	circles2 = [plt.Circle((j+0.2,i+0.2), radius=r*0.8) for r, j, i in zip(R2.flat, x.flat, y.flat)]
	#print(c)
	#print(nc)
	#print("TEST8")
	#print(circles)
#print(c.flatten()) # this is the matrix colors
	col = PatchCollection(circles, array=fc, cmap="PuOr")
	#col = PatchCollection(circles, array=c.flatten(), cmap="PuOr")
	#col2 = PatchCollection(circles2, array=nc.flatten(), cmap="PuOr")
	col2 = PatchCollection(circles2, array=fnc, cmap="PuOr")
	col.set_clim([-1.5,1.5])
	col2.set_clim([-1.5,1.5])
	ax.add_collection(col)
	ax.add_collection(col2)

	#print(xlabels)
	#print(ylabels)
	#print(np.arange(M))
	#print(np.arange(N))
	ax.set(xticks=np.arange(M), yticks=np.arange(N),
	xticklabels=xlabels, yticklabels=ylabels)
	ax.set_xticks(np.arange(M+1)-0.5, minor=True)
	ax.set_yticks(np.arange(N+1)-0.5, minor=True)
	ax.grid(which='minor',linewidth=0.2)
	ax.set_xticklabels(labels=xlabels,rotation = 90,fontsize=5)

	fig.colorbar(col)
	if "base" in prefix:
		fig.set_size_inches(6,3)
	elif "stim" in prefix:
		fig.set_size_inches(12,3)
	plt.tight_layout()
	#plt.subplots_adjust(bottom=0.4)
	#plt.show()
	plt.savefig("heatmap_" + prefix + ".pdf",dpi=300)
	plt.close()
