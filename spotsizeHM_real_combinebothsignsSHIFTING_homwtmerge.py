import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import glob
import pandas


def is_single_line_file(file_path):
	with open(file_path, 'r') as file:
		lines = file.readlines()
		return len(lines) == 1

#Frequency of movement [day]    0.024390      0.024390
#Frequency of movement [night]  0.039604      0.039604
#Location in well [day]         0.157143      0.157143
#Location in well [night]       0.310345      0.310345
#Magnitude of movement [day]    0.025974      0.025974
#Magnitude of movement [night]  0.072100      0.072100
#Seizures [day & night]         0.041667      0.041667
#Sleep [day & night]            0.046875      0.046875
#Frequency of movement [day]    0.474248
#Frequency of movement [night]  0.195755
#Location in well [day]         0.458038
#Location in well [night]       0.493548
#Magnitude of movement [day]   -0.342869
#Magnitude of movement [night] -0.385686
#Seizures [day & night]        -0.296972
#Sleep [day & night]            0.426471

#Significance    SSMD_y
#means2HM_ssmd_stimuli_Kiaa0232_Box12_4_21BC_linearmodel_kiaa0232Box12_wt_vs_kiaa0232Box12_hom_11693279_4294967294
#,Assay Type,SSMD_x,Significance,Percents,SSMD_y
prefixlist = []
filelist = []
#posmeansepsign2HM_ssmd_baseline_replicates_atp1a3a-ko_wt_hom.csv
#posmeansepsign2HM_ssmd_stimuli_pvals_atp1a3a-ko_Box2_07_03_20_GK_cosscz_linearmodel_18-atp1a3aBox2_wt_vs_20-atp1a3aBox2_het_7357888_4294967294.csv
#old way above and below, new above
#for mfile0 in glob.glob("means2HM_ssmd_stimuli_Kiaa0232_Box12_4_21BC_linearmodel_kiaa0232Box12_hetandwt_vs_kiaa0232Box12_hom_11693280_4294967294*"):
for mfile0 in glob.glob("*sepsign*2HM_ssmd**csv"):
	if is_single_line_file(mfile0):
		continue
	if "_linearmodel_" in mfile0:
		prefix = mfile0.split("_linearmodel_")[0].split("_ssmd_")[1]
	else:
		prefix = "_".join(mfile0.split(".")[0].split("_ssmd_")[1].split("_")[:-2])
	filelist.append(mfile0)
#	print(prefix)
	prefixlist.append(prefix)
filelist.sort()
#prefixlist.sort()
#print(filelist)
#s = []
#c = []
#ylabels = []
#N = 0
#print("TEST00",prefix)
s = []
ns = []
c = []
nc = []
ylabels = []
N = 0
for mfile in filelist:
#for mfile in glob.glob("*sepsign*2HM_*csv"):
#for mfile in glob.glob("means2HM_ssmd_stimuli_Kiaa0232_Box12_4_21BC_linearmodel_kiaa0232Box12_hetandwt_vs_kiaa0232Box12_hom_11693280_4294967294*"):
	#posmeansepsign2HM_ssmd_stimuli_standard_replicates_sbno1_wt_hom.csv
	#premfile = "_".join(mfile.split("_")[0:-2])
	if "_linearmodel_" in mfile:
		premfile = mfile.split("_linearmodel_")[0].split("_ssmd_")[1]
	else:
		premfile = "_".join(mfile.split(".")[0].split("_ssmd_")[1].split("_")[:-2])
	#print("TESTXX",premfile,prefix)
	##if prefix != premfile:
	#if prefix not in mfile:
	##	continue
	#print("TESTXX",premfile,prefix)
	print("TEST0",mfile)
	#N = N + 1
	csv = pandas.read_csv(mfile)
	csv = csv.loc[:, ~csv.columns.str.contains('^Unnamed')]
	# below commented is old
	#gene = mfile.split("_linearmodel_")[1].split(".")[0].split("_")[0] + "_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[1] + "_vs_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[4]
	if "_linearmodel_" in mfile:
		gene = mfile.split("_linearmodel_")[0].split("_pvals_")[1].split("_")[0] + "_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[1] + "_vs_" + mfile.split("_linearmodel_")[1].split(".")[0].split("_")[4]
	else:
		gene = mfile.split("_replicates_")[1].split(".")[0].split("_")[0] + "_" + mfile.split("_replicates_")[1].split(".")[0].split("_")[1] + "_vs_" + mfile.split("_replicates_")[1].split(".")[0].split("_")[2]
	#ylabels.append(gene)
	if("neg" in mfile):
		ns.append(np.array(csv["Percents"].tolist())) # Significance is the % of the assays with significance, not sure why it ended up being named that)
		nc.append(np.array(csv["SSMD_y"].tolist()))
	#	gene = gene + "_negative"
		ylabels.append(gene)
	if("pos" in mfile):
		N = N + 1
		s.append(np.array(csv["Percents"].tolist())) # Significance is the % of the assays with significance, not sure why it ended up being named that)
		c.append(np.array(csv["SSMD_y"].tolist()))
		#print(c)
		print("TESTTEST",np.array(csv["Percents"].tolist()).shape)
	#	gene = gene + "_positive"
	print(gene)
	#ylabels.append(gene)
	assays = (len(csv.index))
	M = assays
	xlabels = (csv["Assay Type"].tolist())
s = np.array(s)
ns = np.array(ns)
c = np.array(c)
fc = np.concatenate(c).flatten()
nc = np.array(nc)
fnc = np.concatenate(nc).flatten()
x, y = np.meshgrid(np.arange(M), np.arange(N))
fig, ax = plt.subplots()
R = s/1/2
fs = np.concatenate(s).flatten()
fns = np.concatenate(ns).flatten()
R = fs/1/2
R2 = fns/1/2
circles = [plt.Circle((j-0.2,i-0.2), radius=r*0.8) for r, j, i in zip(R.flat, x.flat, y.flat)]
circles2 = [plt.Circle((j+0.2,i+0.2), radius=r*0.8) for r, j, i in zip(R2.flat, x.flat, y.flat)]
col = PatchCollection(circles, array=fc, cmap="PuOr")
col2 = PatchCollection(circles2, array=fnc, cmap="PuOr")
col.set_clim([-1.5,1.5])
col2.set_clim([-1.5,1.5])
ax.add_collection(col)
ax.add_collection(col2)
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
	fig.set_size_inches(5,7)
elif "stim" in prefix:
	fig.set_size_inches(9,7)
plt.tight_layout()
#plt.subplots_adjust(bottom=0.4)
#plt.show()
plt.savefig("heatmap_" + "asd_nov23_wtvshom_standardSTIMULUS" + ".pdf",dpi=300)
plt.close()
