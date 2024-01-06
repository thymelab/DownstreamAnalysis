import pandas
import numpy as np
import glob

genodict = {}
for sfile in glob.glob("ssmd*"):
	print(sfile)
	csv = pandas.read_csv(sfile)
	csv = csv.loc[:, ~csv.columns.str.contains('^Unnamed')]
	csv = csv.fillna(0)
	csvsig = csv[csv["Significance"] == "Significant"]
	csvpos = csvsig[csvsig["SSMD"] > 0]
	csvneg = csvsig[csvsig["SSMD"] < 0]
	assays_list = list(set(csvsig["Assay Type"].tolist()))
	assays_list0 = list(set(csv["Assay Type"].tolist()))
	nonsigcats = set.difference(set(assays_list0) - set(assays_list))
	means = csvsig.groupby("Assay Type").mean()#agg({'SSMD':custom_mean})#mean(skipna=False)
	posmeans = csvpos.groupby("Assay Type").mean()#agg({'SSMD':custom_mean})#mean(skipna=False)
	negmeans = csvneg.groupby("Assay Type").mean()#agg({'SSMD':custom_mean})#mean(skipna=False)
	means = means.reset_index()
	negmeans = negmeans.reset_index()
	posmeans = posmeans.reset_index()
	sums = csvsig.groupby("Assay Type").count().astype(np.int64)
	possums = csvpos.groupby("Assay Type").count().astype(np.int64)
	negsums = csvneg.groupby("Assay Type").count().astype(np.int64)
	sums = sums.reset_index()
	negsums = negsums.reset_index()
	possums = possums.reset_index()
	sumstotal = csv.groupby("Assay Type").count().astype(np.int64)
	possumstotal = csv.groupby("Assay Type").count().astype(np.int64)
	negsumstotal = csv.groupby("Assay Type").count().astype(np.int64)
	sumstotal = sumstotal.reset_index()
	possumstotal = possumstotal.reset_index()
	negsumstotal = negsumstotal.reset_index()
	
	sums.set_index('Assay Type', inplace=True)
	negsums.set_index('Assay Type', inplace=True)
	possums.set_index('Assay Type', inplace=True)
	sumstotal.set_index('Assay Type', inplace=True)
	possumstotal.set_index('Assay Type', inplace=True)
	negsumstotal.set_index('Assay Type', inplace=True)

	sumstotal['Percents'] = sums['SSMD'] / sumstotal['SSMD']
	possumstotal['Percents'] = possums['SSMD'] / possumstotal['SSMD']
	negsumstotal['Percents'] = negsums['SSMD'] / negsumstotal['SSMD']
	sumstotal = sumstotal.reset_index()
	possumstotal = possumstotal.reset_index()
	negsumstotal = negsumstotal.reset_index()
	percentassays = sumstotal
	pospercentassays = possumstotal
	negpercentassays = negsumstotal
	for ns in nonsigcats:
		percentassays.loc[percentassays["Assay Type"]== ns] = [ns,float(0),float(0),float(0),float(0)]
		pospercentassays.loc[pospercentassays["Assay Type"]== ns] = [ns,float(0),float(0),float(0),float(0)]
		negpercentassays.loc[negpercentassays["Assay Type"]== ns] = [ns,float(0),float(0),float(0),float(0)]
	final = percentassays.merge(means,how='left', left_on="Assay Type", right_on="Assay Type")
	posfinal = pospercentassays.merge(posmeans,how='left', left_on="Assay Type", right_on="Assay Type")
	negfinal = negpercentassays.merge(negmeans,how='left', left_on="Assay Type", right_on="Assay Type")
	final = final.fillna(0)
	posfinal = posfinal.fillna(0)
	negfinal = negfinal.fillna(0)
	posfinal.to_csv("posmeansepsign2HM_" + sfile.split(".")[0] + ".csv")
	negfinal.to_csv("negmeansepsign2HM_" + sfile.split(".")[0] + ".csv")
