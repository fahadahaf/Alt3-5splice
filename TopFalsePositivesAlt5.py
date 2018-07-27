from PyML import *
import numpy as np
import pickle
import os

data_splice = np.loadtxt('Alt5_Splice3.txt',dtype=str)

with open('Alt5_Res_Aggregate_Final_Pattern.pickle') as f:
#with open('Alt5_Res_Individual_Final_Pattern.pickle') as f:
	res = pickle.load(f)
	
res = res[0]
#L_Pos_Main=[]
#L_Neg_Main=[]
l_pos = []


for i in range(0,5):
	ids = np.asarray(res.getPatternID(i))
	pred = np.asarray(res.getPredictedLabels(i))
	given = np.asarray(res.getGivenLabels(i))
	dfunc = res.getDecisionFunction(i)
	dfunc2 = np.argsort(dfunc)

	count = 0

	for j in range(len(dfunc2)-1,-1, -1):
		pos = pred[dfunc2[j]]
		pos2 = given[dfunc2[j]]
		if pos=='1' and pos2=='0' and count<10:
			#print pos,pos2
			l_pos.append(int(ids[dfunc2[j]]))
			#print dfunc2
			count = count+1
			
l_pos = list(set(l_pos))


exPosSplice = data_splice[l_pos]



dirc='Alt5_Top_False_Positives'
os.makedirs(dirc)
np.savetxt(dirc+'/'+'examples_Alt5_splice_false_pos.txt',exPosSplice,fmt='%s')


