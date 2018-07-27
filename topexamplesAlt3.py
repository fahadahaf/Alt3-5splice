from PyML import *
import numpy as np
import pickle
import os

data_exon = np.loadtxt('Alt3_Exon3.txt',dtype=str)
data_intron = np.loadtxt('Alt3_Intron3.txt',dtype=str)
data_splice = np.loadtxt('Alt3_Splice3.txt',dtype=str)

with open('Alt3_Res_Aggregate_Final_Pattern.pickle') as f:
	res = pickle.load(f)
	
res = res[0]
#L_Pos_Main=[]
#L_Neg_Main=[]
l_pos = []
l_neg = []

for i in range(0,5):
	ids = np.asarray(res.getPatternID(i))
	pred = np.asarray(res.getPredictedLabels(i))
	given = np.asarray(res.getGivenLabels(i))
	dfunc = res.getDecisionFunction(i)
	largmax = list(np.argsort(dfunc)[-10:])
	largmax = largmax[::-1]
	largmin = list(np.argsort(dfunc)[:10])
	top_pred_pos = pred[largmax]
	top_given_pos = given[largmax]
	top_pred_neg = pred[largmin]
	top_given_neg = given[largmin]
	#l_pos = []
	#l_neg = []
	for j in range(0,len(top_pred_neg)):
		if top_pred_pos[j]==top_given_pos[j]:
			l_pos.append(int(ids[largmax[j]]))
		if top_pred_neg[j]==top_given_neg[j]:
			l_neg.append(int(ids[largmin[j]]))
	#L_Pos_Main.append(l_pos)
	#L_Neg_Main.append(l_neg)
			
l_pos = list(set(l_pos))
l_neg = list(set(l_neg))

exPosSplice = data_splice[l_pos]
exPosIntron = data_intron[l_pos]
exPosExon = data_exon[l_pos]

exNegSplice = data_splice[l_neg]
exNegIntron = data_intron[l_neg]
exNegExon = data_exon[l_neg]



dirc='Alt3_Top_Positives'
os.makedirs(dirc)
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_splice_pos.txt',exPosSplice,fmt='%s')
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_intron_pos.txt',exPosIntron,fmt='%s')
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_exon_pos.txt',exPosExon,fmt='%s')

dirc='Alt3_Top_Negatives'
os.makedirs(dirc)
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_splice_neg.txt',exNegSplice,fmt='%s')
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_intron_neg.txt',exNegIntron,fmt='%s')
np.savetxt(dirc+'/'+'examples_Alt3_aggregate_exon_neg.txt',exNegExon,fmt='%s')

