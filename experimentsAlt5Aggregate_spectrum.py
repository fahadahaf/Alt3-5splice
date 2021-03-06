from PyML import *
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt

data_exon = np.loadtxt('Alt5_Exon3.txt',dtype=str)
data_intron = np.loadtxt('Alt5_Intron3.txt',dtype=str)
data_splice = np.loadtxt('Alt5_Splice3.txt',dtype=str)

Exon =sequenceData.spectrum_data(data_exon,6)
Intron =sequenceData.spectrum_data(data_intron,6)
#Splice =sequenceData.positional_kmer_data(data_splice,6)
Splice =sequenceData.spectrum_data(data_splice,6)

Exon.attachLabels(Labels('Labels_Alt5_Exon3.txt'))
Intron.attachLabels(Labels('Labels_Alt5_Intron3.txt'))
Splice.attachLabels(Labels('Labels_Alt5_Splice3.txt'))

ids1 = [str(i) for i in range(0,len(Exon))]
ids2 = [str(i) for i in range(0,len(Intron))]
ids3 = [str(i) for i in range(0,len(Splice))]

Exon.labels.patternID = ids1
Intron.labels.patternID = ids2
Splice.labels.patternID = ids3

datas=[Exon,Intron,Splice]
dataAggregate = Aggregate(datas)

classifier1 = SVM()
classifier2 = SVM(ker.Gaussian(gamma = 0.01))
classifier3 = SVM(ker.Polynomial(degree = 3))

res_E_SVM = classifier1.cv(dataAggregate)
res_E_SVM_G = classifier2.cv(dataAggregate)
res_E_SVM_P = classifier3.cv(dataAggregate)

dirc = 'Alt5_Aggregate_Spectrum_Results'
os.makedirs(dirc)

res_Exon = [['Kernel', 'Success_Rate', 'Balanced_Success_Rate', 'AUC']]
res_Exon.append(['Linear', str(res_E_SVM.getSuccessRate()),str(res_E_SVM.getBalancedSuccessRate()),str(res_E_SVM.getROC())])
res_Exon.append(['Gaussian', str(res_E_SVM_G.getSuccessRate()),str(res_E_SVM_G.getBalancedSuccessRate()),str(res_E_SVM_G.getROC())])
res_Exon.append(['Polynomial', str(res_E_SVM_P.getSuccessRate()),str(res_E_SVM_P.getBalancedSuccessRate()),str(res_E_SVM_P.getROC())])

res_E_SVM.plotROC(dirc+'/'+'Linear_ROC.png',show=False)
plt.clf()
res_E_SVM_G.plotROC(dirc+'/'+'Gaussian_ROC.png',show=False)
plt.clf()
res_E_SVM_P.plotROC(dirc+'/'+'Polynomial_ROC.png',show=False)
plt.clf()

np.savetxt(dirc+'/'+'results.txt',res_Exon,fmt='%s')

with open('Alt5_Res_Aggregate_Final_spectrum_Pattern.pickle','w') as f:
	pickle.dump([res_E_SVM,res_E_SVM_G,res_E_SVM_P],f)
