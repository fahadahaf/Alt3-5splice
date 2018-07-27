from PyML import *
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt

data_exon = np.loadtxt('Alt3_Exon3.txt',dtype=str)
data_intron = np.loadtxt('Alt3_Intron3.txt',dtype=str)
data_splice = np.loadtxt('Alt3_Splice3.txt',dtype=str)

Exon =sequenceData.spectrum_data(data_exon,6)
Intron =sequenceData.spectrum_data(data_intron,6)
#Splice =sequenceData.positional_kmer_data(data_splice,6)
Splice = sequenceData.spectrum_data(data_splice,6)

Exon.attachLabels(Labels('Labels_Alt3_Exon3.txt'))
Intron.attachLabels(Labels('Labels_Alt3_Intron3.txt'))
Splice.attachLabels(Labels('Labels_Alt3_Splice3.txt'))

ids1 = [str(i) for i in range(0,len(Exon))]
ids2 = [str(i) for i in range(0,len(Intron))]
ids3 = [str(i) for i in range(0,len(Splice))]

Exon.labels.patternID = ids1
Intron.labels.patternID = ids2
Splice.labels.patternID = ids3

classifier1 = SVM()
classifier2 = SVM(ker.Gaussian(gamma = 0.01))
classifier3 = SVM(ker.Polynomial(degree = 3))

res_E_SVM = classifier1.cv(Exon)
res_E_SVM_G = classifier2.cv(Exon)
res_E_SVM_P = classifier3.cv(Exon)

dirc = 'Alt3_Individual_Spectrum_Results'
os.makedirs(dirc)

res_Exon = [['Kernel', 'Success_Rate', 'Balanced_Success_Rate', 'AUC']]
res_Exon.append(['Linear', str(res_E_SVM.getSuccessRate()),str(res_E_SVM.getBalancedSuccessRate()),str(res_E_SVM.getROC())])
res_Exon.append(['Gaussian', str(res_E_SVM_G.getSuccessRate()),str(res_E_SVM_G.getBalancedSuccessRate()),str(res_E_SVM_G.getROC())])
res_Exon.append(['Polynomial', str(res_E_SVM_P.getSuccessRate()),str(res_E_SVM_P.getBalancedSuccessRate()),str(res_E_SVM_P.getROC())])

res_E_SVM.plotROC(dirc+'/'+'Linear_ROC_Exon.png',show=False)
plt.clf()
res_E_SVM_G.plotROC(dirc+'/'+'Gaussian_ROC_Exon.png',show=False)
plt.clf()
res_E_SVM_P.plotROC(dirc+'/'+'Polynomial_ROC_Exon.png',show=False)
plt.clf()

np.savetxt(dirc+'/'+'results_Exon.txt',res_Exon,fmt='%s')

classifier1 = SVM()
classifier2 = SVM(ker.Gaussian(gamma = 0.01))
classifier3 = SVM(ker.Polynomial(degree = 3))

res_I_SVM = classifier1.cv(Intron)
res_I_SVM_G = classifier2.cv(Intron)
res_I_SVM_P = classifier3.cv(Intron)


res_Intron = [['Kernel', 'Success_Rate', 'Balanced_Success_Rate', 'AUC']]
res_Intron.append(['Linear', str(res_I_SVM.getSuccessRate()),str(res_I_SVM.getBalancedSuccessRate()),str(res_I_SVM.getROC())])
res_Intron.append(['Gaussian', str(res_I_SVM_G.getSuccessRate()),str(res_I_SVM_G.getBalancedSuccessRate()),str(res_I_SVM_G.getROC())])
res_Intron.append(['Polynomial', str(res_I_SVM_P.getSuccessRate()),str(res_I_SVM_P.getBalancedSuccessRate()),str(res_I_SVM_P.getROC())])

res_I_SVM.plotROC(dirc+'/'+'Linear_ROC_Intron.png',show=False)
plt.clf()
res_I_SVM_G.plotROC(dirc+'/'+'Gaussian_ROC_Intron.png',show=False)
plt.clf()
res_I_SVM_P.plotROC(dirc+'/'+'Polynomial_ROC_Intron.png',show=False)
plt.clf()

np.savetxt(dirc+'/'+'results_Intron.txt',res_Intron,fmt='%s')


classifier1 = SVM()
classifier2 = SVM(ker.Gaussian(gamma = 0.01))
classifier3 = SVM(ker.Polynomial(degree = 3))

res_S_SVM = classifier1.cv(Splice)
res_S_SVM_G = classifier2.cv(Splice)
res_S_SVM_P = classifier3.cv(Splice)


res_Splice = [['Kernel', 'Success_Rate', 'Balanced_Success_Rate', 'AUC']]
res_Splice.append(['Linear', str(res_S_SVM.getSuccessRate()),str(res_S_SVM.getBalancedSuccessRate()),str(res_S_SVM.getROC())])
res_Splice.append(['Gaussian', str(res_S_SVM_G.getSuccessRate()),str(res_S_SVM_G.getBalancedSuccessRate()),str(res_S_SVM_G.getROC())])
res_Splice.append(['Polynomial', str(res_S_SVM_P.getSuccessRate()),str(res_S_SVM_P.getBalancedSuccessRate()),str(res_S_SVM_P.getROC())])

res_S_SVM.plotROC(dirc+'/'+'Linear_ROC_Splice.png',show=False)
plt.clf()
res_S_SVM_G.plotROC(dirc+'/'+'Gaussian_ROC_Splice.png',show=False)
plt.clf()
res_S_SVM_P.plotROC(dirc+'/'+'Polynomial_ROC_Splice.png',show=False)
plt.clf()

np.savetxt(dirc+'/'+'results_Splice.txt',res_Splice,fmt='%s')

with open('Alt3_Res_Individual_Final_spectrum_Pattern.pickle','w') as f:
	pickle.dump([res_E_SVM,res_E_SVM_G,res_E_SVM_P,res_I_SVM,res_I_SVM_G,res_I_SVM_P,res_S_SVM,res_S_SVM_G,res_S_SVM_P],f)
