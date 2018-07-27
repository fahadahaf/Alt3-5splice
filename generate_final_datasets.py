import glob
import numpy as np


def generate_trim_dataset(Type):
	files_type = glob.glob(Type+'/*.txt')
	#files_rest = glob.glob('Rest/*.txt')
	
	for i in range (0,len(files_type)):
		filename_type = files_type[i]
		filename_rest = 'Rest/List_rest_'+ filename_type.split('_')[-1]
		#print filename_type
		#print filename_rest
		f = open(filename_type,'rb')
		pos = f.readlines()
		for i in range(0,len(pos)):
			pos[i]=pos[i].strip('\n')

		pos=np.asarray(pos)
		
		f = open(filename_rest,'rb')
		neg = f.readlines()
		for i in range(0,len(neg)):
			neg[i]=neg[i].strip('\n')
		
		neg=np.asarray(neg)
		
		l1=np.ones(len(pos),dtype=int)

		l0=np.zeros(len(neg),dtype=int)

		labels = np.concatenate((l1,l0),axis=0)

		final = np.concatenate((pos,neg),axis=0)
		
		indices = [i for i in range(10000,30000)]

		final = np.delete(final,indices)

		labels = np.delete(labels,indices)
		
		
		np.savetxt('Labels_'+Type+'_'+filename_type.split('_')[-1],labels,fmt='%s')
		np.savetxt(Type+'_'+filename_type.split('_')[-1],final,fmt='%s')
		

#Alt5
generate_trim_dataset('Alt5')

#Alt3
generate_trim_dataset('Alt3')

