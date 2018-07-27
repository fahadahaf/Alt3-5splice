from SpliceGrapher.shared.utils        import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.predict.SpliceSite  import* 
from SpliceGrapher.formats.FastaLoader import FastaLoader
from SpliceGrapher.formats.fasta import *

from SpliceGrapher.formats.GeneModel   import *
#from SpliceGrapher.formats.loadGTF     import *
from SpliceGrapher.formats.fasta       import FastaRecord
from SpliceGrapher.SpliceGraph         import *
from SpliceGrapher.formats.GTFLoader   import *
from SpliceGrapher.formats.loader import loadGeneModels

from glob     import glob
from optparse import OptionParser
import os, sys, warnings, numpy
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#-----------------------------------------------------------------------#

handle = open('a_thaliana.fa','rU')
records = list(SeqIO.parse(handle,'fasta'))
handle.close()

def unify(flist,iden):
	for i in range(0,len(flist)):
		#for j in range(0,len(flist)):
			#if flist[j][0]==flist[i][0] or abs(flist[j][0]-flist[i][0])<=6:
				#flist[j][1]=flist[i][1] if flist[i][1]>flist[j][1] else flist[j][1]
				#flist[j][0]=flist[i][0] if flist[i][0]<flist[j][0] else flist[j][0]
			#elif flist[j][1]==flist[i][1] or abs(flist[j][1]-flist[i][1])<=6:
				#flist[j][0]=flist[i][0] if flist[i][0]<flist[j][0] else flist[j][0]
				#flist[j][1]=flist[i][1] if flist[i][1]>flist[j][1] else flist[j][1]
		
		if iden=='e':
			for j in range(0,100):
				if (i+j) >= len(flist):
					break		
				if flist[i+j][0]==flist[i][0] or abs(flist[i+j][0]-flist[i][0])<=6:
					if flist[i+j][1]<flist[i][1]: #>
						flist[i][1]=flist[i+j][1]
					else:
						flist[i+j][1]=flist[i][1]
				
					if flist[i][0]>flist[i+j][0]: #<
						flist[i+j][0]=flist[i][0]
					else:
						flist[i][0]=flist[i+j][0]
				
			
				elif flist[i+j][1]==flist[i][1] or abs(flist[i+j][1]-flist[i][1])<=6:
					if flist[i+j][1]<flist[i][1]: #>
						flist[i][1]=flist[i+j][1]
					else:
						flist[i+j][1]=flist[i][1]
				
					if flist[i][0]>flist[i+j][0]: #<
						flist[i+j][0]=flist[i][0]
					else:
						flist[i][0]=flist[i+j][0]
		
		elif iden=='s':
			for j in range(0,100):
				if (i+j) >= len(flist):
					break		
				if flist[i+j][0]==flist[i][0] or abs(flist[i+j][0]-flist[i][0])<=6:
					if flist[i+j][0]<flist[i][0]: #>
						flist[i][0]=flist[i+j][0]
						flist[i][1]=flist[i+j][1]
					else:
						flist[i+j][0]=flist[i][0]
						flist[i+j][1]=flist[i][1]
			
				elif flist[i+j][1]==flist[i][1] or abs(flist[i+j][1]-flist[i][1])<=6:
					if flist[i+j][1]<flist[i][1]: #>
						flist[i][1]=flist[i+j][1]
						flist[i][0]=flist[i+j][0]
					else:
						flist[i+j][1]=flist[i][1]
						flist[i+j][0]=flist[i][0]
		else:
			for j in range(0,100):
				if (i+j) >= len(flist):
					break		
				if flist[i+j][0]==flist[i][0] or abs(flist[i+j][0]-flist[i][0])<=6:
					if flist[i+j][1]>flist[i][1]: #<
						flist[i][1]=flist[i+j][1]
					else:
						flist[i+j][1]=flist[i][1]
				
					if flist[i][0]<flist[i+j][0]: #>
						flist[i+j][0]=flist[i][0]
					else:
						flist[i][0]=flist[i+j][0]
				
			
				elif flist[i+j][1]==flist[i][1] or abs(flist[i+j][1]-flist[i][1])<=6:
					if flist[i+j][1]>flist[i][1]: #<
						flist[i][1]=flist[i+j][1]
					else:
						flist[i+j][1]=flist[i][1]
				
					if flist[i][0]<flist[i+j][0]: #>
						flist[i+j][0]=flist[i][0]
					else:
						flist[i][0]=flist[i+j][0]
					
	print i,j
	return flist

def remove_duplicates(clist):
	b = set(tuple(x) for x in clist)
	glist = [list(x) for x in b]
	
	return glist
	
def adjust_lists(exlist,ident):
	#print exlist
	#if ident =='s':
		#return remove_duplicates(exlist)
	#else:
		#return remove_duplicates(unify(exlist,ident))
	return remove_duplicates(unify(exlist,ident))
	
def write_sequences(listfile,Dir, filename):
	if not os.path.exists(Dir): os.makedirs(Dir)
	
	f = open(Dir+'/'+filename,'wb')
	for entry in listfile:
		chrNo = entry[-1]-1
		start = entry[0]-1
		end = entry[1]-1
		chromo_seq = records[chrNo]
		seq=chromo_seq.seq[start:end]
		f.write(seq.tostring()+'\n')
	f.close()
		
def remove_false_positives(list_true, list_false):
	count = 0
	for i in range(0,len(list_true)):
		if (list_true[i] in list_false):
			count +=1
			list_false.remove(list_true[i])
	
	print count
	return list_false
		

def get_balanced_examples(list_Exon,list_Intron,list_Splice):
	list_final_Exon=[]
	list_final_Intron=[]
	list_final_Splice=[]
	for i in range(0,len(list_Exon)):
		flag=False
		if i <= 10000:
			range_j = 0
		else:
			range_j = i-10000
		for j in range(range_j,i+10000):
			if j>=len(list_Intron):
				break
			if abs(list_Exon[i][0]-list_Intron[j][1])==1 or abs(list_Exon[i][1]-list_Intron[j][0])==1:
				if i <= 40000:
					range_k = 0
				else:
					range_k = i-40000
				for k in range(range_k,i+40000):
					if k >= len(list_Splice):
						break
					if abs(list_Exon[i][0]-list_Splice[k][0])<=(20+6) or abs(list_Exon[i][1]-list_Splice[k][0])<=(20+6):
						list_final_Exon.append(list_Exon[i])
						list_final_Intron.append(list_Intron[j])
						list_final_Splice.append(list_Splice[k])
						flag=True
						break
			if flag:
				break
	return (list_final_Exon,list_final_Intron,list_final_Splice)
				
dirc = str(sys.argv[1])		
print dirc
gene_files = glob(dirc+'/*.gff')

List_Alt5_Exon=[]
List_Alt5_Intron=[]
List_Alt5_Splice=[]

List_Alt3_Exon=[]
List_Alt3_Intron=[]
List_Alt3_Splice=[]

List_restC_Exon=[]
List_restC_Intron=[]
List_restC_Splice=[]

List_restP_Exon=[]
List_restP_Intron=[]
List_restP_Splice=[]

gene_files.sort()
count=0
for gfile in gene_files:
	graph = getFirstGraph(gfile)
	if graph.hasAS():
		nodes = graph.resolvedNodes()
		nodes.sort()
		#if len(nodes) >100:
			#print len(nodes)
		for i in range(0,len(nodes)):
			if nodes[i].isAltDonor():
				chromosome = int(nodes[i].chromosome.lower().strip('chr'))
				coordinates_exon = [nodes[i].minpos,nodes[i].maxpos,chromosome]#
				thisnode = nodes[i]
				thisnodechild = thisnode.children[0]
				coordinates_intron = [nodes[i].donorEnd(),thisnodechild.acceptorEnd()]
				coordinates_intron.sort()
				coordinates_intron = [coordinates_intron[0]+1, coordinates_intron[1]-1,chromosome]#
				coordinates_splice = [nodes[i].donorEnd()-20, nodes[i].donorEnd()+20,chromosome]#
				List_Alt5_Exon.append(coordinates_exon)
				List_Alt5_Intron.append(coordinates_intron)
				List_Alt5_Splice.append(coordinates_splice)
				#final_coordinates_alt5 = [coordinates_exon, coordinates_intron, coordinates_splice]
				#List_Alt5.append(final_coordinates_alt5)
				#print 'alt 5', nodes[i].id, final_coordinates_alt5
			
			elif nodes[i].isAltAcceptor():
				coordinates_exon = [nodes[i].minpos,nodes[i].maxpos,chromosome]#
				chromosome = int(nodes[i].chromosome.lower().strip('chr'))
				thisnode = nodes[i]
				thisnodeparent = thisnode.parents[0]
				coordinates_intron = [thisnodeparent.donorEnd(),thisnode.acceptorEnd()]
				coordinates_intron.sort()
				coordinates_intron = [coordinates_intron[0]+1, coordinates_intron[1]-1,chromosome]#
				coordinates_splice = [nodes[i].acceptorEnd()-20, nodes[i].acceptorEnd()+20,chromosome]#
				List_Alt3_Exon.append(coordinates_exon)
				List_Alt3_Intron.append(coordinates_intron)
				List_Alt3_Splice.append(coordinates_splice)
				#final_coordinates_alt3 = [coordinates_exon, coordinates_intron, coordinates_splice]
				#List_Alt3.append(final_coordinates_alt3)
				#print 'alt 3', nodes[i].id, final_coordinates_alt3
			
			else:
				children = nodes[i].children
				parents = nodes[i].parents
				if len(children)!=0:
					chromosome = int(nodes[i].chromosome.lower().strip('chr'))
					coordinates_exon = [nodes[i].minpos,nodes[i].maxpos,chromosome]#
					thisnode = nodes[i]
					thisnodechild = thisnode.children[0]
					coordinates_intron = [nodes[i].donorEnd(),thisnodechild.acceptorEnd()]
					coordinates_intron.sort()
					coordinates_intron = [coordinates_intron[0]+1, coordinates_intron[1]-1,chromosome]#
					coordinates_splice = [nodes[i].donorEnd()-20, nodes[i].donorEnd()+20,chromosome]#
					List_restC_Exon.append(coordinates_exon)
					List_restC_Intron.append(coordinates_intron)
					List_restC_Splice.append(coordinates_splice)
					#final_coordinates_restC = [coordinates_exon, coordinates_intron, coordinates_splice]
					#List_restC.append(final_coordinates_restC)
					#print 'restC', nodes[i].id, final_coordinates_restC
				if len(parents)!=0:
					chromosome = int(nodes[i].chromosome.lower().strip('chr'))
					coordinates_exon = [nodes[i].minpos,nodes[i].maxpos,chromosome]#
					thisnode = nodes[i]
					thisnodeparent = thisnode.parents[0]
					coordinates_intron = [thisnodeparent.donorEnd(),thisnode.acceptorEnd()]
					coordinates_intron.sort()
					coordinates_intron = [coordinates_intron[0]+1, coordinates_intron[1]-1,chromosome]#
					coordinates_splice = [nodes[i].acceptorEnd()-20, nodes[i].acceptorEnd()+20,chromosome]#
					List_restP_Exon.append(coordinates_exon)
					List_restP_Intron.append(coordinates_intron)
					List_restP_Splice.append(coordinates_splice)
					#final_coordinates_restP = [coordinates_exon, coordinates_intron, coordinates_splice]
					#List_restP.append(final_coordinates_restP)
					#print 'restP', nodes[i].id, final_coordinates_restP
	
	
	
	count+=1
	if count%100==0: 
		print count
	#if count==4000:
		#break			

#print List_Alt5_Exon
#print List_Alt5_Intron
#print List_Alt5_Splice


Adj_List_Alt5_Exon=adjust_lists(List_Alt5_Exon,'e') 
Adj_List_Alt5_Intron=adjust_lists(List_Alt5_Intron,'i')
Adj_List_Alt5_Splice=adjust_lists(List_Alt5_Splice,'s') #was 'e' before for all four (below)

Adj_List_Alt3_Exon=adjust_lists(List_Alt3_Exon,'e')
Adj_List_Alt3_Intron=adjust_lists(List_Alt3_Intron,'i')
Adj_List_Alt3_Splice=adjust_lists(List_Alt3_Splice,'s')

Adj_List_restC_Exon=adjust_lists(List_restC_Exon,'e')
Adj_List_restC_Intron=adjust_lists(List_restC_Intron,'i')
Adj_List_restC_Splice=adjust_lists(List_restC_Splice,'s')

Adj_List_restP_Exon=adjust_lists(List_restP_Exon,'e')
Adj_List_restP_Intron=adjust_lists(List_restP_Intron,'i')
Adj_List_restP_Splice=adjust_lists(List_restP_Splice,'s')

List_rest_Exon = Adj_List_restP_Exon + Adj_List_restC_Exon
List_rest_Intron = Adj_List_restP_Intron + Adj_List_restC_Intron
List_rest_Splice = Adj_List_restP_Splice + Adj_List_restC_Splice
			
List_rest_Exon = remove_duplicates(List_rest_Exon)
List_rest_Intron = remove_duplicates(List_rest_Intron)
List_rest_Splice = remove_duplicates(List_rest_Splice)		

print len(Adj_List_Alt3_Exon),len(Adj_List_Alt3_Intron),len(Adj_List_Alt3_Splice)
print len(Adj_List_Alt5_Exon),len(Adj_List_Alt5_Intron),len(Adj_List_Alt5_Splice)
print len(Adj_List_restP_Exon),len(Adj_List_restP_Intron),len(Adj_List_restP_Splice)
print len(Adj_List_restC_Exon),len(Adj_List_restC_Intron),len(Adj_List_restC_Splice)
print len(List_rest_Exon),len(List_rest_Intron),len(List_rest_Splice) 

List_rest_Exon = remove_false_positives(Adj_List_Alt3_Exon,List_rest_Exon)
List_rest_Exon = remove_false_positives(Adj_List_Alt5_Exon,List_rest_Exon)

List_rest_Intron = remove_false_positives(Adj_List_Alt3_Intron,List_rest_Intron)
List_rest_Intron = remove_false_positives(Adj_List_Alt5_Intron,List_rest_Intron)

List_rest_Splice = remove_false_positives(Adj_List_Alt3_Splice,List_rest_Splice)
List_rest_Splice = remove_false_positives(Adj_List_Alt5_Splice,List_rest_Splice)



list.sort(Adj_List_Alt5_Exon)
list.sort(Adj_List_Alt5_Intron)
list.sort(Adj_List_Alt5_Splice)

list.sort(Adj_List_Alt3_Exon)
list.sort(Adj_List_Alt3_Intron)
list.sort(Adj_List_Alt3_Splice)

list.sort(List_rest_Exon)
list.sort(List_rest_Intron)
list.sort(List_rest_Splice)

###################New Code#######################

Adj_List_Alt5_Exon,Adj_List_Alt5_Intron,Adj_List_Alt5_Splice=get_balanced_examples(Adj_List_Alt5_Exon,Adj_List_Alt5_Intron,Adj_List_Alt5_Splice)
Adj_List_Alt3_Exon,Adj_List_Alt3_Intron,Adj_List_Alt3_Splice=get_balanced_examples(Adj_List_Alt3_Exon,Adj_List_Alt3_Intron,Adj_List_Alt3_Splice)
List_rest_Exon,List_rest_Intron,List_rest_Splice=get_balanced_examples(List_rest_Exon,List_rest_Intron,List_rest_Splice)

write_sequences(Adj_List_Alt5_Exon,'Alt5','Adj_List_Alt5_Exon3.txt')
write_sequences(Adj_List_Alt5_Intron,'Alt5','Adj_List_Alt5_Intron3.txt')
write_sequences(Adj_List_Alt5_Splice,'Alt5','Adj_List_Alt5_Splice3.txt')

write_sequences(Adj_List_Alt3_Exon,'Alt3','Adj_List_Alt3_Exon3.txt')
write_sequences(Adj_List_Alt3_Intron,'Alt3','Adj_List_Alt3_Intron3.txt')
write_sequences(Adj_List_Alt3_Splice,'Alt3','Adj_List_Alt3_Splice3.txt')

#write_sequences(Adj_List_restC_Exon,'Adj_List_restC_Exon3.txt')
#write_sequences(Adj_List_restC_Intron,'Adj_List_restC_Intron3.txt')
#write_sequences(Adj_List_restC_Splice,'Adj_List_restC_Splice3.txt')

#write_sequences(Adj_List_restP_Exon,'Adj_List_restP_Exon3.txt')
#write_sequences(Adj_List_restP_Intron,'Adj_List_restP_Intron3.txt')
#write_sequences(Adj_List_restP_Splice,'Adj_List_restP_Splice3.txt')

write_sequences(List_rest_Exon,'Rest','List_rest_Exon3.txt')
write_sequences(List_rest_Intron,'Rest','List_rest_Intron3.txt')
write_sequences(List_rest_Splice,'Rest','List_rest_Splice3.txt')
