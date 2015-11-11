#!/usr/bin/python

import re
import numpy as np
import pandas as pd
from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist


plt.style.use('ggplot')

class CodonBook(object):
	def __init__(self):
		self.codontable={'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu','CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
						 'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met','GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
						 'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser','CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
						 'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
						 'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop','CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
						 'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
						 'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp','CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
						 'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}
		self.sensitivecodons=['GCC','CGA','CGT','CGC','CAA','GAA','GAG','GGC','GGT','ATT','ATC','CTA','CTC','CTT','CCA',
							  'TCC','ACC','GTA','GTC','GTG']
		self.regulatorycodons=['ATT','ATC','CTA','ACC','GTC','GTG']
		self.codonusage=pd.DataFrame([],index=self.codontable.keys())
		self.codonusage['amino_acid']=self.codonusage.index.map(lambda x: self.codontable[x])
		self.aalist=self.codonusage['amino_acid'].unique() # amino acid list

	def add_gene(self,seqrecord,namestartwith='gene='):  # fasta seqrecord object from Biopython
		# parse seqrecord object
		genename=re.search(namestartwith+'[a-z,0-9,\ ,\',\-,\.]+',seqrecord.description,flags=re.IGNORECASE)
		if genename > 0:
			genename=genename.group(0)[len(namestartwith):]
		else:
			genename='noname'	
		seq=str(seqrecord.seq)

		# extract codon info
		triplets=pd.Series([seq[i:i+3] for i in range(3,len(seq)-3,3)])   # split sequences into codons (exclude start and stop codon)
		triplet_counts=triplets.value_counts().to_dict()
		df=pd.DataFrame(triplet_counts.values(),index=triplet_counts.keys(),columns=[genename])

		# add to codonusage dataframe
		self.codonusage=pd.merge(self.codonusage,df,how='outer',left_index=True,right_index=True).fillna(0)
		self.codonusage=self.codonusage.sort_values('amino_acid')

	def normalize_by_aa(self):  # normalize codon frequencies within each amino acid
		df=self.codonusage.copy()
		
		for name, group in df.groupby('amino_acid'):
			subdf=group.drop('amino_acid',axis='columns')
			df.ix[subdf.index,1:]=subdf.apply(lambda x: x/sum(x))

		return df.fillna(0)

	def normalize_by_totalaa(self):
		df=self.codonusage.copy()
		df.ix[:,1:]/df.ix[:,1:].sum()
		
		return df.fillna(0)


	# use this after creating codonusage by add_gene
	def expressed_codon(self,expressiondf):  # expressiondf is a pd.DataFrame and contains gene names (1st column) and counts (2nd column)
		df=expressiondf
		self.expressed_codonusage=self.codonusage.copy()
		
		for gene in df.index:


class SynGeneName(object):
	def __init__(self,gbrecord): # gbrecord is SeqIO.seq object: gbrecord=SeqIO.parse(open('filepath'),'genbank').next()
		self.namedict={}
		for f in gbrecord.features:
			if f.type=='CDS':
				genename=f.qualifiers['gene'][0]
				synname=re.split(';\ *',f.qualifiers['gene_synonym'][0])
				if not self.namedict.has_key(genename):
					synname.append(genename)
					self.namedict.update({name:synname for name in synname})
				

	def lookup(self,primaryname):
		if self.namedict.has_key(primaryname):
			return self.namedict[primaryname]
		else:
			print(primaryname+' does not exist')





class tAI(object):
	def __init__(self):
		self.selective_constraint=pd.DataFrame(index=['I','G','U','C','L'],columns=['T','C','A','G']).fillna(0)
		# initialize selective constraint for wobbling using values from dos Reis et al (2004). 
		self.selective_constraint.ix['G','T']=0.561
		self.selective_constraint.ix['I','C']=0.28
		self.selective_constraint.ix['I','A']=0.9999
		self.selective_constraint.ix['U','G']=0.68
		self.selective_constraint.ix['L','A']=0.89

		self.codondecode_matrix=pd.DataFrame()
		for firstnt in ['T','C','A','G']:
			for secondnt in ['T','C','A','G']:
				for thridnt in ['T','C','A','G']:
					codon=firstnt+secondnt+thridnt
					anticodon=str(Seq(codon).reverse_complement().transcribe())  # anticodons use 'U' here
					self.codondecode_matrix.ix[codon,anticodon]=1

					# wobbling 
					if thridnt=='T':
						wobble_anticodon='G'+anticodon[1:]
						self.codondecode_matrix.ix[codon,wobble_anticodon]=1-self.selective_constraint.ix['G','T']
					elif thridnt=='C':
						wobble_anticodon='A'+anticodon[1:]  # 'I' on tRNA is derived from 'A'
						self.codondecode_matrix.ix[codon,wobble_anticodon]=1-self.selective_constraint.ix['I','C']
					elif thridnt=='A':
						wobble_anticodon='A'+anticodon[1:]  # 'I' on tRNA is derived from 'A'
						self.codondecode_matrix.ix[codon,wobble_anticodon]=1-self.selective_constraint.ix['I','A']
					elif thridnt=='G':
						wobble_anticodon='U'+anticodon[1:] 
						self.codondecode_matrix.ix[codon,wobble_anticodon]=1-self.selective_constraint.ix['U','G']
					
		# special wobbling case: 'ATA' (Ile) decoded by 'LAU' ('C' is modified to become 'L')
		self.codondecode_matrix.ix['ATA','LAU']=1-self.selective_constraint.ix['L','A']
		self.codondecode_matrix=self.codondecode_matrix.fillna(0)

		# deal with stop codons
		self.codondecode_matrix.ix['TGA',:]=0
		self.codondecode_matrix.ix['TAA',:]=0
		self.codondecode_matrix.ix['TAG',:]=0



	def weights_from_tRNAgenecopy(self,gbrecord):  # gbrecord is SeqIO.seq object: gbrecord=SeqIO.parse(open('filepath'),'genbank').next()
		self.tRNAtable=pd.DataFrame(index=self.codondecode_matrix.columns,columns=['amino_acid','count']).fillna(0)
		for f in gbrecord.features:
			if f.type=='tRNA' and (f.qualifiers['product'][0]!='tRNA-OTHER' and f.qualifiers['product'][0]!='tRNA-Sec'): 
				anticodon=f.qualifiers['note'][0][11:14]
				aa_name=f.qualifiers['product'][0][-3:]
				
				if aa_name=='Ile' and anticodon=='CAU':
					self.tRNAtable.ix['LAU','amino_acid']=aa_name
					self.tRNAtable.ix['LAU','count']+=1
				else:
					self.tRNAtable.ix[anticodon,'amino_acid']=aa_name
					self.tRNAtable.ix[anticodon,'count']+=1

		w=self.codondecode_matrix.dot(self.tRNAtable['count'])
		w=w/w.max()

		return w

	def weights_from_tRNAabundance(self,trnadf): # trnadf should be pd.DataFrame and include anticodon (1st column) and abundances (2nd column) 
		self.tRNAtable=pd.DataFrame(index=self.codondecode_matrix.columns,columns=['amino_acid']).fillna(0)
		self.tRNAtable=pd.merge(self.tRNAtable,trnadf,how='outer',left_index=True,right_index=True).reindex(self.codondecode_matrix.columns)

		w=self.codondecode_matrix.dot(self.tRNAtable['ave'].fillna(0))
		w=w/w.max()

		return w
		
	
	def lookup_wobbleanticodon(self,codon):
		return self.codondecode_matrix.ix[codon,:][self.codondecode_matrix.ix[codon,:]>0]







	


if __name__ == "__main__":

	# load fasta
	filepath='sequence_data/lambda_cds.txt'
	record=SeqIO.parse(open(filepath,'rU'),'fasta')
	book=CodonBook()

	for feature in record:
		book.add_gene(feature)

	# df=book.normalize_by_totalaa().ix[:,1:].transpose()
	# df=book.codonusage.ix[book.sensitivecodons,1:].transpose()
	# Compute and plot first dendrogram.
	fig = plt.figure(figsize=(18,12))
	ax1 = fig.add_axes([0.05,0.1,0.2,0.8])
	y=linkage(df,method='ward')
	z=dendrogram(y,orientation='right')
	ax1.set_xticks([])
	ax1.set_yticks([])


	# Plot distance matrix.
	axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
	idx1 = z['leaves']
	# idx2 = Z2['leaves']
	df = df.ix[idx1,:]
	# D = D.ix[:,idx2]
	im = axmatrix.matshow(df, aspect='auto', origin='lower', cmap=matplotlib.cm.YlOrRd)
	axmatrix.set_xticks(range(len(df.columns)))
	axmatrix.set_yticks(range(len(df.index)))
	axmatrix.set_xticklabels(df.columns,rotation=90)
	axmatrix.set_yticklabels(df.index)

	# Plot colorbar.
	axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
	plt.colorbar(im,cax=axcolor)
	fig.savefig('cluster.pdf')
	# fig.show()








