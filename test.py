#!/usr/bin/python

import re
import numpy as np
import pandas as pd
from pybedtools import BedTool
from Bio import SeqIO
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist


plt.style.use('ggplot')

class CodonUsageBook(object):
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

	def add_gene(self,df):  # df contains codons as index and their counts on a column
		self.codonusage=pd.merge(self.codonusage,df,how='outer',left_index=True,right_index=True).fillna(0)
		self.codonusage=self.codonusage.sort_values('amino_acid')

	def normalize_by_aa(self):  # normalize codon frequencies within each amino acid
		df=self.codonusage.copy()
		
		for name, group in df.groupby('amino_acid'):
			subdf=group.drop('amino_acid',axis='columns')
			df.ix[subdf.index,1:]=subdf.apply(lambda x: x/sum(x))

		return df.fillna(0)

	def normalize_by_totalaa(self):
		return self.codonusage.ix[:,1:]/self.codonusage.ix[:,1:].sum()


class CodonUsagePerGene(object):
	def __init__(self,seqrecord): # fasta seqrecord object
		self.seqrecord=seqrecord
		self.seq=seqrecord.seq.tostring()
		self.description=re.split(' \[|\] \[|\]',self.seqrecord.description)
		self.name=self.description[1]

	def codon_usage(self):
		triplets=pd.Series([self.seq[i:i+3] for i in range(3,len(self.seq),3)])   # split sequences into codons (exclude start codon)
		triplet_counts=triplets.value_counts().to_dict()
		return pd.DataFrame(triplet_counts.values(),index=triplet_counts.keys(),columns=[self.name])


 



if __name__ == "__main__":

# load fasta
	filepath='../sequence_data/BW25113_cds.txt'
	record=SeqIO.parse(open(filepath,'rU'),'fasta')
	book=CodonUsageBook()

	for feature in record:
		gene=CodonUsagePerGene(feature)
		book.add_gene(gene.codon_usage())

	df=book.normalize_by_totalaa()

	plt.pcolor(df.ix[:,1:],cmap=matplotlib.cm.Blues,vmin=0,vmax=df.ix[:,1:].max().max())



	D=df.ix[book.regulatorycodons,1:].transpose()
	# D=book.codonusage.ix[:,1:].transpose()
	# Compute and plot first dendrogram.
	fig = plt.figure(figsize=(18,12))
	ax1 = fig.add_axes([0.05,0.1,0.2,0.8])
	Y = linkage(D, method='ward')
	Z1 = dendrogram(Y, orientation='right')
	ax1.set_xticks([])
	ax1.set_yticks([])

	# Compute and plot second dendrogram.
	# ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
	# Y = linkage(D.transpose(), method='weighted')
	# Z2 = dendrogram(Y)
	# ax2.set_xticks([])
	# ax2.set_yticks([])

	# Plot distance matrix.
	axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
	idx1 = Z1['leaves']
	# idx2 = Z2['leaves']
	D = D.ix[idx1,:]
	# D = D.ix[:,idx2]
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=matplotlib.cm.YlGnBu)
	axmatrix.set_xticks(range(len(D.columns)))
	axmatrix.set_yticks(range(len(D.index)))
	axmatrix.set_xticklabels(D.columns,rotation=90)
	axmatrix.set_yticklabels(D.index.map(lambda x:x[5:]))

	# Plot colorbar.
	# axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
	# fig.savefig('cluster.pdf')
	fig.show()








