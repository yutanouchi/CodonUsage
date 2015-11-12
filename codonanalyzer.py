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
		self.locustag_dict={} # look-up dictionary between primary gene names and locus tag
		self.synnames_dict={}  # look-up dictionary from syn gene names to locus tag

	def add_genes(self,seqrecord):  
	# seqrecord is a genbank SeqIO.seq object from Biopython. Genes are lostes by locus tag (e.g. b0001)
		# parse seqrecord object
		seq=seqrecord.seq

		for f in seqrecord.features:
			if f.type=='CDS': # only look at CDSs
				strand=f.strand   # 1: positive, -1: negative
				start=f.location.start.position
				end=f.location.end.position
				locus_tag=f.qualifiers['locus_tag'][0]
				
				if f.qualifiers.has_key('gene'):
					genename=f.qualifiers['gene'][0] # primary gene name
				else:
					genename=locus_tag

				if f.qualifiers.has_key('gene_synonym'):
					synnames=re.split(';\ *',f.qualifiers['gene_synonym'][0])
				else:
					synnames=[]

				# update dict. the last element of synnames is the primary gene name
				# note that there are multiple copies of ins genes (e.g. 4 copies of insB1)
				# and genename-to-locus_tag dict will be overwritten and leavs only the last entry  
				self.locustag_dict.update({genename.lower():locus_tag})	
				self.locustag_dict.update({locus_tag.lower():genename})
				self.synnames_dict.update({name.lower():locus_tag for name in synnames})

				# fetch cds sequence
				cds_seq=seq[start:end]
				if strand<0: # if cds is on negative strand
					cds_seq=cds_seq.reverse_complement()
				cds_seq=str(cds_seq)
		
				# extract codon info
				triplets=pd.Series([cds_seq[i:i+3] for i in range(3,len(cds_seq)-3,3)])   # split sequences into codons (exclude start and stop codon)
				triplet_counts=triplets.value_counts().to_dict()
				df=pd.DataFrame(triplet_counts.values(),index=triplet_counts.keys(),columns=[locus_tag])

				# add to codonusage dataframe
				self.codonusage=pd.merge(self.codonusage,df,how='outer',left_index=True,right_index=True).fillna(0)
				self.codonusage=self.codonusage.sort_values('amino_acid')

	def lookup_locustag(self,genename):
		if self.locustag_dict.has_key(genename.lower()):
			return self.locustag_dict[genename.lower()]
		else:
			# print(genename+' does not exist in locustag dict')
			return ''

	def lookup_synnames(self,genename):
		if self.synnames_dict.has_key(genename.lower()):
			return self.synnames_dict[genename.lower()]
		else:
			# print(genename+' does not exist in synnames dict')
			return ''

	def normalize_by_aa(self):  # normalize codon frequencies within each amino acid
		df=self.codonusage.copy()
		
		for name, group in df.groupby('amino_acid'):
			subdf=group.drop('amino_acid',axis='columns')
			df.ix[subdf.index,1:]=subdf.apply(lambda x: x/sum(x))

		return df.fillna(0)

	def normalize_by_totalaa(self):
		df=self.codonusage.copy()
		df.ix[:,1:]=df.ix[:,1:]/df.ix[:,1:].sum()
		
		return df.fillna(0)

	# use this after creating codonusage by add_genes
	def codon_expression(self,expressiondf):  
	# expressiondf is a pd.DataFrame and contains gene names (1st column) and counts (2nd column)
		if len(self.codonusage.columns)>2:   # check if genes have been added to codonusage
			df=expressiondf
			codon_expression=self.codonusage.copy()
			locustag_list=self.codonusage.columns[1:]
			expression_counts=pd.DataFrame(index=['count'],columns=locustag_list).fillna(0)
			genenotfound=[]
			
			for gene in df.index:  # loop over genes in the expression data
				locus_tag=self.lookup_locustag(gene)
				is_genein=locustag_list.isin([locus_tag])
				if any(is_genein):  # see if the 'gene' is listed in .expressed_codonusage
					expression_counts.ix['count',locus_tag]=df.ix[gene]
				else:
					locus_tag=self.lookup_synnames(gene)
					is_genein=locustag_list.isin([locus_tag])
					if any(is_genein): # see if the syn name of 'gene' is listed in .expressed_codonusage
						expression_counts.ix['count',locus_tag]=df.ix[gene]

					else:
						genenotfound.append(gene)
					
			print genenotfound
			codon_expression.ix[:,1]=self.codonusage.ix[:,1:]*expression_counts.ix['count',:]

			return codon_expression.sum(axis='columns')

		else:
			print('No gene has been added. Perform add_genes method first')




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



	def weights_from_tRNAgenecopy(self,gbrecord):  
	# gbrecord is SeqIO.seq object: gbrecord=SeqIO.parse(open('filepath'),'genbank').next()
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

		weights=self.codondecode_matrix.dot(self.tRNAtable['count'])
		weights=weights/weights.max()

		return weights

	def weights_from_tRNAabundance(self,trnadf): 
	# trnadf should be pd.DataFrame and include anticodon (1st column) and abundances (2nd column) 
		self.tRNAtable=pd.DataFrame(index=self.codondecode_matrix.columns,columns=['amino_acid']).fillna(0)
		self.tRNAtable=pd.merge(self.tRNAtable,trnadf,how='outer',left_index=True,right_index=True).reindex(self.codondecode_matrix.columns)

		weights=self.codondecode_matrix.dot(self.tRNAtable['ave'].fillna(0))
		weights=weights/weights.max()

		return weights
		
	
	def lookup_wobbleanticodon(self,codon):
		return self.codondecode_matrix.ix[codon,:][self.codondecode_matrix.ix[codon,:]>0]







	


if __name__ == "__main__":

	
	# load fasta
	

	# make name table
	# import ipdb; ipdb.set_trace()#

	




	import ipdb; ipdb.set_trace()#

	








