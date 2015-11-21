#!/usr/bin/python

import re
import types
import abc
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
	__codontable={'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu','CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
				'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met','GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
				'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser','CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
				'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
				'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop','CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
				'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
				'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp','CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
				'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}
	__reversecodontable={'ALA':['GCT','GCC','GCA','GCG'],'ARG':['CGT','CGC','CGA','CGG','AGA','AGG'],'ASN':['AAT','AAC'],
					   'CYS':['TGT','TGC'],'GLN':['CAA','CAG'],'GLU':['GAA','GAG'],'GLY':['GGT','GGC','GGA','GGG'],
					   'HIS':['CAT','CAC'],'ILE':['ATT','ATC','ATA'],'MET':['ATG'],'LEU':['TTA','TTG','CTT','CTC','CTA','CTG'],
					   'LYS':['AAA','AAG'],'Phe':['TTT','TTC'],'PRO':['CCT','CCC','CCA','CCG'],'SER':['TCT','TCC','TCA','TCG','AGT','AGC'],
					   'THR':['ACT','ACC','ACA','ACG'],'TRP':['TGG'],'TYR':['TAT','TAC'],'VAL':['GTT','GTC','GTA','GTG'],
					'Stop':['TAA','TGA','TAG']}
	__aminoacidtable={'ALA':'A','A':'Ala','ARG':'R','R':'Arg','ASN':'N','N':'Asn','ASP':'D','D':'Asp','CYS':'C','C':'Cys',
					'GLN':'Q','Q':'Gln','GLU':'E','E':'Glu','GLY':'G','G':'Gly','HIS':'H','H':'His','ILE':'I','I':'Ile',
					'MET':'M','M':'Met','LEU':'L','L':'Leu','LYS':'K','K':'Lys','PHE':'F','F':'Phe','PRO':'P','P':'Pro',
					'SER':'S','S':'Ser','THR':'T','T':'Thr','TRP':'W','W':'Trp','TYR':'Y','Y':'Tyr','VAL':'V','V':'Val','STOP':'Stop'}
	
	def __init__(self):
		
		self.codonusage=pd.DataFrame([],index=self.__codontable.keys())
		self.codonusage['amino_acid']=self.codonusage.index.map(lambda x: self.__codontable[x])
		self.aalist=self.codonusage['amino_acid'].unique() # amino acid list
		self.aausage=pd.DataFrame([],index=self.aalist)
		self.locustag_dict={} # look-up dictionary between primary gene names and locus tag
		self.synnames_dict={}  # look-up dictionary from syn gene names to locus tag

	@classmethod
	def lookup_codon(cls,arg):
		
		if cls.__codontable.has_key(arg.upper()):
			return cls.__codontable[arg.upper()]
		elif cls.__reversecodontable.has_key(arg.upper()):
			return cls.__reversecodontable[arg.upper()]
		else:
			return ''

	def add_genes(self,seqrecord):  
	# seqrecord is a genbank SeqIO.seq object from Biopython. Genes are identified by locus tag (e.g. b0001)
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

				if f.qualifiers.has_key('protein_id'):
					protein_id=f.qualifiers['protein_id'][0]  # protein id
				else:
					protein_id=genename

				if f.qualifiers.has_key('gene_synonym'):
					synnames=re.split(';\ *',f.qualifiers['gene_synonym'][0])
				else:
					synnames=[]

				# update dict. 
				# Note that there are multiple copies of ins genes (e.g. 4 copies of insB1)
				# and genename-to-locus_tag dict will be overwritten and leavs only the last entry  
				self.locustag_dict.update({genename.lower():locus_tag})	
				self.locustag_dict.update({locus_tag.lower():genename})
				self.locustag_dict.update({protein_id.lower():locus_tag})
				self.synnames_dict.update({name.lower():locus_tag for name in synnames})

				# fetch cds sequence
				cds_seq=seq[start:end]
				if strand<0: # if cds is on negative strand
					cds_seq=cds_seq.reverse_complement()
				cds_seq=str(cds_seq)
		
				# extract codon and amino acid info
				triplets=pd.Series([cds_seq[i:i+3] for i in range(3,len(cds_seq)-3,3)])   # split sequences into codons (exclude start and stop codon)
				aa=triplets.apply(lambda x:self.__codontable[x])
				triplet_counts=triplets.value_counts().to_dict()
				aa_counts=aa.value_counts().to_dict()
				df=pd.DataFrame(triplet_counts.values(),index=triplet_counts.keys(),columns=[locus_tag])
				df2=pd.DataFrame(aa_counts.values(),index=aa_counts.keys(),columns=[locus_tag])

				# add to codonusage and aa dataframes
				self.codonusage=pd.merge(self.codonusage,df,how='outer',left_index=True,right_index=True).fillna(0)
				# self.codonusage=self.codonusage.sort_values('amino_acid')
				self.aausage=pd.merge(self.aausage,df2,how='outer',left_index=True,right_index=True).fillna(0)


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

	# calculate the codon usage in terms of gene expression based on expression data
	# use this after creating codonusage by add_genes
	def codon_expression(self,expressiondf):  
	# expressiondf is a pd.DataFrame and contains gene names (1st column) and counts (2nd column)
		if len(self.codonusage.columns)>2:   # check if genes have been added to codonusage
			codon_expression=self.codonusage.copy()
			locustag_list=self.codonusage.columns[1:]
			expression_counts=pd.DataFrame(index=['count'],columns=locustag_list).fillna(0)
			genenotfound=[]
			
			for gene in expressiondf.index:  # loop over genes in the expression data
				locus_tag=self.lookup_locustag(gene)
				is_genein=locustag_list.isin([locus_tag])
				if any(is_genein):  # see if the 'gene' is listed in codon_expression
					expression_counts.ix['count',locus_tag]=expressiondf.ix[gene]
				else:
					locus_tag=self.lookup_synnames(gene)
					is_genein=locustag_list.isin([locus_tag])
					if any(is_genein): # see if the syn name of 'gene' is listed in codon_expression
						expression_counts.ix['count',locus_tag]=expressiondf.ix[gene]

					else:
						genenotfound.append(gene)
					
			print genenotfound
			codon_expression.ix[:,1:]=self.codonusage.ix[:,1:].apply(lambda x:expression_counts.ix['count',:]*x,axis='columns')
			
			return codon_expression.sum(axis='columns')

		else:
			print('No gene has been added. Perform add_genes method first')




class GenomeWalker(object):

	def __init__(self,seqrecord):
	# seqrecord is a genbank SeqIO.seq object from Biopython. 
		self.seq=seqrecprd.seq
		self.genomesize=len(seq)

	def walk_and_find(self,bedobj,window=200000):
	# bedobj is a BedTool object that contains features of interest (e.g., tRNA)
		arm=np.round(window/2)

		for origin in range(0,self.genomesize):  # origin is zero-based
			left=origin-arm  # zero-based
			right=origin+1+arm  # one-based
			
			if left < 0:

			if right > genomesize:

			else:
				find_feature(left,right)


		def find_feature(self,left,right)




class TransEff(object):
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
	# gbrecord is SeqIO.seq genbank object that contains genomic info: gbrecord=SeqIO.parse(open('filepath'),'genbank').next()
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
		return self.codondecode_matrix.ix[codon.upper(),:][self.codondecode_matrix.ix[codon.upper(),:]>0]







	


if __name__ == "__main__":

	


	import ipdb; ipdb.set_trace()#

	








