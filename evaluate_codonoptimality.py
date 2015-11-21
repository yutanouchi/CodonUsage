#!/usr/bin/python


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import numpy as np
import codonanalyzer

plt.style.use('ggplot')


# load data
colifilepath='sequence_data/MG1655.gb'
colirecord=SeqIO.parse(open(colifilepath,'rU'),'genbank').next()
colibook=codonanalyzer.CodonBook()
colibook.add_genes(colirecord)

lambdafilepath='sequence_data/lambda.gb'
lambdarecord=SeqIO.parse(open(lambdafilepath,'rU'),'genbank').next()
lambdabook=codonanalyzer.CodonBook()
lambdabook.add_genes(lambdarecord)


# load expression data
proteincounts=pd.read_csv('lit_data/ProteinCounts_Li_2014.csv',index_col='Gene')
demand=colibook.codon_expression(proteincounts['MOPS complete'])
demand=demand/demand.max()

# calculate tRNA supply
coli_te=codonanalyzer.TransEff()
trna_abundance=pd.read_csv('lit_data/tRNA_Dong.csv',index_col=0)
supply=coli_te.weights_from_tRNAabundance(trna_abundance)

# translation efficiency index nte
nte=supply/demand
nte=nte/nte.max()

# phage codon usage


import ipdb; ipdb.set_trace()#




def randomize_codon(aaseq,dict):
	








