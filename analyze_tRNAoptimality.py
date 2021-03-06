#!/usr/bin/python


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import codonanalyzer

plt.style.use('ggplot')

# load data
colifilepath='sequence_data/MG1655.gb'
colirecord=SeqIO.parse(open(colifilepath,'rU'),'genbank').next()
colibook=codonanalyzer.CodonBook()
colibook.add_genes(colirecord)


# load expression data
proteincounts=pd.read_csv('lit_data/ProteinCounts_Li_2014.csv',index_col='Gene')
demand=colibook.codon_expression(proteincounts['MOPS minimal'])
demand=demand/demand.max()
sorted_demand=demand.sort_values()

# calculate tRNA supply
coli_te=codonanalyzer.TransEff()
trna_abundance=pd.read_csv('lit_data/tRNA_Dong.csv',index_col=0)
supply=coli_te.weights_from_tRNAabundance(trna_abundance)
sorted_supply=supply.sort_values()

nte=supply/demand
nte=nte/nte.max()
nte['amino_acid']=nte.index.map(lambda x:colibook.codontable[x])
nte_sorted=nte.sort_values()

import ipdb; ipdb.set_trace()#
fig, axes=plt.subplots(3,1, figsize=(12,8))
nte_sorted.plot(kind='bar', ax=axes[0])
axes[0].set_xticklabels([])
supply[nte_sorted.index].plot(kind='bar', ax=axes[1])
axes[1].set_xticklabels([])
demand[nte_sorted.index].plot(kind='bar', ax=axes[2])


# fig, axes=plt.subplots(3,1, figsize=(12,8))
# nte_sorted.plot(kind='bar', ax=axes[0])
# axes[0].set_xticklabels([])
# supply[nte_sorted.index].plot(kind='bar', ax=axes[1])
# axes[1].set_xticklabels([])
# demand[nte_sorted.index].plot(kind='bar', ax=axes[2])
# fig.savefig('MG1655_tRNAoptimality_from_tRNAabundance.pdf')
# fig.show()




import ipdb; ipdb.set_trace()#


