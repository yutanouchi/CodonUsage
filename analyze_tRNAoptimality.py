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
filepath='sequence_data/MG1655.gb'
record=SeqIO.parse(open(filepath,'rU'),'genbank').next()
book=codonanalyzer.CodonBook()
book.add_genes(record)

# load expression data
proteincounts=pd.read_csv('lit_data/ProteinCounts_Li_2014.csv',index_col='Gene')
demand=book.codon_expression(proteincounts['MOPS complete'])

# calculate tRNA supply
supply=codonanalyzer.tAI()
trna_abundance=pd.read_csv('lit_data/tRNA_Dong.csv',index_col=0)
weights=supply.weights_from_tRNAabundance(trna_abundance)

demand=demand/demand.max()
df=pd.merge(weights.to_frame(name='supply'),demand.to_frame(name='demand'),how='outer',left_index=True,right_index=True)
df.ix[:,'ratio']=df.ix[:,'supply']/df.ix[:,'demand']
df_sorted=df.ix[df['ratio'].sort_values().index,:]

fig, axes=plt.subplots(3,1, figsize](15,10))
df_sorted['ratio'].plot(kind='bar', ax=axes[0])
axes[0].set_xticklabels([])
df_sorted['supply'].plot(kind='bar', ax=axes[1])
axes[1].set_xticklabels([])
df_sorted['demand'].plot(kind='bar', ax=axes[2])
fig.show()




import ipdb; ipdb.set_trace()#


